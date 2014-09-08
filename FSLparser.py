'''Export of FSL results into NI-DM

@author: Camille Maumet <c.m.j.maumet@warwick.ac.uk>
@copyright: University of Warwick 2013-2014
'''

from HTMLParser import HTMLParser
from htmlentitydefs import name2codepoint
import re
from prov.model import ProvBundle, ProvRecord, ProvEntity
# import prov.model.graph
import os
import numpy as np
import nibabel as nib
from NIDMStat import NIDMStat


''' Parse an FSL result directory to extract the pieces information stored in NI-DM (for statistical results)
'''
class FSL_NIDM():

    def __init__(self, *args, **kwargs):
        self.feat_dir = kwargs.pop('feat_dir')
        self.export_dir = os.path.join(self.feat_dir, 'nidm')
        
        self.design_file = os.path.join(self.feat_dir, 'design.fsf');
        self.find_reference_space();

        self.nidm = NIDMStat(export_dir=self.export_dir, standard_space=self.standard_space, custom_standard=self.custom_standard);

        
        self.parse_feat_dir()

    # Main function: parse a feat directory and build the corresponding NI-DM graph
    def parse_feat_dir(self):
        self.add_report_file(os.path.join(self.feat_dir, 'report_poststats.html'))
        self.add_model_fitting()
        self.maskFile = os.path.join(self.feat_dir, 'mask.nii.gz')
        self.add_search_space()

        # Find parameter estimates
        for file in os.listdir(os.path.join(self.feat_dir, 'stats')):
            if file.startswith("pe"):
                if file.endswith(".nii.gz"):
                    s = re.compile('pe\d+')
                    penum = s.search(file)
                    penum = penum.group()
                    penum = penum.replace('pe', '')
                    self.add_parameter_estimate(os.path.join(self.feat_dir, 'stats', file), penum)

        # Find excursion sets (in a given feat directory we have one excursion set per contrast)
        for file in os.listdir(self.feat_dir):
            if file.startswith("thresh_zstat"):
                if file.endswith(".nii.gz"):
                    s = re.compile('zstat\d+')
                    zstatnum = s.search(file)
                    zstatnum = zstatnum.group()
                    statnum = zstatnum.replace('zstat', '')
                    self.add_contrast(statnum)
                    self.add_clusters_peaks(statnum)        

    # Add model fitting, residuals map
    def add_model_fitting(self):
        residuals_file = os.path.join(self.feat_dir, 'stats', 'sigmasquareds.nii.gz')
        # FIXME: Check if there is an alternative file to use here
        if not os.path.isfile(residuals_file):
            residuals_file = None;

        grand_mean_file = os.path.join(self.feat_dir, 'mean_func.nii.gz')

        # FIXME: Check if there is an alternative file to use here
        if not os.path.isfile(grand_mean_file):
            grand_mean_file = None;

        design_matrix_file = open(os.path.join(self.feat_dir, 'design.mat'), 'r')
        design_matrix = np.loadtxt(design_matrix_file, skiprows=5, ndmin=2)

        self.nidm.create_model_fitting(residuals_file, grand_mean_file, design_matrix)

    # For a parameter estimate, create the parameter estimate map emtity
    def add_parameter_estimate(self, pe_file, pe_num):
        self.nidm.create_parameter_estimate(pe_file, pe_num)

    # Find reference space 
    def find_reference_space(self):
        designFile = open(self.design_file, 'r')
        designTxt = designFile.read()
        standard_space_search = re.compile(r'.*set fmri\(regstandard_yn\) (?P<isStandard>[\d]+).*')
        extracted_data = standard_space_search.search(designTxt) 
        self.standard_space = bool(extracted_data.group('isStandard'))

        if self.standard_space:
            standard_space_search = re.compile(r'.*set fmri\(alternateReference_yn\) (?P<isCustom>[\d]+).*')
            extracted_data = standard_space_search.search(designTxt) 
            if not extracted_data is None:
                self.custom_standard = (extracted_data.group('isCustom') == "1");
            else:
                standard_space_search = re.compile(r'.*set fmri\(regstandard\) (?P<regstd>.+).*')
                extracted_data = standard_space_search.search(designTxt) 
                if not extracted_data is None:
                    self.custom_standard = True;
        else:
            self.custom_standard = False;


    # For a given contrast, create the contrast map, contrast variance map, contrast and statistical map emtities
    def add_contrast(self, contrast_num):
        contrast_file = os.path.join(self.feat_dir, 'stats', 'cope'+str(contrast_num)+'.nii.gz')
        varcontrast_file = os.path.join(self.feat_dir, 'stats', 'varcope'+str(contrast_num)+'.nii.gz')
        stat_map_file = os.path.join(self.feat_dir, 'stats', 'tstat'+str(contrast_num)+'.nii.gz')
        z_stat_map_file = os.path.join(self.feat_dir, 'stats', 'zstat'+str(contrast_num)+'.nii.gz')

        # Get contrast name and contrast weights from design.fsf file
        designFile = open(self.design_file, 'r')
        designTxt = designFile.read()
        contrast_name_search = re.compile(r'.*set fmri\(conname_real\.'+contrast_num+'\) "(?P<contrastName>[\w\s><]+)".*')
        extracted_data = contrast_name_search.search(designTxt) 

        contrast_weight_search = re.compile(r'.*set fmri\(con_real'+contrast_num+'\.\d+\) (?P<contrastWeight>\d+)')
        contrastWeights = str(re.findall(contrast_weight_search, designTxt)).replace("'", '')

        # FIXME: to do only once (and not each time we load a new contrast)
        dofFile = open(os.path.join(self.feat_dir, 'stats', 'dof'), 'r')
        dof = float(dofFile.read())

        self.nidm.create_contrast_map(contrast_file, varcontrast_file, stat_map_file, z_stat_map_file,
            extracted_data.group('contrastName').strip(), contrast_num, dof, contrastWeights)

    # Create the search space entity generated by an inference activity
    def add_search_space(self):
        search_space_file = os.path.join(self.feat_dir, 'mask.nii.gz')
        smoothnessFile = os.path.join(self.feat_dir, 'stats', 'smoothness')

        # Load DLH, VOLUME and RESELS
        smoothness = np.loadtxt(smoothnessFile, usecols=[1])
        self.nidm.create_search_space(search_space_file=search_space_file, search_volume=int(smoothness[1]), resel_size_in_voxels=float(smoothness[2]), dlh=float(smoothness[0]))

    # Create the thresholding information for an inference activity (height threshold and extent threshold)
    def add_report_file(self, report_file):
        self.reportFile = report_file
        parser = MyFSLReportParser();
        file = open(report_file, 'r')
        parser.feed(file.read());

        self.nidm.create_thresholds( voxel_threshold=parser.get_voxel_thresh_value(), 
            voxel_p_uncorr=parser.get_voxel_p_uncorr(), 
            voxel_p_corr=parser.get_voxel_p_corr(), 
            extent=parser.get_extent_value(),
            extent_p_uncorr=parser.get_extent_p_uncorr(), 
            extent_p_corr=parser.get_extent_p_corr())

        self.nidm.create_software(parser.get_feat_version())

    # Create excursion set, clusters and peaks entities
    def add_clusters_peaks(self, stat_num):
        visualisation = os.path.join(self.feat_dir, 'rendered_thresh_zstat'+stat_num+'.png')

        # Excursion set
        zFileImg = os.path.join(self.feat_dir, 'thresh_zstat'+stat_num+'.nii.gz')
        self.nidm.create_excursion_set(excursion_set_file=zFileImg, stat_num=stat_num, visualisation=visualisation)

        # Cluster list (positions in voxels)
        cluster_file = os.path.join(self.feat_dir, 'cluster_zstat'+stat_num+'.txt')
        if not os.path.isfile(cluster_file):
            cluster_file = None;
        else:
            cluster_table = np.loadtxt(cluster_file, skiprows=1, ndmin=2)


        # Cluster list (positions in mm)
        cluster_std_file = os.path.join(self.feat_dir, 'cluster_zstat'+stat_num+'_std.txt')
        if not os.path.isfile(cluster_std_file):
            cluster_std_file = None;
            # cluster_std_table = np.zeros_like(cluster_table)*float('nan')
        else:
            cluster_std_table = np.loadtxt(cluster_std_file, skiprows=1, ndmin=2)
        
        # if not cluster_file:
            # cluster_table = np.zeros_like(cluster_std_table)*float('nan')            

        
        
        # Peaks
        peak_file = os.path.join(self.feat_dir, 'lmax_zstat'+stat_num+'.txt')
        if not os.path.isfile(peak_file):
            peak_file = None;
        else:
            peak_table = np.loadtxt(peak_file, skiprows=1, ndmin=2)

        peak_std_file = os.path.join(self.feat_dir, 'lmax_zstat'+stat_num+'_std.txt')
        if not os.path.isfile(peak_std_file):
            peak_std_file = None;
            # peak_std_table = np.zeros_like(peak_table)*float('nan')       
        else:
            peak_std_table = np.loadtxt(peak_std_file, skiprows=1, ndmin=2)       
        
        if (cluster_file is not None) and (cluster_std_file is not None):
            clusters_join_table = np.column_stack((cluster_table, cluster_std_table))
            for cluster_row in clusters_join_table: 
                self.nidm.create_cluster(id=int(cluster_row[0]), size=int(cluster_row[1]), pFWER=float(cluster_row[2]),
                    x=float(cluster_row[8]),y=float(cluster_row[9]),z=float(cluster_row[10]),
                    x_std=float(cluster_row[24]),y_std=float(cluster_row[25]),z_std=float(cluster_row[26]),
                    stat_num=stat_num)
        elif (cluster_file is not None):
            for cluster_row in cluster_table: 
                self.nidm.create_cluster(id=int(cluster_row[0]), size=int(cluster_row[1]), pFWER=float(cluster_row[2]),
                    x=float(cluster_row[8]),y=float(cluster_row[9]),z=float(cluster_row[10]), stat_num=stat_num)
        elif (cluster_std_file is not None):
            for cluster_row in cluster_std_table: 
                self.nidm.create_cluster(id=int(cluster_row[0]), size=int(cluster_row[1]), pFWER=float(cluster_row[2]),
                    x_std=float(cluster_row[8]),y_std=float(cluster_row[9]),z_std=float(cluster_row[10]), stat_num=stat_num)
        
        prev_cluster = -1
        if (peak_file is not None) and (peak_std_file is not None):
            peaks_join_table = np.column_stack((peak_table, peak_std_table))
            for peak_row in peaks_join_table:    
                cluster_id = int(peak_row[0])

                if not cluster_id == prev_cluster:
                    peakIndex = 1;

                self.nidm.create_peak(id=int(peakIndex), x=int(peak_row[2]), y=int(peak_row[3]), z=int(peak_row[4]), 
                    x_std=float(peak_row[7]), y_std=float(peak_row[8]), z_std=float(peak_row[9]),
                    equivZ=float(peak_row[1]), cluster_id=cluster_id, stat_num=stat_num, max_peak=(peakIndex==1))
                prev_cluster = cluster_id

                peakIndex = peakIndex + 1
        elif (peak_file is not None):
            for peak_row in peak_table:    
                cluster_id = int(peak_row[0])

                if not cluster_id == prev_cluster:
                    peakIndex = 1;

                self.nidm.create_peak(id=int(peakIndex), x=int(peak_row[2]), y=int(peak_row[3]), z=int(peak_row[4]),equivZ=float(peak_row[1]), cluster_id=cluster_id, stat_num=stat_num, max_peak=(peakIndex==1))
                prev_cluster = cluster_id

                peakIndex = peakIndex + 1
        elif (peak_std_file is not None):
            for peak_row in peak_std_table:    
                cluster_id = int(peak_row[0])

                if not cluster_id == prev_cluster:
                    peakIndex = 1;

                self.nidm.create_peak(id=int(peakIndex), x_std=int(peak_row[2]), y_std=int(peak_row[3]), z_std=int(peak_row[4]), equivZ=float(peak_row[1]), cluster_id=cluster_id, stat_num=stat_num, max_peak=(peakIndex==1))
                prev_cluster = cluster_id

                peakIndex = peakIndex + 1

        
    # Create a graph as a provn and a json serialisations
    def save_prov_to_files(self):
        self.nidm.save_prov_to_files()

'''HTML parser for FSL report files: extract the thresholding information

'''
# TODO: check if the thresholding information is stored elsewhere in FSL files
class MyFSLReportParser(HTMLParser):

    def __init__(self, *args, **kwargs):
        HTMLParser.__init__(self, *args, **kwargs)
        self.descriptions = []
        self.inside_a_element = 0
        self.hyperlinks = []
        self.found_intro = False;
        self.feat_version = ''
        self.pValue = []
        self.threshType = ''

    def handle_starttag(self, tag, attrs):
        if tag == "a":
            for name, value in attrs:
                if name == "href":
                    self.hyperlinks.append(value)
                    self.inside_a_element = 1

    def handle_endtag(self, tag):
        if tag == 'a':
            self.inside_a_element = 0
    def handle_data(self, data):
        if self.inside_a_element:
            self.descriptions.append(data)
        elif not self.found_intro:
            # Look for p-value, type of thresholding and feat version in introductory text
            patternVoxelThresh = re.compile(r'.*Version (?P<featversion>\d+\.\d+),.* thresholded using (?P<threshtype>.*) thresholding .* P=(?P<pvalue>\d+\.\d+)')

            extracted_data = patternVoxelThresh.search(data) 
            
            if extracted_data is not None:
                self.feat_version = extracted_data.group('featversion')
                self.voxel_thresh_value = None;
                self.voxel_p_corr = float(extracted_data.group('pvalue'))
                self.voxel_p_uncorr = None
                self.extent_value = 0;
                self.extent_p_corr = 1
                self.extent_p_uncorr = 1
                # self.threshType = extracted_data.group('threshtype')
                self.found_intro = True;
            else:
                patternClusterThresh = re.compile(r'.*Version (?P<featversion>\d+\.\d+),.* thresholded using (?P<threshtype>.*) determined by Z\>(?P<zvalue>\d+\.\d+) and a .* P=(?P<pvalue>\d+\.\d+) .*')
                extracted_data = patternClusterThresh.search(data) 

                if extracted_data is not None:
                    self.feat_version = extracted_data.group('featversion')
                    self.voxel_thresh_value = float(extracted_data.group('zvalue'))
                    self.voxel_p_corr = None
                    self.voxel_p_uncorr = None
                    self.extent_value = None;
                    self.extent_p_corr = float(extracted_data.group('pvalue'));
                    self.extent_p_uncorr = None
                    # self.threshType = extracted_data.group('threshtype')
                    self.found_intro = True;
    def get_feat_version(self):
        return self.feat_version

    def get_voxel_thresh_value(self):
        return self.voxel_thresh_value

    def get_voxel_p_corr(self):
        return self.voxel_p_corr

    def get_voxel_p_uncorr(self):
        return self.voxel_p_uncorr

    def get_extent_value(self):
        return self.extent_value

    def get_extent_p_corr(self):
        return self.extent_p_corr

    def get_extent_p_uncorr(self):
        return self.extent_p_uncorr


