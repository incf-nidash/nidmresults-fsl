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
from constants import *
from To_NIDM import To_NIDM

''' Parse an FSL result directory to extract the pieces information stored in NIDM-Results
'''
class FSLparser(To_NIDM, object):

    def __init__(self, *args, **kwargs):
        
        self.feat_dir = kwargs.pop('feat_dir')
        
        self.export_dir = os.path.join(self.feat_dir, 'nidm')
        
        design_file = os.path.join(self.feat_dir, 'design.fsf');
        design_file_open = open(design_file, 'r')
        self.design_txt = design_file_open.read()
        
        report_file = os.path.join(self.feat_dir, 'report_poststats.html')
        self.reportParser = MyFSLReportParser();
        file = open(report_file, 'r')
        self.reportParser.feed(file.read());

        self.feat_version = self.reportParser.feat_version

        super(FSLparser, self).__init__(*args, **kwargs)

        # # self.nidm.create_thresholds(**reportParser.threshold)
        # self.nidm.create_software(reportParser.feat_version)

        # self.parse_feat_dir()

    # Main function: parse a feat directory and build the corresponding NIDM graph
    def parse_feat_dir(self):
        # self.add_report_file(os.path.join(self.feat_dir, 'report_poststats.html'))
        
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

    def find_threshold(self):
        return self.reportParser.threshold

    def find_residuals_file(self):
        residuals_file = os.path.join(self.feat_dir, 'stats', 'sigmasquareds.nii.gz')
        # FIXME: Check if there is an alternative file to use here
        if not os.path.isfile(residuals_file):
            residuals_file = None;
        return residuals_file

    def find_grand_mean_file(self):
        grand_mean_file = os.path.join(self.feat_dir, 'mean_func.nii.gz')

        # FIXME: Check if there is an alternative file to use here
        if not os.path.isfile(grand_mean_file):
            grand_mean_file = None;
        return grand_mean_file

    def find_design_matrix(self):
        design_matrix_file = open(os.path.join(self.feat_dir, 'design.mat'), 'r')
        design_matrix = np.loadtxt(design_matrix_file, skiprows=5, ndmin=2)
        return design_matrix

    def find_mask_file(self):
        mask_file = os.path.join(self.feat_dir, 'mask.nii.gz')
        return mask_file

    # Add model fitting, residuals map
    def find_noise_model(self):
        fmri_level_search = re.compile(r'.*set fmri\(level\) (?P<fmrilevel>\d+).*')
        fmri_level_found = fmri_level_search.search(self.design_txt)
        fmri_level = int(fmri_level_found.group('fmrilevel'))
        first_level = (fmri_level == 1)

        if first_level:
            variance_homo = True
            dependance = SERIALLY_CORR
            variance_spatial = SPATIALLY_LOCAL
            dependance_spatial = SPATIALLY_REGUL
        else:
            variance_homo = False
            dependance = INDEPEDENT_CORR
            variance_spatial = SPATIALLY_LOCAL
            dependance_spatial = None

        error_distribution = GAUSSIAN_DISTRIBUTION

        return dict(error_distribution=error_distribution,
            variance_homo=variance_homo, dependance=dependance,
            variance_spatial=variance_spatial, dependance_spatial=dependance_spatial)

    # For a parameter estimate, create the parameter estimate map emtity
    def add_parameter_estimate(self, pe_file, pe_num):
        self.nidm.create_parameter_estimate(pe_file, pe_num)

    # Find reference space 
    def find_coordinate_system(self):
        standard_space_search = re.compile(r'.*set fmri\(regstandard_yn\) (?P<isStandard>[\d]+).*')
        extracted_data = standard_space_search.search(self.design_txt) 
        standard_space = bool(extracted_data.group('isStandard'))

        if standard_space:
            standard_space_search = re.compile(r'.*set fmri\(alternateReference_yn\) (?P<isCustom>[\d]+).*')
            extracted_data = standard_space_search.search(self.design_txt) 
            if not extracted_data is None:
                custom_standard = (extracted_data.group('isCustom') == "1");
            else:
                standard_space_search = re.compile(r'.*set fmri\(regstandard\) (?P<regstd>.+).*')
                extracted_data = standard_space_search.search(self.design_txt) 
                if not extracted_data is None:
                    custom_standard = True;
        else:
            custom_standard = False;

        # As in https://github.com/ni-/ni-dm/issues/52 (not accepted yet)
        if not standard_space:
            coordinate_system = NIDM['SubjectSpace'];
        else:
            if not custom_standard:
                coordinate_system = NIDM['IcbmMni152NonLinear6thGenerationCoordinateSystem'];
            else:
                coordinate_system = NIDM['StandarizedSpace'];

        return coordinate_system


    # For a given contrast, create the contrast map, contrast variance map, contrast and statistical map emtities
    def add_contrast(self, contrast_num):
        contrast_file = os.path.join(self.feat_dir, 'stats', 'cope'+str(contrast_num)+'.nii.gz')
        varcontrast_file = os.path.join(self.feat_dir, 'stats', 'varcope'+str(contrast_num)+'.nii.gz')
        stat_map_file = os.path.join(self.feat_dir, 'stats', 'tstat'+str(contrast_num)+'.nii.gz')
        z_stat_map_file = os.path.join(self.feat_dir, 'stats', 'zstat'+str(contrast_num)+'.nii.gz')

        # Get contrast name and contrast weights from design.fsf file
        contrast_name_search = re.compile(r'.*set fmri\(conname_real\.'+contrast_num+'\) "(?P<contrastName>[\w\s><]+)".*')
        extracted_data = contrast_name_search.search(self.design_txt) 

        contrast_weight_search = re.compile(r'.*set fmri\(con_real'+contrast_num+'\.\d+\) (?P<contrastWeight>\d+)')
        contrastWeights = str(re.findall(contrast_weight_search, self.design_txt)).replace("'", '')

        # FIXME: to do only once (and not each time we load a new contrast)
        dofFile = open(os.path.join(self.feat_dir, 'stats', 'dof'), 'r')
        dof = float(dofFile.read())

        # Find connectivity criterion
        # FIXME: maybe not always "4"?
        feat_post_log_file = open(os.path.join(self.feat_dir, 'logs', 'feat4_post'), 'r')
        feat_post_log = feat_post_log_file.read()
        connectivity_search = re.compile(r'.* --connectivity=(?P<connectivity>\d+)+ .*')
        connectivity = int(connectivity_search.search(feat_post_log).group('connectivity'))

        peak_dist_search = re.compile(r'.* --peakdist=(?P<peakdist>\d+)+ .*')
        peak_dist_found = peak_dist_search.search(feat_post_log)
        if peak_dist_found:
            peak_dist = float(peak_dist_found.group('peakdist'))
        else:
            # If not specified, default value is zero (cf. http://fsl.fmrib.ox.ac.uk/fsl/fslwiki/Cluster)
            peak_dist = 0.0

        num_peak_search = re.compile(r'.* --num=(?P<numpeak>\d+)+ .*')
        num_peak_found = num_peak_search.search(feat_post_log)
        if num_peak_found:
            num_peak = int(num_peak_found.group('numpeak'))
        else:
            num_peak_search = re.compile(r'.* -n=(?P<numpeak>\d+)+ .*')
            num_peak_found = num_peak_search.search(feat_post_log)
            if num_peak_found:
                num_peak = int(num_peak_found.group('numpeak'))
            else:
                # If not specified, default value is inf? (cf. http://fsl.fmrib.ox.ac.uk/fsl/fslwiki/Cluster)
                # Is it ok to say no limit with -1 (as for Inf we would need float...)
                num_peak = -1


        # # FIXME deal with the case in which we are contrast masking by more than one contrast
        # contrast_masking_search = re.compile(r'.*set fmri\(conmask'+contrast_num+'_(?P<maskingconnum>\d+)\) (?P<domask>\d+).*')
        # contrast_masking_found = contrast_masking_search.search(self.design_txt)
        # do_contrast_masking = float(contrast_masking_found.group('domask'))
        # if do_contrast_masking:
        #     contrast_masking_num = contrast_masking_found.group('maskingconnum')
        #     contrast_masking_file = 
        # else:
        #     contrast_masking_num = None
        # FIXME: We need an example with more than one contrast to code contrast masking        
        contrast_masking_file = self.maskFile

        self.nidm.create_contrast_map(contrast_file, varcontrast_file, stat_map_file, z_stat_map_file,
            extracted_data.group('contrastName').strip(), contrast_num, dof, contrastWeights, 
            connectivity, peak_dist, num_peak, contrast_masking_file)

    # Create the search space entity generated by an inference activity
    def add_search_space(self):
        # FIXME this needs to be estimated
        search_space_file = os.path.join(self.feat_dir, 'mask.nii.gz')
        smoothnessFile = os.path.join(self.feat_dir, 'stats', 'smoothness')

        # Load DLH, VOLUME and RESELS
        smoothness = np.loadtxt(smoothnessFile, usecols=[1])
        self.nidm.create_search_space(search_space_file=search_space_file, search_volume=int(smoothness[1]), resel_size_in_voxels=float(smoothness[2]), dlh=float(smoothness[0]))


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

                # Though peak coordinates in voxels are integer, we use a float type to comply with the rdfs:range
                self.nidm.create_peak(id=int(peakIndex), x=float(peak_row[2]), y=float(peak_row[3]), z=float(peak_row[4]), 
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
        self.threshold = dict()

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

            voxel_threshold = None

            extracted_data = patternVoxelThresh.search(data) 
            if extracted_data is not None:
                self.feat_version = extracted_data.group('featversion')
                voxel_threshold = None;
                voxel_p_corr = float(extracted_data.group('pvalue'))
                voxel_p_uncorr = None
                extent_value = 0;
                extent_p_corr = 1
                extent_p_uncorr = 1
                # self.threshType = extracted_data.group('threshtype')
                self.found_intro = True;
            else:
                patternClusterThresh = re.compile(r'.*Version (?P<featversion>\d+\.\d+),.* thresholded using (?P<threshtype>.*) determined by Z\>(?P<zvalue>\d+\.\d+) and a .* P=(?P<pvalue>\d+\.\d+) .*')
                extracted_data = patternClusterThresh.search(data) 

                if extracted_data is not None:
                    self.feat_version = extracted_data.group('featversion')
                    voxel_threshold = float(extracted_data.group('zvalue'))
                    voxel_p_corr = None
                    voxel_p_uncorr = None
                    extent_value = None;
                    extent_p_corr = float(extracted_data.group('pvalue'));
                    extent_p_uncorr = None
                    # self.threshType = extracted_data.group('threshtype')
                    self.found_intro = True;

            if voxel_threshold:
                self.threshold = dict(voxel_threshold=voxel_threshold,
                    voxel_p_corr=voxel_p_corr, voxel_p_uncorr=voxel_p_uncorr, 
                    extent=extent_value, extent_p_corr=extent_p_corr,
                    extent_p_uncorr=extent_p_uncorr)



