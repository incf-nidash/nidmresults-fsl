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
from nidm_exporter.nidm_exporter import NIDMExporter
from nidm_exporter.objects.constants import *
from nidm_exporter.objects.modelfitting import *
from nidm_exporter.objects.contrast import *
from nidm_exporter.objects.inference import *
from fslobjects import *

''' Parse an FSL result directory to extract the pieces information stored in NIDM-Results
'''
class FSLtoNIDMExporter(NIDMExporter, object):

    def __init__(self, *args, **kwargs):
        self.feat_dir = kwargs.pop('feat_dir')        
        self.export_dir = os.path.join(self.feat_dir, 'nidm')
        self.design_file = os.path.join(self.feat_dir, 'design.fsf');
        self.report_file = os.path.join(self.feat_dir, 'report_poststats.html')
        self.feat_post_log_file = open(os.path.join(self.feat_dir, 'logs', 'feat4_post'), 'r')
        self.elements_to_export = list()
        self.activities_to_export = list()
        self.version = kwargs.pop('version')   
        # self.software = "FSL"
        self.coordinate_space_id = 1
        self.coordinate_system = None
        self.contrast_names_by_num = dict()

    def parse(self):        
        design_file_open = open(self.design_file, 'r')
        self.design_txt = design_file_open.read()
        
        self.reportParser = MyFSLReportParser();
        file = open(self.report_file, 'r')
        self.reportParser.feed(file.read());

        self.feat_post_log = self.feat_post_log_file.read()

        
        # self.nidm = NIDMStat(software="FSL", version="0.2.0", results_artefacts=self.feat_dir, export_dir=self.export_dir)

        if not self.coordinate_system:
            self._get_coordinate_system()

        super(FSLtoNIDMExporter, self).parse()
        

        # # self.nidm.create_thresholds(**reportParser.threshold)
        # self.nidm.create_software(reportParser.feat_version)

        # self.parse_feat_dir()

    def find_software(self):
        software = Software(feat_version = self.reportParser.feat_version)

        return software

    def find_model_fitting(self):
        design_matrix = self._get_design_matrix()
        data = self._get_data()
        error_model = self._get_error_model()

        rms_map = self._get_residual_mean_squares_map()
        param_estimates = self._get_parameter_estimate_maps()
        mask_map = self._get_mask_map()
        grand_mean_map = self._get_grand_mean(mask_map.file)

        activity = self._get_model_parameters_estimations(error_model)

        model_fitting = ModelFitting(activity, design_matrix, data, error_model, param_estimates, 
            rms_map, mask_map, grand_mean_map)

        return list([model_fitting])


    def _get_design_matrix(self):
        design_matrix_file = open(os.path.join(self.feat_dir, 'design.mat'), 'r')
        design_matrix_values = np.loadtxt(design_matrix_file, skiprows=5, ndmin=2)

        # In FSL, we have a single matrix per feat directory so a 
        # single ModelParamEstimation activity
        design_matrix = DesignMatrix(design_matrix_values, self.export_dir)

        return design_matrix

    def _get_data(self):
        data = Data(self.export_dir)
        # In FSL, we have a data entity per feat directory
        return data

    def _get_error_model(self):
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

        error_model = ErrorModel(error_distribution, variance_homo, 
            variance_spatial, dependance, dependance_spatial)

        # In FSL we have a single design per feat directory
        return error_model

    def _get_parameter_estimate_maps(self):
        # Find parameter estimates
        parameter_estimates = list()
        for filename in os.listdir(os.path.join(self.feat_dir, 'stats')):
            if filename.startswith("pe"):
                if filename.endswith(".nii.gz"):
                    s = re.compile('pe\d+')
                    penum = s.search(filename)
                    penum = penum.group()
                    penum = penum.replace('pe', '')
                    full_path_file = os.path.join(self.feat_dir, 'stats', filename)
                    parameter_estimate = ParameterEstimateMap(full_path_file, penum, 
                        self.coordinate_space_id, self.coordinate_system)
                    self.coordinate_space_id = self.coordinate_space_id + 1
                    # self.add_parameter_estimate(os.path.join(self.feat_dir, 'stats', filename), penum)
                    parameter_estimates.append(parameter_estimate)

        return parameter_estimates

    def _get_grand_mean(self, mask_file):
        grand_mean_file = os.path.join(self.feat_dir, 'mean_func.nii.gz')

        # FIXME: Check if there is an alternative file to use here
        if not os.path.isfile(grand_mean_file):
            grand_mean = None;
        else:
            grand_mean = GrandMeanMap(grand_mean_file, mask_file, self.coordinate_system, 
                self.coordinate_space_id, self.export_dir)
        
        return grand_mean

    def _get_coordinate_system(self):
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

        self.coordinate_system = coordinate_system

    def _get_residual_mean_squares_map(self):
        residuals_file = os.path.join(self.feat_dir, 'stats', 'sigmasquareds.nii.gz')
        # FIXME: Check if there is an alternative file to use here
        if not os.path.isfile(residuals_file):
            residuals = None;

        rms_map = ResidualMeanSquares(self.export_dir, residuals_file, 
            self.coordinate_system, self.coordinate_space_id)
        self.coordinate_space_id = self.coordinate_space_id + 1

        return rms_map

    def _get_mask_map(self):
        mask_file = os.path.join(self.feat_dir, 'mask.nii.gz')
        mask_map = MaskMap(self.export_dir, mask_file,
            self.coordinate_system, self.coordinate_space_id)
        self.coordinate_space_id = self.coordinate_space_id + 1
        return mask_map       


    # def find_contrast_weights(self):
        # con_weights = ContrastWeights(contrast_num, contrast_name, contrast_weights)

    # Retreive estimated contrasts. Return a list of Contrast objects.
    def find_contrasts(self):
        # In FSL there is a single model fitting activity per feat directory
        model_fitting_id = self.model_fittings[0].activity.id

        dofFile = open(os.path.join(self.feat_dir, 'stats', 'dof'), 'r')
        dof = float(dofFile.read())

        contrasts = dict()
        for filename in os.listdir(self.feat_dir):
            if filename.startswith("thresh_zstat"):
                if filename.endswith(".nii.gz"):
                    s = re.compile('zstat\d+')
                    zstatnum = s.search(filename)
                    zstatnum = zstatnum.group()
                    statnum = zstatnum.replace('zstat', '')

                    contrast_num = statnum

                    # Get contrast name
                    name_search = re.compile(r'.*set fmri\(conname_real\.'+statnum+\
                                                        '\) "(?P<contrastName>[\w\s><]+)".*')
                    extracted_data = name_search.search(self.design_txt) 
                    contrast_name = extracted_data.group('contrastName').strip()
                    self.contrast_names_by_num[statnum] = contrast_name

                    # Contrast estimation activity
                    estimation = ContrastEstimation(contrast_num, contrast_name, self.software.id)

                    # Get contrast weights
                    weight_search = re.compile(r'.*set fmri\(con_real'+statnum+\
                        '\.\d+\) (?P<contrastWeight>\d+)')
                    contrast_weights = str(re.findall(weight_search, self.design_txt)).replace("'", '')
                    weights = ContrastWeights(contrast_num, contrast_name, contrast_weights)
                   
                    # Find which betas were used to compute the contrast
                    pe_ids = list()
                    pe_index = 1
                    contrast_weights = contrast_weights.replace(' ', '').replace('[', '').replace(']', '').split(',')
                    for beta_index in contrast_weights:
                        if int(beta_index) == 1:
                            for model_fitting in self.model_fittings:
                                for pe in model_fitting.param_estimates:
                                    s = re.compile('pe\d+')
                                    pe_num = s.search(pe.file)
                                    pe_num = pe_num.group()
                                    pe_num = pe_num.replace('pe', '')
                                    if pe_num == pe_index:
                                        pe_ids.append(pe.id)
                            # self.provBundle.used(NIIRI['contrast_estimation_id_'+contrast_num], NIIRI['beta_map_id_'+str(peIndex)])
                        pe_index += 1;

                    # Convert to immutable tuple to be used as key
                    pe_ids = tuple(pe_ids)

                    # Get contrast, standard error and statistic files
                    con_file = os.path.join(self.feat_dir, 'stats', 'cope'+str(statnum)+'.nii.gz')
                    contrast_map = ContrastMap(con_file, contrast_num, contrast_name, self.coordinate_system, 
                        self.coordinate_space_id, self.export_dir)
                    self.coordinate_space_id += 1

                    varcontrast_file = os.path.join(self.feat_dir, 'stats', 'varcope'+str(statnum)+'.nii.gz')# In FSL varcope are exported
                    is_variance = True
                    std_err_map = ContrastStdErrMap(contrast_num, varcontrast_file, is_variance, 
                        self.coordinate_system, 
                        self.coordinate_space_id, self.export_dir)
                    self.coordinate_space_id += 1

                    stat_file = os.path.join(self.feat_dir, 'stats', 'tstat'+str(statnum)+'.nii.gz')
                    stat_map = StatisticMap(stat_file, 'T', contrast_num, contrast_name, dof, 
                        self.coordinate_system, self.coordinate_space_id, self.export_dir)
                    self.coordinate_space_id += 1

                    z_stat_file = os.path.join(self.feat_dir, 'stats', 'zstat'+str(statnum)+'.nii.gz')
                    z_stat_map = StatisticMap(z_stat_file, 'Z', contrast_num, contrast_name, dof, 
                        self.coordinate_system, self.coordinate_space_id, self.export_dir)
                    self.coordinate_space_id += 1

                    con = Contrast(contrast_num, contrast_name, weights, estimation, 
                        contrast_map, std_err_map, stat_map, z_stat_map)

                    contrasts.setdefault((model_fitting_id, pe_ids), list()).append(con)

        return contrasts

    def _search_in_fsf(self, regexp):
        info_search = re.compile(regexp)
        info_found = info_search.search(self.design_txt)
        info = info_found.group('info')
        return info

    # Retreive inference processes along with peaks and clusters. Return a list of Inference objects.
    def find_inferences(self):
        inferences = dict()
        # Find excursion sets (in a given feat directory we have one excursion set per contrast)
        for filename in os.listdir(self.feat_dir):
            if filename.startswith("thresh_zstat"):
                if filename.endswith(".nii.gz"):
                    s = re.compile('zstat\d+')
                    zstatnum = s.search(filename)
                    zstatnum = zstatnum.group()
                    stat_num = zstatnum.replace('zstat', '')

                    # Find corresponding contrast estimation activity
                    for contrasts in self.contrasts.values():
                        for contrast in contrasts:
                            s = re.compile('cope\d+')
                            con_num = s.search(contrast.contrast_map.file)
                            con_num = con_num.group()
                            con_num = con_num.replace('cope', '')
                            if con_num == stat_num:
                                contrast_id = contrast.estimation.id

                    # Inference activity
                    inference_act = InferenceActivity(stat_num, self.contrast_names_by_num[stat_num])

                    # Excursion set
                    visualisation = os.path.join(self.feat_dir, 'rendered_thresh_zstat'+stat_num+'.png')
                    zFileImg = os.path.join(self.feat_dir, 'thresh_zstat'+stat_num+'.nii.gz')
                    exc_set = ExcursionSet(zFileImg, stat_num, visualisation, self.coordinate_system, 
                        self.coordinate_space_id, self.export_dir)

                    # Height Threshold
                    prob_thresh = float(self._search_in_fsf(r'.*set fmri\(prob_thresh\) (?P<info>\d+\.?\d+).*'))
                    z_thresh = float(self._search_in_fsf(r'.*set fmri\(z_thresh\) (?P<info>\d+\.?\d+).*'))
                    thresh_type = self._search_in_fsf(r'.*set fmri\(thresh\) (?P<info>\d+).*')

                    # FIXME: deal with 0 = no thresh?
                    voxel_uncorr = (thresh_type == 1)
                    voxel_corr = (thresh_type == 2)
                    cluster_thresh = (thresh_type == 3)
                    
                    stat_threshold = None
                    p_corr_threshold = None
                    p_uncorr_threshold = None
                    if voxel_uncorr:
                        p_uncorr_threshold = prob_thresh
                    elif voxel_corr:
                        p_corr_threshold = prob_thresh
                    else:
                        stat_threshold = z_thresh
                        extent_p_corr = prob_thresh

                    height_thresh = HeightThreshold(stat_threshold, p_corr_threshold, p_uncorr_threshold)

                    # Extent Threshold
                    extent_thresh = ExtentThreshold(p_corr=extent_p_corr)

                    # Clusters (and associated peaks)
                    clusters = self._get_clusters_peaks(stat_num)

                    # Peak and Cluster Definition Criteria
                    peak_criteria = PeakCriteria(stat_num, self._get_num_peaks(), self._get_peak_dist())
                    clus_criteria = ClusterCriteria(stat_num, self._get_connectivity())

                    # Display mask
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
                    contrast_masking_file = self.find_mask_file()
                    display_mask = DisplayMaskMap(stat_num, contrast_masking_file, self.coordinate_system, self.coordinate_space_id, self.export_dir)
                    self.coordinate_space_id + 1


                    # Search space
                    search_space = self._get_search_space()

                    inference = Inference(inference_act, height_thresh, extent_thresh, 
                        peak_criteria, clus_criteria, display_mask, exc_set, clusters, 
                        search_space, self.software.id)

                    if contrast_id in inferences:
                        inferences[contrast_id].append(inference)
                    else:
                        inferences[contrast_id] = list([inference])

        return inferences


    # def find_excursion_sets(self):
    #     for filename in os.listdir(self.feat_dir):
    #         if filename.startswith("thresh_zstat"):
    #             if filename.endswith(".nii.gz"):
    #                 s = re.compile('zstat\d+')
    #                 zstatnum = s.search(filename)
    #                 zstatnum = zstatnum.group()
    #                 stat_num = zstatnum.replace('zstat', '')

                    # self.nidm.create_excursion_set(excursion_set_file=zFileImg, stat_num=stat_num, visualisation=visualisation)


    # Main function: parse a feat directory and build the corresponding NIDM graph
    def parse_feat_dir(self):
        # self.add_report_file(os.path.join(self.feat_dir, 'report_poststats.html'))
        
        
        self.add_search_space()




    def find_mask_file(self):
        mask_file = os.path.join(self.feat_dir, 'mask.nii.gz')
        return mask_file



    # For a parameter estimate, create the parameter estimate map emtity
    def add_parameter_estimate(self, pe_file, pe_num):
        self.nidm.create_parameter_estimate(pe_file, pe_num)


    def _get_num_peaks(self):
        num_peak_search = re.compile(r'.* --num=(?P<numpeak>\d+)+ .*')
        num_peak_found = num_peak_search.search(self.feat_post_log)
        if num_peak_found:
            num_peak = int(num_peak_found.group('numpeak'))
        else:
            num_peak_search = re.compile(r'.* -n=(?P<numpeak>\d+)+ .*')
            num_peak_found = num_peak_search.search(self.feat_post_log)
            if num_peak_found:
                num_peak = int(num_peak_found.group('numpeak'))
            else:
                # If not specified, default value is inf? (cf. http://fsl.fmrib.ox.ac.uk/fsl/fslwiki/Cluster)
                # Is it ok to say no limit with -1 (as for Inf we would need float...)
                num_peak = -1
        return num_peak

    def _get_peak_dist(self):
        peak_dist_search = re.compile(r'.* --peakdist=(?P<peakdist>\d+)+ .*')
        peak_dist_found = peak_dist_search.search(self.feat_post_log)
        if peak_dist_found:
            peak_dist = float(peak_dist_found.group('peakdist'))
        else:
            # If not specified, default value is zero (cf. http://fsl.fmrib.ox.ac.uk/fsl/fslwiki/Cluster)
            peak_dist = 0.0

        return peak_dist

    def _get_connectivity(self):
        # Find connectivity criterion
        # FIXME: maybe not always "4"?
        connectivity_search = re.compile(r'.* --connectivity=(?P<connectivity>\d+)+ .*')
        connectivity = int(connectivity_search.search(self.feat_post_log).group('connectivity'))

        return connectivity

    # Create the search space entity generated by an inference activity
    def _get_search_space(self):
        # FIXME this needs to be estimated
        search_space_file = os.path.join(self.feat_dir, 'mask.nii.gz')
        smoothnessFile = os.path.join(self.feat_dir, 'stats', 'smoothness')

        # Load DLH, VOLUME and RESELS
        smoothness = np.loadtxt(smoothnessFile, usecols=[1])

        search_space = SearchSpace(search_space_file=search_space_file, search_volume=int(smoothness[1]),
         resel_size_in_voxels=float(smoothness[2]), dlh=float(smoothness[0]),
         coordinate_system=self.coordinate_system, coordinate_space_id=self.coordinate_space_id,
         export_dir=self.export_dir)

        return search_space


    # Create excursion set, clusters and peaks entities
    def _get_clusters_peaks(self, stat_num):

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

        peaks = dict()
        prev_cluster = -1
        if (peak_file is not None) and (peak_std_file is not None):

            peaks_join_table = np.column_stack((peak_table, peak_std_table))
            for peak_row in peaks_join_table:    
                cluster_id = int(peak_row[0])

                if not cluster_id == prev_cluster:
                    # First peak in this cluster
                    peakIndex = 1;
                # Though peak coordinates in voxels are integer, we use a float type to comply with the rdfs:range
                peak = Peak(peak_index=int(peakIndex), x=float(peak_row[2]), y=float(peak_row[3]), z=float(peak_row[4]), 
                    x_std=float(peak_row[7]), y_std=float(peak_row[8]), z_std=float(peak_row[9]),
                    equiv_z=float(peak_row[1]), cluster_index=cluster_id, stat_num=stat_num, max_peak=(peakIndex==1))
                if cluster_id in peaks:
                    peaks[cluster_id].append(peak)
                else:
                    peaks[cluster_id] = list([peak])
                
                prev_cluster = cluster_id

                peakIndex = peakIndex + 1
        elif (peak_file is not None):
            for peak_row in peak_table:    
                cluster_id = int(peak_row[0])

                if not cluster_id == prev_cluster:
                    peakIndex = 1;

                peak = Peak(peak_index=int(peakIndex), x=int(peak_row[2]), y=int(peak_row[3]), z=int(peak_row[4]),
                    equiv_z=float(peak_row[1]), cluster_index=cluster_id, stat_num=stat_num, max_peak=(peakIndex==1))
                if cluster_id in peaks:
                    peaks[cluster_id].append(peak)
                else:
                    peaks[cluster_id] = list([peak])             

                prev_cluster = cluster_id

                peakIndex = peakIndex + 1
        elif (peak_std_file is not None):
            for peak_row in peak_std_table:    
                cluster_id = int(peak_row[0])

                if not cluster_id == prev_cluster:
                    peakIndex = 1;

                peak = Peak(peak_index=int(peakIndex), x_std=int(peak_row[2]), y_std=int(peak_row[3]), 
                    z_std=int(peak_row[4]), equiv_z=float(peak_row[1]), cluster_index=cluster_id, 
                    stat_num=stat_num, max_peak=(peakIndex==1))
                if cluster_id in peaks:
                    peaks[cluster_id].append(peak)
                else:
                    peaks[cluster_id] = list([peak])
                prev_cluster = cluster_id

                peakIndex = peakIndex + 1     
        
        clusters = list()

        if (cluster_file is not None) and (cluster_std_file is not None):
            clusters_join_table = np.column_stack((cluster_table, cluster_std_table))
            for cluster_row in clusters_join_table: 
                cluster_id = int(cluster_row[0])
                size = int(cluster_row[1])
                pFWER = float(cluster_row[2])
                x = float(cluster_row[8])
                y = float(cluster_row[9])
                z = float(cluster_row[10])
                x_std = float(cluster_row[24])
                y_std = float(cluster_row[25])
                z_std = float(cluster_row[26])
                clusters.append(Cluster(cluster_id=cluster_id, size=size, pFWER=pFWER, peaks=peaks[cluster_id],
                    x=x,y=y,z=z,x_std=x_std,y_std=y_std,z_std=z_std))
                # clusters[cluster_row[0]] = Cluster(cluster_id=int(cluster_row[0]), size=int(cluster_row[1]), pFWER=float(cluster_row[2]),
                #     x=float(cluster_row[8]),y=float(cluster_row[9]),z=float(cluster_row[10]),
                #     x_std=float(cluster_row[24]),y_std=float(cluster_row[25]),z_std=float(cluster_row[26]))
        elif (cluster_file is not None):
            for cluster_row in cluster_table: 
                cluster_id = int(cluster_row[0])
                size = int(cluster_row[1])
                pFWER = float(cluster_row[2])
                x = float(cluster_row[8])
                y = float(cluster_row[9])
                z = float(cluster_row[10])
                x_std = None
                y_std = None
                z_std = None
                clusters.append(Cluster(cluster_id=cluster_id, size=size, pFWER=pFWER, peaks=peaks[cluster_id],
                    x=x,y=y,z=z,x_std=x_std,y_std=y_std,z_std=z_std))
                # clusters[cluster_row[0]] = Cluster(id=int(cluster_row[0]), size=int(cluster_row[1]), pFWER=float(cluster_row[2]),
                    # x=float(cluster_row[8]),y=float(cluster_row[9]),z=float(cluster_row[10]))
        elif (cluster_std_file is not None):
            for cluster_row in cluster_std_table: 
                cluster_id = int(cluster_row[0])
                size = int(cluster_row[1])
                pFWER = float(cluster_row[2])
                x_std = float(cluster_row[8])
                y_std = float(cluster_row[9])
                z_std = float(cluster_row[10])
                x = None
                y = None
                z = None
                clusters.append(Cluster(cluster_id=cluster_id, size=size, pFWER=pFWER, peaks=peaks[cluster_id],
                    x=x,y=y,z=z,x_std=x_std,y_std=y_std,z_std=z_std))
                # clusters[cluster_row[0]] = Cluster(id=int(cluster_row[0]), size=int(cluster_row[1]), pFWER=float(cluster_row[2]),
                #     x_std=float(cluster_row[8]),y_std=float(cluster_row[9]),z_std=float(cluster_row[10]))

        return clusters
        
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



