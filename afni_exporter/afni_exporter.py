"""
Export neuroimaging results created with feat in AFNI following NIDM-Results 
specification.

@author: Rick Reynolds/Camille Maumet
@copyright: 
"""

import re
import os
import glob
import numpy as np
from exporter.exporter import NIDMExporter
from exporter.objects.constants import *
from exporter.objects.modelfitting import *
from exporter.objects.contrast import *
from exporter.objects.inference import *
import objects.afni_objects as afniobjs

class AFNItoNIDMExporter(NIDMExporter, object):
    """ 
    Parse an AFNI group result to extract the pieces information to be 
    stored in NIDM-Results and generate a NIDM-Results export.
    """

    def __init__(self, dset, csim_dset, p_uncor=0.01, p_cor=0.05,
                 nidm_ver="0.2.0"):
        self.stat_dset = dset
        self.clust_dset = csim_dset
        self.p_uncor = p_uncor
        self.p_cor = p_cor

        self.ind_contr = 0
        self.ind_stat  = 1

        self.afni_dir = os.path.dirname(self.stat_dset)

        nidm_dirs = glob.glob(os.path.join(self.afni_dir, 'nidm*'))
        if nidm_dirs:
            if nidm_dirs[-1] == os.path.join(self.afni_dir, 'nidm'):
                export_dir_num = 1
            else:
                m = re.search('(?<=nidm_).*', nidm_dirs[-1])
                export_dir_num = int(m.group(0))+1

            self.export_dir = os.path.join(self.afni_dir, \
                'nidm'+"_{0:0>4}".format(export_dir_num))
        else:
            self.export_dir = os.path.join(self.afni_dir, 'nidm')

        self.design_txt = None
        self.coordinate_system = None

    def parse(self):
        """ 
        Parse an AFNI result directory to extract the pieces information to be 
        stored in NIDM-Results.
        """    
        # rcr - ponder
        # design_file_open = open(self.design_file, 'r')
        # self.design_txt = design_file_open.read()
        

        # Retreive coordinate space used for current analysis
        if not self.coordinate_system:
            self._get_coordinate_system()

        super(AFNItoNIDMExporter, self).parse()

    def _add_namespaces(self):
        """ 
        Overload of parent _add_namespaces to add AFNI namespace.
        """
        super(AFNItoNIDMExporter, self).__init__()
        self.doc.add_namespace(AFNI)

    def _find_software(self):
        """ 
        Return an object of type Software describing the version of AFNI used to
        compute the current analysis.
        """
        #version_re = r'.*set fmri\(version\) (?P<info>\d+\.?\d+).*'
        #afni_version = self._search_in_fsf(version_re)

        version = os.system("afni -ver")   # get a string
        software = afniobjs.Software(version=version)

        return software

    def _find_model_fitting(self):
        """ 
        Parse AFNI result directory to retreive model fitting information. 
        Return a list of objects of type ModelFitting.
        """
        design_matrix = self._get_design_matrix()
        data = self._get_data()
        error_model = self._get_error_model()

        rms_map = self._get_residual_mean_squares_map()
        param_estimates = self._get_param_estimate_maps()
        mask_map = self._get_mask_map()
        grand_mean_map = self._get_grand_mean(mask_map.file)

        activity = self._get_model_parameters_estimations(error_model)

        model_fitting = ModelFitting(activity, design_matrix, data, 
            error_model, param_estimates, rms_map, mask_map, grand_mean_map)

        return list([model_fitting])

    def _find_contrasts(self):
        """ 
        Parse AFNI result directory to retreive information about contrasts. 
        Return a dictionary of (key, value) pairs where key is a tuple 
        containing the identifier of a ModelParametersEstimation object and a 
        tuple of identifiers of ParameterEstimateMap objects, and value is an 
        object of type Contrast.
        """
        # There is a single Model Parameters Estimations activity per feat 
        # directory, all contrasts will therefore depend on it.
        mf_id = self.model_fittings[0].activity.id

        # Degrees of freedom
        # FIXME: check what happens when more than one contrast is performed
        dof_file = open(os.path.join(self.afni_dir, 'stats', 'dof'), 'r')
        dof = float(dof_file.read())

        exc_sets = glob.glob(os.path.join(self.afni_dir, 'thresh_z*.nii.gz'))

        contrasts = dict()
        for filename in exc_sets:
            s = re.compile('zf?stat\d+')
            zstatnum = s.search(filename)
            zstatnum = zstatnum.group()

            if zstatnum.startswith("zstat"):
                stat_type = "T"
                con_num = zstatnum.replace('zstat', '')
            elif zstatnum.startswith("zfstat"):
                stat_type = "F"
                con_num = zstatnum.replace('zfstat', '')

            # If more than one excursion set is reported, we need to use 
            # an index in the file names of the file exported in nidm
            if len(exc_sets) > 1:
                stat_num = "_"+stat_type.upper()+"{0:0>3}".format(con_num)
            else:
                stat_num = ""    

            # Contrast name
            name_re = r'.*set fmri\(conname_real\.'+con_num+\
                                '\) "(?P<info>[\w\s><]+)".*'
            contrast_name = self._search_in_fsf(name_re)
            self.contrast_names_by_num[con_num] = contrast_name

            # Contrast estimation activity
            estimation = ContrastEstimation(con_num, contrast_name)

            # Contrast weights
            weights_re = r'.*set fmri\(con_real'+con_num+\
                                '\.\d+\) (?P<info>\d+)'
            weight_search = re.compile(weights_re)
            contrast_weights = str(re.findall(weight_search, 
                self.design_txt)).replace("'", '')

            weights = ContrastWeights(stat_num, contrast_name, 
                contrast_weights, stat_type)
           
            # Find which parameter estimates were used to compute the 
            # contrast
            pe_ids = list()
            pe_index = 1
            contrast_weights = contrast_weights.replace(' ', '')
            contrast_weights = contrast_weights.replace('[', '')
            contrast_weights = contrast_weights.replace(']', '')
            contrast_weights = contrast_weights.split(',')

            # Whenever a "1" is found in contrast_weights, the 
            # parameter estimate map identified by the corresponding 
            # index is in use
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
                pe_index += 1;

            # Convert to immutable tuple to be used as key
            pe_ids = tuple(pe_ids)

            # Contrast Map
            con_file = os.path.join(self.afni_dir, 
                'stats', 'cope'+str(con_num)+'.nii.gz')
            contrast_map = ContrastMap(con_file, stat_num, 
                contrast_name, self.coordinate_system, 
                self.coordinate_space_id, self.export_dir)
            self.coordinate_space_id += 1

            # Contrast Variance and Standard Error Maps
            varcontrast_file = os.path.join(self.afni_dir, 
                'stats', 'varcope'+str(con_num)+'.nii.gz')
            is_variance = True
            std_err_map = ContrastStdErrMap(stat_num, 
                varcontrast_file, is_variance, self.coordinate_system, 
                self.coordinate_space_id, self.export_dir)
            self.coordinate_space_id += 2

            # Statistic Map
            stat_file = os.path.join(self.afni_dir, 
                'stats', stat_type.lower()+'stat'+str(con_num)+'.nii.gz')
            stat_map = StatisticMap(stat_file, stat_type, stat_num, 
                contrast_name, dof, self.coordinate_system, 
                self.coordinate_space_id, self.export_dir)
            self.coordinate_space_id += 1

            # Z-Statistic Map
            z_stat_file = os.path.join(self.afni_dir, 
                'stats', 'zstat'+str(con_num)+'.nii.gz')
            z_stat_map = StatisticMap(z_stat_file, 'Z', stat_num, 
                contrast_name, dof, self.coordinate_system, 
                self.coordinate_space_id, self.export_dir)
            self.coordinate_space_id += 1

            con = Contrast(con_num, contrast_name, weights, estimation, 
                contrast_map, std_err_map, stat_map, z_stat_map)

            contrasts.setdefault((mf_id, pe_ids), list()).append(con)

        return contrasts

    def _find_inferences(self):
        """ 
        Parse AFNI result directory to retreive information about inference 
        along with peaks and clusters. Return a dictionary of (key, value) 
        pairs where key is the identifier of a ContrastEstimation object and 
        value is an object of type Inference.
        """
        inferences = dict()

        exc_sets = glob.glob(os.path.join(self.afni_dir, 'thresh_z*.nii.gz'))

        # Find excursion sets (in a given feat directory we have one excursion 
        # set per contrast)
        for filename in exc_sets:
            s = re.compile('zf?stat\d+')
            zstatnum = s.search(filename)
            zstatnum = zstatnum.group()
            if zstatnum.startswith("zstat"):
                stat_type = "T"
                stat_num = zstatnum.replace('zstat', '')
            elif zstatnum.startswith("zfstat"):
                stat_type = "F"
                stat_num = zstatnum.replace('zfstat', '')

            # If more than one excursion set is reported, we need to use 
            # an index in the file names of the file exported in nidm
            if len(exc_sets) > 1:
                stat_num_t = "_"+stat_type.upper()+"{0:0>3}".format(stat_num)
            else:
                stat_num_t = ""  

            # Find corresponding contrast estimation activity
            for contrasts in self.contrasts.values():
                for contrast in contrasts:
                    s = re.compile('cope\d+')
                    con_num = s.search(contrast.contrast_map.file)
                    con_num = con_num.group()
                    con_num = con_num.replace('cope', '')
                    if con_num == stat_num:
                        con_id = contrast.estimation.id

            # Inference activity
            inference_act = InferenceActivity(stat_num, 
                self.contrast_names_by_num[stat_num])

            # Excursion set
            visualisation = os.path.join(self.afni_dir, 
                'rendered_thresh_zstat'+stat_num+'.png')
            zFileImg = os.path.join(self.afni_dir, 
                'thresh_zstat'+stat_num+'.nii.gz')
            exc_set = ExcursionSet(zFileImg, stat_num_t, visualisation, 
                self.coordinate_system, self.coordinate_space_id, 
                self.export_dir)
            self.coordinate_space_id += 1

            # Height Threshold
            prob_re = r'.*set fmri\(prob_thresh\) (?P<info>\d+\.?\d+).*'
            z_re = r'.*set fmri\(z_thresh\) (?P<info>\d+\.?\d+).*'
            type_re = r'.*set fmri\(thresh\) (?P<info>\d+).*'

            prob_thresh = float(self._search_in_fsf(prob_re))
            z_thresh = float(self._search_in_fsf(z_re))
            thresh_type = self._search_in_fsf(type_re)

            # FIXME: deal with 0 = no thresh?
            voxel_uncorr = (thresh_type == 1)
            voxel_corr = (thresh_type == 2)
            # cluster_thresh = (thresh_type == 3)
            
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

            height_thresh = HeightThreshold(stat_threshold, 
                p_corr_threshold, p_uncorr_threshold)

            # Extent Threshold
            extent_thresh = ExtentThreshold(p_corr=extent_p_corr)

            # Clusters (and associated peaks)
            clusters = self._get_clusters_peaks(stat_num)

            # Peak and Cluster Definition Criteria
            peak_criteria = PeakCriteria(stat_num, 
                self._get_num_peaks(), self._get_peak_dist())
            clus_criteria = ClusterCriteria(stat_num, 
                self._get_connectivity())

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
            contrast_masking_file = self._get_display_mask()
            display_mask = DisplayMaskMap(stat_num, 
                contrast_masking_file, self.coordinate_system, 
                self.coordinate_space_id, self.export_dir)
            self.coordinate_space_id += 1

            # Search space
            search_space = self._get_search_space()

            inference = Inference(inference_act, height_thresh, 
                extent_thresh, peak_criteria, clus_criteria, 
                display_mask, exc_set, clusters, search_space, 
                self.software.id)

            inferences.setdefault(con_id, list()).append(inference)

        return inferences

    def _get_design_matrix(self):
        """ 
        Parse AFNI result directory to retreive information about the design 
        matrix. Return an object of type DesignMatrix.
        """
        # rcr - DesignMatrix looks FSLish, may need alteration
        # design_matrix = DesignMatrix(None, None)
        design_matrix = None
        return design_matrix

    def _get_data(self):
        """ 
        Parse AFNI result directory to retreive information about the data. 
        Return an object of type Data. 
        """
        grand_mean_scaling = False
        target_intensity = 1.0  # rcr - maybe get from first level
        data = Data(grand_mean_scaling, target_intensity)
        return data

    def _get_error_model(self):
        """ 
        Parse AFNI result directory to retreive information about the error 
        model. Return an object of type ErrorModel. 
        """
        # fmri_level_re = r'.*set fmri\(level\) (?P<info>\d+).*'
        fmri_level = 2 # group, for now
        self.first_level = (fmri_level == 1)

        if self.first_level:    # rcr - update
            variance_homo = True
            dependance = SERIALLY_CORR
            variance_spatial = SPATIALLY_LOCAL
            dependance_spatial = SPATIALLY_LOCAL
        else:
            variance_homo = True                # rcr - for ttest
            dependance = INDEPEDENT_CORR
            variance_spatial = SPATIALLY_LOCAL
            dependance_spatial = None

        error_distribution = GAUSSIAN_DISTRIBUTION
        error_model = ErrorModel(error_distribution, variance_homo, 
            variance_spatial, dependance, dependance_spatial)
        return error_model

    def _get_residual_mean_squares_map(self):
        """ 
        Parse AFNI result directory to retreive information about the residual 
        mean squares map. Return an object of type ResidualMeanSquares. 
        """

        # rcr - finish
        # rms_map = ResidualMeanSquares(self.export_dir, residuals_file, 
        #     self.coordinate_system, self.coordinate_space_id)
        rms_map = None

        return rms_map

    def _get_param_estimate_maps(self):
        """ 
        Parse AFNI result directory to retreive information about the parameter
        estimates. Return a list of objects of type ParameterEstimateMap. 
        """

        # create temporary set of NIFTI volumes of betas output by ttest
        # (for A-B, return one for A and for B)
        # then delete the temp files
        param_estimates = list()
        # rcr - ponder penum (beta index, 1-based?)
        for filename, ind in enumerate(os.listdir(os.path.join(self.afni_dir, 'stats'))):
            # if AFNI format, convert to NIFTI as temp file
            full_path_file = NIFTI_NAME # rcr
            penum = ind 
            param_estimate = ParameterEstimateMap(full_path_file, 
                penum, self.coordinate_space_id, 
                self.coordinate_system)
            self.coordinate_space_id = self.coordinate_space_id + 1
            param_estimates.append(param_estimate)
        return param_estimates

    def _get_mask_map(self):
        """ 
        Parse AFNI result directory to retreive information about the mask 
        created as part of Model Parameters Estimation. Return an object of 
        type MaskMap. 
        """
        # this is infered mask
        # (elsewhere: would want to know custom_mask (input vs implied))

        mask_file = os.path.join(self.afni_dir, 'mask.nii.gz')
        mask_map = MaskMap(self.export_dir, mask_file,
            self.coordinate_system, self.coordinate_space_id)
        self.coordinate_space_id = self.coordinate_space_id + 1
        return mask_map       

    def _get_grand_mean(self, mask_file):
        """ 
        Parse AFNI result directory to retreive information about the grand 
        mean map. Return an object of type GrandMeanMap. 
        """
        grand_mean_file = os.path.join(self.afni_dir, 'mean_func.nii.gz')

        # FIXME: Check if there is an alternative file to use here (maybe)
        # depending on AFNI version
        if not os.path.isfile(grand_mean_file):
            grand_mean = None;
        else:
            grand_mean = GrandMeanMap(grand_mean_file, mask_file, 
                self.coordinate_system, self.coordinate_space_id, 
                self.export_dir)
            self.coordinate_space_id = self.coordinate_space_id + 1
        
        return grand_mean

    def _get_coordinate_system(self):
        """ 
        Parse AFNI result directory to retreive information about the 
        coordinate system used in the current analysis (dependent on the 
        template).
        """
        space_re = r'.*set fmri\(regstandard_yn\) (?P<info>[\d]+).*'
        standard_space = bool(self._search_in_fsf(space_re))

        if standard_space:
            custom_re = r'.*set fmri\(alternateReference_yn\) (?P<info>[\d]+).*'
            custom_space = (self._search_in_fsf(custom_re) == "1")

            if custom_space is not None:
                custom_standard = (custom_space == "1");
            else:
                custom_re = r'.*set fmri\(regstandard\) (?P<info>.+).*'
                custom_space = self._search_in_fsf(custom_re)

                if custom_space is not None:
                    custom_standard = True;
        else:
            custom_standard = False;

        if not standard_space:
            coordinate_system = NIDM['SubjectSpace'];
        else:
            if not custom_standard:
                coordinate_system = \
                    NIDM['IcbmMni152NonLinear6thGenerationCoordinateSystem'];
            else:
                coordinate_system = NIDM['StandarizedSpace'];

        self.coordinate_system = coordinate_system

    def _search_in_fsf(self, regexp):
        """ 
        Look for information matching regular expression 'regexp' in the design
        file of the current study.
        """
        # rcr - what are we looking for, group labels?
        if self.design_txt != None:
           info_search = re.compile(regexp)
           info_found = info_search.search(self.design_txt)
           info = info_found.group('info')
        else:
           info = ''
        return info

    def _get_display_mask(self):
        """ 
        Parse AFNI result directory to retreive information about display mask. 
        Return an object of type ResidualMeanSquares. 
        """
        # FIXME this should be updated with actual contrast masking file
        mask_file = os.path.join(self.afni_dir, 'mask.nii.gz')
        return mask_file

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
                # FIXME: for now omitted if not explicitely defined
                num_peak = None
        return num_peak

    def _get_peak_dist(self):
        peak_dist_search = re.compile(r'.* --peakdist=(?P<peakdist>\d+)+ .*')
        peak_dist_found = peak_dist_search.search(self.feat_post_log)
        if peak_dist_found:
            peak_dist = float(peak_dist_found.group('peakdist'))
        else:
            peak_dist = 0.0

        return peak_dist

    def _get_connectivity(self):
        """ 
        Parse AFNI result directory to retreive peak connectivity within a 
        cluster. 
        """
        conn_re = r'.* --connectivity=(?P<connectivity>\d+)+ .*'
        connectivity_search = re.compile(conn_re)
        connectivity = int(connectivity_search.search(self.feat_post_log).group('connectivity'))

        return connectivity

    def _get_search_space(self):
        """ 
        Parse AFNI result directory to retreive information about the search 
        space. Return an object of type SearchSpace. 
        """
        # FIXME this needs to be estimated
        search_space_file = os.path.join(self.afni_dir, 'mask.nii.gz')
        smoothness_file = os.path.join(self.afni_dir, 'stats', 'smoothness')

        # Load DLH, VOLUME and RESELS
        smoothness = np.loadtxt(smoothness_file, usecols=[1])

        search_space = SearchSpace(search_space_file=search_space_file, 
            search_volume=int(smoothness[1]),
            resel_size_in_voxels=float(smoothness[2]), 
            dlh=float(smoothness[0]), 
            random_field_stationarity=True,
            coordinate_system=self.coordinate_system, 
            coordinate_space_id=self.coordinate_space_id,
            export_dir=self.export_dir)
        self.coordinate_space_id += 1

        return search_space

    def _get_clusters_peaks(self, stat_num):
        """ 
        Parse AFNI result directory to retreive information about the clusters 
        and peaks declared significant for statistic 'stat_num'. Return a list
        of Cluster objects.
        """

        # Cluster list (positions in voxels)
        cluster_file = os.path.join(self.afni_dir, 
            'cluster_zstat'+stat_num+'.txt')
        if not os.path.isfile(cluster_file):
            cluster_file = None;
        else:
            cluster_table = np.loadtxt(cluster_file, skiprows=1, ndmin=2)

        # Cluster list (positions in mm)
        cluster_std_file = os.path.join(self.afni_dir, 
            'cluster_zstat'+stat_num+'_std.txt')
        if not os.path.isfile(cluster_std_file):
            cluster_std_file = None;
            # cluster_std_table = np.zeros_like(cluster_table)*float('nan')
        else:
            cluster_std_table = np.loadtxt(cluster_std_file, skiprows=1, 
                ndmin=2) 
        
        # Peaks
        peak_file = os.path.join(self.afni_dir, 'lmax_zstat'+stat_num+'.txt')
        if not os.path.isfile(peak_file):
            peak_file = None;
        else:
            peak_table = np.loadtxt(peak_file, skiprows=1, ndmin=2)

        peak_std_file = os.path.join(self.afni_dir, 
            'lmax_zstat'+stat_num+'_std.txt')
        if not os.path.isfile(peak_std_file):
            peak_std_file = None;  
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

                # Though peak coordinates in voxels are integer, we use a float 
                # type to comply with the rdfs:range
                peak = Peak(peak_index=int(peakIndex), x=float(peak_row[2]), 
                    y=float(peak_row[3]), z=float(peak_row[4]), 
                    x_std=float(peak_row[7]), y_std=float(peak_row[8]), 
                    z_std=float(peak_row[9]), equiv_z=float(peak_row[1]), 
                    cluster_index=cluster_id, stat_num=stat_num, 
                    max_peak=(peakIndex==1))
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

                peak = Peak(peak_index=int(peakIndex), x=int(peak_row[2]), 
                    y=int(peak_row[3]), z=int(peak_row[4]),
                    equiv_z=float(peak_row[1]), cluster_index=cluster_id, 
                    stat_num=stat_num, max_peak=(peakIndex==1))
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

                peak = Peak(peak_index=int(peakIndex), x_std=int(peak_row[2]), 
                    y_std=int(peak_row[3]), z_std=int(peak_row[4]), 
                    equiv_z=float(peak_row[1]), cluster_index=cluster_id, 
                    stat_num=stat_num, max_peak=(peakIndex==1))
                if cluster_id in peaks:
                    peaks[cluster_id].append(peak)
                else:
                    peaks[cluster_id] = list([peak])
                prev_cluster = cluster_id

                peakIndex = peakIndex + 1     
        
        clusters = list()

        if (cluster_file is not None) and (cluster_std_file is not None):
            clusters_join_table = np.column_stack((cluster_table, 
                cluster_std_table))
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
                clusters.append(Cluster(cluster_num=cluster_id, size=size, 
                    pFWER=pFWER, peaks=peaks[cluster_id], x=x, y=y, z=z,
                    x_std=x_std, y_std=y_std, z_std=z_std))

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
                clusters.append(Cluster(cluster_num=cluster_id, size=size, 
                    pFWER=pFWER, peaks=peaks[cluster_id], x=x, y=y, z=z,
                    x_std=x_std,y_std=y_std,z_std=z_std))
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
                clusters.append(Cluster(cluster_num=cluster_id, size=size, 
                    pFWER=pFWER, peaks=peaks[cluster_id], x=x, y=y, z=z,
                    x_std=x_std,y_std=y_std,z_std=z_std))

        return clusters
