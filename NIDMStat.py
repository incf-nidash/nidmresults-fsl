'''Python implementation of NI-DM (for statistical results)

@author: Camille Maumet <c.m.j.maumet@warwick.ac.uk>
@copyright: University of Warwick 2013-2014
'''

from prov.model import ProvBundle, ProvDocument, Namespace, ProvRecord, ProvEntity, Identifier
from prov.model import PROV
import os
import numpy as np
import nibabel as nib
import hashlib
import shutil


NIDM = Namespace('nidm', "http://www.incf.org/ns/nidash/nidm#")
NIIRI = Namespace("niiri", "http://iri.nidash.org/")
CRYPTO = Namespace("crypto", "http://id.loc.gov/vocabulary/preservation/cryptographicHashFunctions#")
FSL = Namespace("fsl", "http://www.incf.org/ns/nidash/fsl#")
DCT = Namespace("dct", "http://purl.org/dc/terms/")

''' Create a NIDM export file and copy related nifti files
'''
class NIDMStat():

    def __init__(self, *args, **kwargs):
        # FIXME: Merge coordinateSpace entities if identical?
        # FIXME: Use actual URIs instead

        # Directory in which export will be stored
        self.export_dir = kwargs.pop('export_dir')
        if not os.path.exists(self.export_dir):
            os.makedirs(self.export_dir)
        # FIXME: Do something if export dir already exist

        # Keep track of the number of coordinateSpace entity already generated
        self.coordinateSpaceId = 0

        # Create namespaces
        self.provDocument = ProvDocument();
        self.provDocument.add_namespace("neurolex", "http://neurolex.org/wiki/")
        self.provDocument.add_namespace(FSL)
        self.provDocument.add_namespace(NIDM)
        self.provDocument.add_namespace(NIIRI)
        self.provDocument.add_namespace(CRYPTO)

        self.provBundle = ProvBundle(identifier=NIIRI['fsl_results_id'])
        
        self.standard_space = kwargs.pop('standard_space')
        self.custom_standard = kwargs.pop('custom_standard')
       
        # FIXME: Check one-tailed or two-tailed test and get test type from data
               
        # FIXME: is this really empty? If so, should be deleted form the model
        # g.entity(NIIRI['stat_image_properties_id'], other_attributes=( 
        #     (PROV['type'], FSL['statisticImageProperties']), 
        #     (PROV['label'], 'Statistical image properties')))
        
    def create_software(self, feat_version):
        # Add software agent: FSL
        self.provBundle.agent(NIIRI['software_id'], other_attributes=( 
            (PROV['type'], NIDM['FSL']), 
            (PROV['type'], PROV['SoftwareAgent']),
            (PROV['label'],'FSL'),
            # FIXME find FSL software version
            (FSL['featVersion'], feat_version) ))
        
    def create_thresholds(self, *args, **kwargs):
        voxel_threshold = kwargs.pop('voxel_threshold')
        voxel_p_uncorr = kwargs.pop('voxel_p_uncorr')
        voxel_p_corr = kwargs.pop('voxel_p_corr')
        thresh_desc = ""
        if not voxel_threshold is None:
            thresh_desc = "Z>"+str(voxel_threshold)
            user_threshold_type = NIDM['zStatistic']
        elif not voxel_p_uncorr is None:
            thresh_desc = "p<"+str(voxel_p_uncorr)+" uncorr."
            user_threshold_type = NIDM['pValueUncorrected']
        elif not voxel_p_corr is None:
            thresh_desc = "p<"+str(voxel_p_corr)+" (GRF)"
            user_threshold_type = NIDM['pValueFWER']

        # FIXME: Do we want to calculate an uncorrected p equivalent to the Z thresh? 
        # FIXME: Do we want/Can we find a corrected p equivalent to the Z thresh? 
        heightThreshAllFields = {
            PROV['type']: NIDM['HeightThreshold'], PROV['label']: "Height Threshold: "+thresh_desc,
            NIDM['userSpecifiedThresholdType']: user_threshold_type , PROV['value']: voxel_threshold,
            NIDM['pValueUncorrected']: voxel_p_uncorr, NIDM['pValueFWER']: voxel_p_corr
            }
        self.provBundle.entity(NIIRI['height_threshold_id'], other_attributes=dict((k,v) for k,v in heightThreshAllFields.iteritems() if v is not None))

        extent = kwargs.pop('extent')
        extent_p_uncorr = kwargs.pop('extent_p_uncorr')
        extent_p_corr = kwargs.pop('extent_p_corr')
        thresh_desc = ""
        if not extent is None:
            thresh_desc = "k>"+str(extent)
            user_threshold_type = NIDM['clusterSizeInVoxels']
        elif not extent_p_uncorr is None:
            thresh_desc = "p<"+str(extent_p_uncorr)+" uncorr."
            user_threshold_type = NIDM['pValueUncorrected']
        elif not extent_p_corr is None:
            thresh_desc = "p<"+str(extent_p_corr)+" corr."
            user_threshold_type = NIDM['pValueFWER']
        extent_thresh_all_fields = {
            PROV['type']: NIDM['ExtentThreshold'], PROV['label']: "Extent Threshold: "+thresh_desc, NIDM['clusterSizeInVoxels']: extent,
            NIDM['userSpecifiedThresholdType']: user_threshold_type, NIDM['pValueUncorrected']: extent_p_uncorr, NIDM['pValueFWER']: extent_p_corr
        }
        self.provBundle.entity(NIIRI['extent_threshold_id'], other_attributes=dict((k,v) for k,v in extent_thresh_all_fields.iteritems() if v is not None))

    def create_coordinate(self, coordinate_id, label_id, x=None, y=None, z=None, x_std=None, y_std=None, z_std=None):
        # We can not have this in the dictionnary because we want to keep the duplicate prov:type attribute
        typeAndLabelAttributes = [(PROV['type'],PROV['Location']),
            (PROV['type'], NIDM['Coordinate']),
            (PROV['label'], "Coordinate "+label_id)]

        coordinateAttributes = {
            FSL['coordinate1InVoxels'] : x,
            FSL['coordinate2InVoxels'] : y,
            FSL['coordinate3InVoxels'] : z,
            NIDM['coordinate1'] : x_std,
            NIDM['coordinate2'] : y_std,
            NIDM['coordinate3'] : z_std
            };

        self.provBundle.entity(coordinate_id, 
            other_attributes=typeAndLabelAttributes+list(dict((k,v) for k,v in coordinateAttributes.iteritems() if not v is None).items()))

    def get_sha_sum(self, nifti_file):
        nifti_img = nib.load(nifti_file)
        data = nifti_img.get_data()
        # Fix needed as in https://github.com/pymc-devs/pymc/issues/327
        if not data.flags["C_CONTIGUOUS"]:
          data = np.ascontiguousarray(data)
        return hashlib.sha512(data).hexdigest()

    def create_cluster(self, stat_num, id, size, pFWER, *args, **kwargs):
        clusterIndex = id
        stat_num = stat_num

        # FIXME deal with multiple contrasts
        cluster_id = clusterIndex

        self.provBundle.entity(NIIRI['cluster_000'+str(cluster_id)], other_attributes=( 
                             (PROV['type'] , NIDM['Cluster']), 
                             (PROV['label'], "Cluster 000"+str(cluster_id)),
                             (NIDM['clusterSizeInVoxels'], size),
                             (NIDM['pValueFWER'], pFWER )))
        self.provBundle.wasDerivedFrom(NIIRI['cluster_000'+str(cluster_id)], NIIRI['excursion_set_id_'+str(stat_num)])

        self.create_coordinate(NIIRI['COG_coordinate_000'+str(cluster_id)], '000'+str(cluster_id),**kwargs)

        self.provBundle.entity(NIIRI['center_of_gravity_'+str(cluster_id)], other_attributes=( 
                     (PROV['type'] , FSL['CenterOfGravity']), 
                     (PROV['label'], "Center of gravity "+str(cluster_id)),
                     (PROV['location'] , NIIRI['COG_coordinate_000'+str(cluster_id)]))   )
        self.provBundle.wasDerivedFrom(NIIRI['center_of_gravity_'+str(cluster_id)], NIIRI['cluster_000'+str(cluster_id)])

    def create_peak(self, id, cluster_id, equivZ, stat_num, max_peak, *args, **kwargs):
        peakIndex = id
        clusterIndex = cluster_id
        stat_num = stat_num

        # FIXME: Currently assumes less than 10 clusters per contrast
        cluster_id = clusterIndex

        # FIXME: Currently assumes less than 100 peaks 
        peakUniqueId = '000'+str(clusterIndex)+'_'+str(peakIndex)

        self.create_coordinate(NIIRI['coordinate_'+str(peakUniqueId)], str(peakUniqueId), **kwargs)

        other_attributes = [ 
            (PROV['type'] , NIDM['Peak']), 
            (PROV['label'] , "Peak "+str(peakUniqueId)), 
            (NIDM['equivalentZStatistic'], equivZ), 
            (PROV['location'] , NIIRI['coordinate_'+str(peakUniqueId)])]

        if max_peak:
            other_attributes.insert(0, (PROV['type'], FSL['ClusterMaximumStatistic']))

        self.provBundle.entity(NIIRI['peak_'+str(peakUniqueId)], other_attributes=other_attributes)
        self.provBundle.wasDerivedFrom(NIIRI['peak_'+str(peakUniqueId)], NIIRI['cluster_000'+str(cluster_id)])

    def create_model_fitting(self, residuals_file, grand_mean_file, design_matrix):
        # Copy residuals map in export directory
        if not residuals_file is None:
            original_file = residuals_file
            residuals_file = os.path.join(self.export_dir, 'ResidualMeanSquares.nii.gz')
            shutil.copy(original_file, residuals_file)
            path, residuals_filename = os.path.split(residuals_file)
            residuals_file = os.path.join(self.export_dir,residuals_filename)  

        # Create "Model Parameter estimation" activity
        self.provBundle.activity(NIIRI['model_parameters_estimation_id'], other_attributes=( 
            (PROV['type'], NIDM['ModelParametersEstimation']),(PROV['label'], "Model Parameters Estimation")))
        self.provBundle.wasAssociatedWith(NIIRI['model_parameters_estimation_id'], NIIRI['software_id'])

        # Create "Data" entity
        # FIXME: grand mean scaling?
        # FIXME: medianIntensity
        self.provBundle.entity(NIIRI['data_id'],
            other_attributes=(  (PROV['type'],NIDM['Data']),
                                (PROV['type'],PROV['Collection']),
                                (PROV['label'],"Data"),
                                (NIDM['grandMeanScaling'], True),
                                (NIDM['targetIntensity'], 10000.0)))
        self.provBundle.used(NIIRI['model_parameters_estimation_id'], NIIRI['data_id'])

        if not residuals_file is None:
            # Create "residuals map" entity
            self.provBundle.entity(NIIRI['residual_mean_squares_map_id'], 
                other_attributes=( (PROV['type'],NIDM['ResidualMeanSquaresMap'],), 
                                   (DCT['format'], NIDM['Nifti1Gzip']), 
                                   (PROV['location'], Identifier("file://./"+residuals_filename) ),
                                   (PROV['label'],"Residual Mean Squares Map" ),
                                   (NIDM['originalFileName'],residuals_filename ),
                                   (CRYPTO['sha512'], self.get_sha_sum(residuals_file)),
                                   (NIDM['atCoordinateSpace'], self.create_coordinate_space(residuals_file))))
            self.provBundle.wasGeneratedBy(NIIRI['residual_mean_squares_map_id'], NIIRI['model_parameters_estimation_id'])  
        else:
            # FIXME: Replace with log
            print "No residual file"

        # Grand Mean Map entity
        original_file = grand_mean_file
        grand_mean_file = os.path.join(self.export_dir, 'GrandMean.nii.gz')
        shutil.copy(original_file, grand_mean_file)
        path, grand_mean_filename = os.path.split(grand_mean_file)
        grand_mean_file = os.path.join(self.export_dir,grand_mean_filename)  

        # TODO: nidm:maskedMedian "115"^^xsd:int ;
        self.provBundle.entity(NIIRI['grand_mean_map_id'], 
            other_attributes=( (PROV['type'],NIDM['GrandMeanMap']), 
                               (DCT['format'], NIDM['Nifti1Gzip']), 
                               (PROV['label'],"Grand Mean Map"), 
                               (NIDM['originalFileName'], grand_mean_filename),
                               (NIDM['atCoordinateSpace'], self.create_coordinate_space(grand_mean_file)),
                               (CRYPTO['sha512'], self.get_sha_sum(grand_mean_file)),
                               (PROV['location'], Identifier("file://./"+grand_mean_filename))))      
        self.provBundle.wasGeneratedBy(NIIRI['grand_mean_map_id'], NIIRI['model_parameters_estimation_id'],)
        
        # Create cvs file containing design matrix
        design_matrix_csv = 'design_matrix.csv'
        np.savetxt(os.path.join(self.export_dir, design_matrix_csv), np.asarray(design_matrix), delimiter=",")

        # Create "design matrix" entity
        self.provBundle.entity(NIIRI['design_matrix_id'], 
            other_attributes=( (PROV['type'],NIDM['DesignMatrix']), 
                               (PROV['label'],"Design Matrix"), 
                               # (NIDM['originalFileName'],design_matrix_csv ),
                               (PROV['location'], Identifier("file://./"+design_matrix_csv))))       
        self.provBundle.used(NIIRI['model_parameters_estimation_id'], NIIRI['design_matrix_id'])

    # Generate prov for contrast map
    def create_parameter_estimate(self, pe_file, pe_num):
        # Copy parameter estimate map in export directory
        # shutil.copy(pe_file, self.export_dir)
        path, pe_filename = os.path.split(pe_file)
        # pe_file = os.path.join(self.export_dir,pe_filename)       

        # Parameter estimate entity
        self.provBundle.entity(NIIRI['beta_map_id_'+str(pe_num)], 
            other_attributes=( (PROV['type'], NIDM['ParameterEstimateMap']), 
                               # (DCT['format'], NIDM['Nifti1Gzip']), 
                               # (PROV['location'], Identifier("file://./"+pe_filename)),
                               (NIDM['originalFileName'], pe_filename), 
                               (NIDM['atCoordinateSpace'], self.create_coordinate_space(pe_file)),
                               # (CRYPTO['sha512'], self.get_sha_sum(pe_file)),
                               (PROV['label'], "Parameter estimate "+str(pe_num))))
        
        self.provBundle.wasGeneratedBy(NIIRI['beta_map_id_'+str(pe_num)], NIIRI['model_parameters_estimation_id'])  

    # Generate prov for contrast map
    def create_contrast_map(self, cope_file, var_cope_file, stat_file, z_stat_file, contrast_name, contrast_num, dof, contrastWeights):
        # Contrast id entity
        # FIXME: Get contrast weights
        # FIXME: Deal with F weights
        self.provBundle.entity(NIIRI['contrast_id_'+contrast_num], 
            other_attributes=( (PROV['type'], NIDM['ContrastWeights']), 
                               (NIDM['statisticType'], NIDM['TStatistic']),
                               (PROV['label'], "Contrast Weights: "+contrast_name), 
                               (NIDM['contrastName'], contrast_name),
                               (PROV['value'], contrastWeights)))

        # Create related activities
        self.provBundle.activity(NIIRI['contrast_estimation_id_'+contrast_num], other_attributes=( 
            (PROV['type'], NIDM['ContrastEstimation']),
            (PROV['label'], "Contrast estimation: "+contrast_name)))

        # Find which betas were used to compute the contrast
        contrastWeights = contrastWeights.replace(' ', '').replace('[', '').replace(']', '').split(',')
        peIndex = 1
        for betaIndex in contrastWeights:
            if int(betaIndex) == 1:
                self.provBundle.used(NIIRI['contrast_estimation_id_'+contrast_num], NIIRI['beta_map_id_'+str(peIndex)])
            peIndex += 1;

        # Copy contrast map in export directory
        original_file = cope_file
        cope_file = os.path.join(self.export_dir, 'Contrast.nii.gz')
        shutil.copy(original_file, cope_file)
        path, cope_filename = os.path.split(cope_file)
        cope_file = os.path.join(self.export_dir,cope_filename)  

        # Contrast Map entity
        path, filename = os.path.split(cope_file)
        self.provBundle.entity('niiri:'+'contrast_map_id_'+contrast_num, other_attributes=( 
            (PROV['type'], NIDM['ContrastMap']), 
            (DCT['format'], NIDM['Nifti1Gzip']), 
            (NIDM['atCoordinateSpace'], self.create_coordinate_space(cope_file)),
            (PROV['location'], Identifier("file://./"+cope_filename)),
            (NIDM['originalFileName'], cope_filename),
            (NIDM['contrastName'], contrast_name),
            (CRYPTO['sha512'], self.get_sha_sum(cope_file)),
            (PROV['label'], "Contrast Map: "+contrast_name)))
        
        self.provBundle.wasGeneratedBy(NIIRI['contrast_map_id_'+contrast_num], NIIRI['contrast_estimation_id_'+contrast_num])
        self.provBundle.wasAssociatedWith(NIIRI['contrast_estimation_id_'+contrast_num], NIIRI['software_id'])

        # Copy contrast variance map in export directory
        shutil.copy(var_cope_file, self.export_dir)
        path, var_cope_filename = os.path.split(var_cope_file)
        var_cope_file = os.path.join(self.export_dir,var_cope_filename)  

        # Contrast Variance Map entity
        self.provBundle.entity('niiri:'+'contrast_variance_map_id_'+contrast_num, other_attributes=( 
            (PROV['type'], NIDM['Map']), 
            # (DCT['format'], NIDM['Nifti1Gzip']), 
            (NIDM['atCoordinateSpace'], self.create_coordinate_space(var_cope_file)),
            # (PROV['location'], Identifier("file://./"+var_cope_filename)),
            (CRYPTO['sha512'], self.get_sha_sum(var_cope_file)),
            (NIDM['originalFileName'], var_cope_filename)))
            # (PROV['label'], "Contrast Variance Map "+contrast_num)))
        

        # Create standard error map from contrast variance map
        var_cope_img = nib.load(var_cope_file)
        contrast_variance = var_cope_img.get_data()

        standard_error_img = nib.Nifti1Image(np.sqrt(contrast_variance), var_cope_img.get_qform())
        standard_error_file = var_cope_file.replace('var', 'sqrt_var')
        nib.save(standard_error_img, standard_error_file)

        path, filename = os.path.split(standard_error_file)
        self.provBundle.entity('niiri:'+'contrast_standard_error_map_id_'+contrast_num, other_attributes=( 
            (PROV['type'], NIDM['ContrastStandardErrorMap']), 
            (DCT['format'], NIDM['Nifti1Gzip']), 
            (NIDM['atCoordinateSpace'], self.create_coordinate_space(standard_error_file)),
            (PROV['location'], Identifier("file://./"+filename)),
            (CRYPTO['sha512'], self.get_sha_sum(standard_error_file)),
            # (NIDM['originalFileName'], filename),
            (PROV['label'], "Contrast Standard Error Map")))
        
        self.provBundle.wasGeneratedBy(NIIRI['contrast_standard_error_map_id_'+contrast_num], NIIRI['contrast_estimation_id_'+contrast_num])
        self.provBundle.wasDerivedFrom(NIIRI['contrast_standard_error_map_id_'+contrast_num], NIIRI['contrast_variance_map_id_'+contrast_num])

        
        # FIXME: Remove TODOs

        # Copy Z-Statistical map in export directory
        original_file = z_stat_file
        z_stat_file = os.path.join(self.export_dir, 'ZStatistic.nii.gz')
        shutil.copy(original_file, z_stat_file)
        path, z_stat_filename = os.path.split(z_stat_file)
        z_stat_file = os.path.join(self.export_dir,z_stat_filename)       

        # Create "Z-Statistic Map" entity
        self.provBundle.entity(NIIRI['z_statistic_map_id_'+contrast_num ],
            other_attributes=(  (PROV['type'], NIDM['StatisticMap']), 
                                (DCT['format'], NIDM['Nifti1Gzip']), 
                                (PROV['label'], "Z-Statistic Map: "+contrast_name) ,
                                (PROV['location'], Identifier("file://./"+z_stat_filename)),
                                (NIDM['statisticType'], NIDM['ZStatistic']), 
                                (NIDM['contrastName'], contrast_name),
                                (NIDM['originalFileName'], z_stat_filename),
                                (CRYPTO['sha512'], self.get_sha_sum(z_stat_file)),
                                (NIDM['atCoordinateSpace'], self.create_coordinate_space(z_stat_file)),
                                ) )

        # Copy Statistical map in export directory
        original_file = stat_file
        stat_file = os.path.join(self.export_dir, 'TStatistic.nii.gz')
        shutil.copy(original_file, stat_file)
        path, stat_filename = os.path.split(stat_file)
        stat_file = os.path.join(self.export_dir,stat_filename)     

        # Create "Statistic Map" entity
        # FIXME: Deal with other than t-contrast maps: dof + statisticType
        self.provBundle.entity(NIIRI['statistic_map_id_'+contrast_num ],
            other_attributes=(  (PROV['type'], NIDM['StatisticMap']), 
                                (DCT['format'], NIDM['Nifti1Gzip']), 
                                (PROV['label'], "Statistic Map: "+contrast_name) ,
                                (PROV['location'], Identifier("file://./"+stat_filename)),
                                (NIDM['statisticType'], NIDM['TStatistic']), 
                                (NIDM['originalFileName'], stat_filename),
                                (NIDM['contrastName'], contrast_name),
                                (NIDM['errorDegreesOfFreedom'], dof),
                                (NIDM['effectDegreesOfFreedom'], 1.0),
                                (CRYPTO['sha512'], self.get_sha_sum(stat_file)),
                                (NIDM['atCoordinateSpace'], self.create_coordinate_space(stat_file)),
                                ) )
        
        self.provBundle.wasGeneratedBy(NIIRI['statistic_map_id_'+contrast_num], NIIRI['contrast_estimation_id_'+contrast_num])
               
        self.provBundle.wasGeneratedBy(NIIRI['z_statistic_map_id_'+contrast_num], NIIRI['contrast_estimation_id_'+contrast_num])
        self.provBundle.used(NIIRI['contrast_estimation_id_'+contrast_num], NIIRI['residual_mean_squares_map_id'])
        self.provBundle.used(NIIRI['contrast_estimation_id_'+contrast_num], NIIRI['design_matrix_id'])
        self.provBundle.used(NIIRI['contrast_estimation_id_'+contrast_num], NIIRI['contrast_id_'+contrast_num])


        # In FSL we have a single thresholding (extent, height) applied to all contrasts 
        # FIXME: Deal with two-tailed inference?
        self.provBundle.activity(NIIRI['inference_id_'+contrast_num], 
            other_attributes=( (PROV['type'], NIDM['Inference']), 
                               (PROV['label'] , "Inference: "+contrast_name),
                               (NIDM['hasAlternativeHypothesis'] , NIDM['OneTailedTest'])))
        self.provBundle.used(NIIRI['inference_id_'+contrast_num], NIIRI['height_threshold_id'])
        self.provBundle.used(NIIRI['inference_id_'+contrast_num], NIIRI['extent_threshold_id'])
        self.provBundle.used(NIIRI['inference_id_'+contrast_num], NIIRI['z_statistic_map_id_'+contrast_num])

        self.provBundle.wasAssociatedWith(NIIRI['inference_id_'+contrast_num], NIIRI['software_id'])

        self.provBundle.wasGeneratedBy(NIIRI['search_space_id'], NIIRI['inference_id_'+contrast_num])
        # self.provBundle.wasGeneratedBy(NIIRI['stat_image_properties_id'], NIIRI['inference_id_'+contrast_num])

    # Generate prov for a coordinate space entity 
    def create_coordinate_space(self, niftiFile):
        self.coordinateSpaceId = self.coordinateSpaceId + 1
        thresImg = nib.load(niftiFile)
        thresImgHdr = thresImg.get_header()

        numDim = len(thresImg.shape)

        # As in https://github.com/ni-/ni-dm/issues/52 (not accepted yet)
        if not self.standard_space:
            coordinateSystem = NIDM['SubjectSpace'];
        else:
            if not self.custom_standard:
                coordinateSystem = NIDM['IcbmMni152NonLinear6thGenerationCoordinateSystem'];
            else:
                coordinateSystem = NIDM['StandarizedSpace'];

        mydict = { 
            PROV['type']: NIDM['CoordinateSpace'], 
            NIDM['dimensionsInVoxels']: str(thresImg.shape).replace('(', '[').replace(')', ']'),
            NIDM['numberOfDimensions']: numDim,
            NIDM['voxelToWorldMapping']: '%s'%', '.join(str(thresImg.get_qform()).strip('()').replace('. ', '').split()).replace('[,', '[').replace('\n', ''),
            # FIXME: How to get the coordinate system? default for FSL?
            NIDM['inWorldCoordinateSystem']: coordinateSystem,           
            # FIXME: this gives mm, sec => what is wrong: FSL file, nibabel, other?
            # NIDM['voxelUnits']: '[%s]'%str(thresImgHdr.get_xyzt_units()).strip('()'),
            NIDM['voxelUnits']: "['mm', 'mm', 'mm']",
            NIDM['voxelSize']: '[%s]'%', '.join(map(str, thresImgHdr['pixdim'][1:(numDim+1)])),
            PROV['label']: "Coordinate space "+str(self.coordinateSpaceId)}

        self.provBundle.entity(NIIRI['coordinate_space_id_'+str(self.coordinateSpaceId)], other_attributes=mydict)
        return NIIRI['coordinate_space_id_'+str(self.coordinateSpaceId)]

    # Generate prov for search space entity generated by the inference activity
    def create_search_space(self, search_space_file, search_volume, resel_size_in_voxels, dlh):
        # Copy "Mask map" in export directory
        original_search_space_file = search_space_file
        search_space_file = os.path.join(self.export_dir, 'SearchSpace.nii.gz')
        shutil.copy(original_search_space_file, search_space_file)
        path, search_space_filename = os.path.split(search_space_file)
        search_space_file = os.path.join(self.export_dir,search_space_filename)   
        
        # Crate "Mask map" entity
        self.provBundle.entity(NIIRI['search_space_id'], other_attributes=( 
                (PROV['label'], "Search Space Map"), 
                (DCT['format'], NIDM['Nifti1Gzip']), 
                (PROV['type'], NIDM['SearchSpaceMap']), 
                (PROV['location'], Identifier("file://./"+search_space_filename)),
                (NIDM['originalFileName'], search_space_filename),
                (NIDM['atCoordinateSpace'], self.create_coordinate_space(search_space_file)),
                (FSL['searchVolumeInVoxels'], search_volume),
                (CRYPTO['sha512'], self.get_sha_sum(search_space_file)),
                (FSL['reselSizeInVoxels'], resel_size_in_voxels),
                (FSL['dlh'], dlh)))
        
        

    def create_excursion_set(self, excursion_set_file, stat_num, visualisation):
        # Copy "Excursion set map" in export directory
        original_excursion_set_file = excursion_set_file
        excursion_set_file = os.path.join(self.export_dir, 'ExcursionSet.nii.gz')
        shutil.copy(original_excursion_set_file, excursion_set_file)
        path, excursion_set_filename = os.path.split(excursion_set_file)
        excursion_set_file = os.path.join(self.export_dir,excursion_set_filename)   

        # Copy visualisation of excursion set in export directory
        shutil.copy(visualisation, self.export_dir)
        path, visu_filename = os.path.split(visualisation)

        # Create "Excursion set" entity
        self.provBundle.entity(NIIRI['excursion_set_id_'+str(stat_num)], other_attributes=( 
            (PROV['type'], NIDM['ExcursionSet']), 
            (DCT['format'], NIDM['Nifti1Gzip']), 
            (PROV['location'], Identifier("file://./"+excursion_set_filename)),
            (NIDM['originalFileName'], excursion_set_filename),
            (NIDM['atCoordinateSpace'], self.create_coordinate_space(excursion_set_file)),
            (PROV['label'], "Excursion Set"),
            (NIDM['visualisation'], NIIRI['excursion_set_png_id_'+str(stat_num)]),
            (CRYPTO['sha512'], self.get_sha_sum(excursion_set_file)),
            ))

        # Create "png visualisation of Excursion set" entity
        self.provBundle.entity(NIIRI['excursion_set_png_id_'+str(stat_num)], other_attributes=( 
            (PROV['type'], NIDM['Image']), 
            (PROV['location'], Identifier("file://./"+visu_filename)),
            (DCT['format'], "image/png"),
            ))
        self.provBundle.wasGeneratedBy(NIIRI['excursion_set_id_'+str(stat_num)], NIIRI['inference_id_'+str(stat_num)])

    def save_prov_to_files(self, showattributes=False):
        self.provDocument.add_bundle(self.provBundle)

        suffixName = ''
        if showattributes is False:
            suffixName = '_without_attributes'

        # jsondata = self.provBundle.get_provjson(indent=4)
        # JSONfile = open(os.path.join(self.export_dir, 'nidm.json'), 'w');
        # JSONfile.write(jsondata)
        PROVNfile = open(os.path.join(self.export_dir, 'nidm.provn'), 'w');
        PROVNfile.write(self.provDocument.get_provn(4))