'''Python implementation of NI-DM (for statistical results)

@author: Camille Maumet <c.m.j.maumet@warwick.ac.uk>
@copyright: University of Warwick 2013-2014
'''

from prov.model import ProvBundle, Namespace, ProvRecord, ProvExceptionCannotUnifyAttribute, graph, ProvEntity, Identifier
import prov.model.graph
from prov.model import PROV
import os
import numpy as np
import nibabel as nib
import hashlib
import shutil


NIDM = Namespace('nidm', "http://www.incf.org/ns/nidash/nidm#")
NIIRI = Namespace("niiri", "http://iri.nidash.org/")
CRYPTO = Namespace("crypto", "http://www.w3.org/2000/10/swap/crypto#")
FSL = Namespace("fsl", "http://www.incf.org/ns/nidash/fsl#")

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
        self.provBundle = ProvBundle()
        self.provBundle.add_namespace("neurolex", "http://neurolex.org/wiki/")
        self.provBundle.add_namespace(FSL)
        self.provBundle.add_namespace(NIDM)
        self.provBundle.add_namespace(NIIRI)
        self.provBundle.add_namespace(CRYPTO)
       
       
        # FIXME: Check one-tailed or two-tailed test and get test type from data
        # FIXME: We want to be able to add for than one inference activity for on graph -> create a function for that
        
        
        # # FIXME: is this really empty? If so, should be deleted
        # g.entity(NIIRI['stat_image_properties_id'], other_attributes=( 
        #     (PROV['type'], FSL['statisticImageProperties']), 
        #     (PROV['label'], 'Statistical image properties')))
        
    def create_software(self, feat_version):
        # Add software agent: FSL
        self.provBundle.agent(NIIRI['software_id'], other_attributes=( 
            (PROV['type'], NIDM['Fsl']), 
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

    def create_coordinate(self, coordinate_id, label_id, x, y, z, x_std, y_std, z_std):
        self.provBundle.entity(coordinate_id, other_attributes=( 
            (PROV['type'] , PROV['Location']), 
            (PROV['type'] , NIDM['Coordinate']),
            (PROV['label'] , "Coordinate "+label_id),
            # FIXME: Set coordinate system
            (NIDM['coordinateSystem'] , NIDM['mniCoordinateSystem']),
            (NIDM['coordinate1'] , x),
            (NIDM['coordinate2'] , y),
            (NIDM['coordinate3'] , z),
            (NIDM['coordinate1InUnits'] , x_std),
            (NIDM['coordinate2InUnits'] , y_std),
            (NIDM['coordinate3InUnits'] , z_std)
            ))

    def get_sha_sum(self, nifti_file):
        nifti_img = nib.load(nifti_file)
        return hashlib.sha224(nifti_img.get_data()).hexdigest()

    def create_cluster(self, *args, **kwargs):
        clusterIndex = kwargs.pop('id')
        stat_num = int(kwargs.pop('stat_num'))

        # FIXME deal with multiple contrasts
        cluster_id = clusterIndex

        self.provBundle.entity(NIIRI['cluster_000'+str(cluster_id)], other_attributes=( 
                             (PROV['type'] , NIDM['ClusterLevelStatistic']), 
                             (PROV['label'], "Cluster Level Statistic: 000"+str(cluster_id)),
                             (NIDM['clusterSizeInVoxels'], kwargs.pop('size')),
                             (NIDM['pValueFWER'], kwargs.pop('pFWER') )))
        self.provBundle.wasDerivedFrom(NIIRI['cluster_000'+str(cluster_id)], NIIRI['excursion_set_id_'+str(stat_num)])

        self.create_coordinate(NIIRI['COG_coordinate_000'+str(cluster_id)], '000'+str(cluster_id),kwargs.pop('COG1'), kwargs.pop('COG2'), kwargs.pop('COG3'), kwargs.pop('COG1_std'), kwargs.pop('COG2_std'), kwargs.pop('COG3_std'))

        self.provBundle.entity(NIIRI['center_of_gravity_'+str(cluster_id)], other_attributes=( 
                     (PROV['type'] , FSL['CenterOfGravity']), 
                     (PROV['label'], "Center of gravity "+str(cluster_id)),
                     (PROV['location'] , NIIRI['COG_coordinate_000'+str(cluster_id)]))   )
        self.provBundle.wasDerivedFrom(NIIRI['center_of_gravity_'+str(cluster_id)], NIIRI['cluster_000'+str(cluster_id)])

    def create_peak(self, *args, **kwargs):
        peakIndex = kwargs.pop('id')
        clusterIndex = kwargs.pop('cluster_id')
        stat_num = int(kwargs.pop('stat_num'))

        # FIXME: Currently assumes less than 10 clusters per contrast
        cluster_id = clusterIndex

        # FIXME: Currently assumes less than 100 peaks 
        peakUniqueId = '000'+str(clusterIndex)+'_'+str(peakIndex)

        self.create_coordinate(NIIRI['coordinate_'+str(peakUniqueId)], str(peakUniqueId), kwargs.pop('x'), kwargs.pop('y'), kwargs.pop('z'), kwargs.pop('std_x'), kwargs.pop('std_y'), kwargs.pop('std_z'))

        self.provBundle.entity(NIIRI['peak_'+str(peakUniqueId)], other_attributes=( 
            (PROV['type'] , NIDM['PeakLevelStatistic']), 
            (PROV['label'] , "Peak "+str(peakUniqueId)), 
            (NIDM['equivalentZStatistic'], kwargs.pop('equivZ')),
            (PROV['location'] , NIIRI['coordinate_'+str(peakUniqueId)]))         )
        self.provBundle.wasDerivedFrom(NIIRI['peak_'+str(peakUniqueId)], NIIRI['cluster_000'+str(cluster_id)])

    def create_model_fitting(self, residualsFile, design_matrix):
        # Create cvs file containing design matrix
        design_matrix_csv = 'design_matrix.csv'
        np.savetxt(os.path.join(self.export_dir, design_matrix_csv), np.asarray(design_matrix), delimiter=",")

        path, filename = os.path.split(residualsFile)
        self.provBundle.entity(NIIRI['residual_mean_squares_map_id'], 
            other_attributes=( (PROV['type'],NIDM['ResidualMeanSquaresMap'],), 
                               (PROV['location'], Identifier("file://./stats/"+filename) ),
                               (PROV['label'],"Residual Mean Squares Map" ),
                               (NIDM['fileName'],filename ),
                               (CRYPTO['sha'], self.get_sha_sum(residualsFile)),
                               (NIDM['coordinateSpace'], self.create_coordinate_space(residualsFile))))
        
        self.provBundle.entity(NIIRI['design_matrix_id'], 
            other_attributes=( (PROV['type'],NIDM['DesignMatrix']), 
                               (PROV['label'],"Design Matrix"), 
                               (NIDM['fileName'],design_matrix_csv ),
                               (PROV['location'], Identifier("file://./"+design_matrix_csv))))
        self.provBundle.activity(NIIRI['model_parameters_estimation_id'], other_attributes=( 
            (PROV['type'], NIDM['ModelParametersEstimation']),(PROV['label'], "Model Parameters Estimation")))
        self.provBundle.used(NIIRI['model_parameters_estimation_id'], NIIRI['design_matrix_id'])
        self.provBundle.wasAssociatedWith(NIIRI['model_parameters_estimation_id'], NIIRI['software_id'])
        
        self.provBundle.wasGeneratedBy(NIIRI['residual_mean_squares_map_id'], NIIRI['model_parameters_estimation_id'])  

    # Generate prov for contrast map
    def create_parameter_estimate(self, pe_file, pe_num):
        # Copy parameter estimate map in export directory
        shutil.copy(pe_file, self.export_dir)
        path, pe_filename = os.path.split(pe_file)
        pe_file = os.path.join(self.export_dir,pe_filename)       

        # Parameter estimate entity
        self.provBundle.entity(NIIRI['beta_map_id_'+str(pe_num)], 
            other_attributes=( (PROV['type'], NIDM['BetaMap']), 
                               (PROV['location'], Identifier("file://./"+pe_filename)),
                               (NIDM['fileName'], pe_filename), 
                               (NIDM['coordinateSpace'], self.create_coordinate_space(pe_file)),
                               (CRYPTO['sha'], self.get_sha_sum(pe_file)),
                               (PROV['label'], "Parameter estimate "+str(pe_num))))
        
        self.provBundle.wasGeneratedBy(NIIRI['beta_map_id_'+str(pe_num)], NIIRI['model_parameters_estimation_id'])  

    # Generate prov for contrast map
    def create_contrast_map(self, copeFile, varCopeFile, statFile, zStatFile, contrastName, contrast_num, dof, contrastWeights):
        # Contrast id entity
        # FIXME: Get contrast weights
        self.provBundle.entity(NIIRI['contrast_id_'+contrast_num], 
            other_attributes=( (PROV['type'], NIDM['TContrast']), 
                               (PROV['label'], "T Contrast: "+contrastName), 
                               (NIDM['contrastName'], contrastName),
                               (NIDM['contrastWeights'], contrastWeights)))

        # Create related activities
        self.provBundle.activity(NIIRI['contrast_estimation_id_'+contrast_num], other_attributes=( 
            (PROV['type'], NIDM['ContrastEstimation']),
            (PROV['label'], "Contrast estimation: "+contrastName)))

        # Contrast Map entity
        path, filename = os.path.split(copeFile)
        self.provBundle.entity('niiri:'+'contrast_map_id_'+contrast_num, other_attributes=( 
            (PROV['type'], NIDM['ContrastMap']), 
            (NIDM['coordinateSpace'], self.create_coordinate_space(copeFile)),
            (PROV['location'], Identifier("file://./stats/"+filename)),
            (NIDM['fileName'], filename),
            (NIDM['contrastName'], contrastName),
            (CRYPTO['sha'], self.get_sha_sum(copeFile)),
            (PROV['label'], "Contrast Map: "+contrastName)))
        
        self.provBundle.wasGeneratedBy(NIIRI['contrast_map_id_'+contrast_num], NIIRI['contrast_estimation_id_'+contrast_num])
        self.provBundle.wasAssociatedWith(NIIRI['contrast_estimation_id_'+contrast_num], NIIRI['software_id'])

        # Contrast Variance Map entity
        path, filename = os.path.split(varCopeFile)

        self.provBundle.entity('niiri:'+'contrast_variance_map_id_'+contrast_num, other_attributes=( 
            (PROV['type'], FSL['VarCope']), 
            (NIDM['coordinateSpace'], self.create_coordinate_space(varCopeFile)),
            (PROV['location'], Identifier("file://./stats/"+filename)),
            (CRYPTO['sha'], self.get_sha_sum(varCopeFile)),
            (NIDM['fileName'], filename),
            (PROV['label'], "Contrast Variance Map "+contrast_num)))
        
        self.provBundle.wasGeneratedBy(NIIRI['contrast_variance_map_id_'+contrast_num], NIIRI['contrast_estimation_id_'+contrast_num])

        # Create standard error map from contrast variance map
        varCopeImg = nib.load(varCopeFile)
        contrastVariance = varCopeImg.get_data()

        standardErrorImg = nib.Nifti1Image(np.sqrt(contrastVariance), varCopeImg.get_qform())
        standardErrorFile = varCopeFile.replace('var', 'sqrt_var')
        nib.save(standardErrorImg, standardErrorFile)

        path, filename = os.path.split(standardErrorFile)
        self.provBundle.entity('niiri:'+'contrast_standard_error_map_id_'+contrast_num, other_attributes=( 
            (PROV['type'], NIDM['ContrastStandardErrorMap']), 
            (NIDM['coordinateSpace'], self.create_coordinate_space(standardErrorFile)),
            (PROV['location'], Identifier("file://./stats/"+filename)),
            (CRYPTO['sha'], self.get_sha_sum(standardErrorFile)),
            (NIDM['fileName'], filename),
            (PROV['label'], "Contrast Standard Error Map")))
        
        self.provBundle.wasDerivedFrom(NIIRI['contrast_standard_error_map_id_'+contrast_num], NIIRI['contrast_variance_map_id_'+contrast_num])

        
        # FIXME: Remove TODOs
        # FIXME: Add sha sum

        # Z-Statistical Map entity
        path, filename = os.path.split(zStatFile)
        self.provBundle.entity(NIIRI['z_statistical_map_id_'+contrast_num ],
            other_attributes=(  (PROV['type'], FSL['ZStatisticalMap']), 
                                (PROV['label'], "Z-statistical Map: "+contrastName) ,
                                (PROV['location'], Identifier("file://./stats/"+filename)),
                                (NIDM['contrastName'], contrastName),
                                (NIDM['fileName'], filename),
                                (CRYPTO['sha'], self.get_sha_sum(zStatFile)),
                                (NIDM['coordinateSpace'], self.create_coordinate_space(zStatFile)),
                                ) )

        # Statistical Map entity
        path, filename = os.path.split(statFile)
        # FIXME: Deal with other than t-contrast maps: dof
        self.provBundle.entity(NIIRI['statistical_map_id_'+contrast_num ],
            other_attributes=(  (PROV['type'], NIDM['TStatisticalMap']), 
                                (PROV['label'], "Statistical Map: "+contrastName) ,
                                (PROV['location'], Identifier("file://./stats/"+filename)),
                                (NIDM['fileName'], filename),
                                (NIDM['contrastName'], contrastName),
                                (NIDM['errorDegreesOfFreedom'], dof),
                                (NIDM['effectDegreesOfFreedom'], 1.0),
                                (CRYPTO['sha'], self.get_sha_sum(statFile)),
                                (NIDM['coordinateSpace'], self.create_coordinate_space(statFile)),
                                ) )
        
        self.provBundle.wasGeneratedBy(NIIRI['statistical_map_id_'+contrast_num], NIIRI['contrast_estimation_id_'+contrast_num])
               
        self.provBundle.wasGeneratedBy(NIIRI['z_statistical_map_id_'+contrast_num], NIIRI['contrast_estimation_id_'+contrast_num])
        self.provBundle.used(NIIRI['contrast_estimation_id_'+contrast_num], NIIRI['residual_mean_squares_map_id'])
        self.provBundle.used(NIIRI['contrast_estimation_id_'+contrast_num], NIIRI['design_matrix_id'])
        self.provBundle.used(NIIRI['contrast_estimation_id_'+contrast_num], NIIRI['contrast_id_'+contrast_num])


        # In FSL we have a single thresholding (extent, height) applied to all contrasts 
        # FIXME: Deal with two-tailed inference?
        self.provBundle.activity(NIIRI['inference_id_'+contrast_num], 
            other_attributes=( (PROV['type'], NIDM['InferenceOneTailed']), 
                               (PROV['label'] , "Inference: "+contrastName)))
        self.provBundle.used(NIIRI['inference_id_'+contrast_num], NIIRI['height_threshold_id'])
        self.provBundle.used(NIIRI['inference_id_'+contrast_num], NIIRI['extent_threshold_id'])
        self.provBundle.used(NIIRI['inference_id_'+contrast_num], NIIRI['z_statistical_map_id_'+contrast_num])

        self.provBundle.wasAssociatedWith(NIIRI['inference_id_'+contrast_num], NIIRI['software_id'])

        self.provBundle.wasGeneratedBy(NIIRI['search_space_id'], NIIRI['inference_id_'+contrast_num])
        # self.provBundle.wasGeneratedBy(NIIRI['stat_image_properties_id'], NIIRI['inference_id_'+contrast_num])

    # Generate prov for a coordinate space entity 
    def create_coordinate_space(self, niftiFile):
        self.coordinateSpaceId = self.coordinateSpaceId + 1
        thresImg = nib.load(niftiFile)
        thresImgHdr = thresImg.get_header()

        numDim = len(thresImg.shape)

        mydict = { 
            PROV['type']: NIDM['CoordinateSpace'], 
            NIDM['dimensions']: str(thresImg.shape).replace('(', '[').replace(')', ']'),
            NIDM['numberOfDimensions']: numDim,
            NIDM['voxelToWorldMapping']: '%s'%', '.join(str(thresImg.get_qform()).strip('()').replace('. ', '').split()).replace('[,', '[').replace('\n', ''),
            # FIXME: How to get the coordinate system? default for FSL?
            NIDM['coordinateSystem']: NIDM['mniCoordinateSystem'],           
            # FIXME: this gives mm, sec => what is wrong: FSL file, nibabel, other?
            # NIDM['voxelUnits']: '[%s]'%str(thresImgHdr.get_xyzt_units()).strip('()'),
            NIDM['voxelUnits']: "['mm', 'mm', 'mm']",
            NIDM['voxelSize']: '[%s]'%', '.join(map(str, thresImgHdr['pixdim'][1:(numDim+1)])),
            PROV['label']: "Coordinate space "+str(self.coordinateSpaceId)}

        self.provBundle.entity(NIIRI['coordinate_space_id_'+str(self.coordinateSpaceId)], other_attributes=mydict)
        return NIIRI['coordinate_space_id_'+str(self.coordinateSpaceId)]

    # Generate prov for search space entity generated by the inference activity
    def create_search_space(self, search_space_file, search_volume, resel_size_in_voxels, dlh):
        path, filename = os.path.split(search_space_file)

        self.provBundle.entity(NIIRI['search_space_id'], other_attributes=( 
                (PROV['label'], "Search Space"), 
                (PROV['type'], NIDM['Mask']), 
                (PROV['location'], Identifier("file://./"+filename)),
                (NIDM['fileName'], filename),
                (NIDM['coordinateSpace'], self.create_coordinate_space(search_space_file)),
                (NIDM['searchVolumeInVoxels'], search_volume),
                (CRYPTO['sha'], self.get_sha_sum(search_space_file)),
                (FSL['reselSizeInVoxels'], resel_size_in_voxels),
                (FSL['dlh'], dlh)))
        
        

    def create_excursion_set(self, excursion_set_file, stat_num, visualisation):
        zFileImg = excursion_set_file
        path, filename = os.path.split(zFileImg)

        # Copy visualisation of excursion set in export directory
        shutil.copy(visualisation, self.export_dir)
        path, visu_filename = os.path.split(visualisation)

        self.provBundle.entity(NIIRI['excursion_set_id_'+str(stat_num)], other_attributes=( 
            (PROV['type'], NIDM['ExcursionSet']), 
            (PROV['location'], Identifier("file://./"+filename)),
            (NIDM['fileName'], filename),
            (NIDM['coordinateSpace'], self.create_coordinate_space(zFileImg)),
            (PROV['label'], "Excursion Set"),
            (NIDM['visualisation'], Identifier("file://./"+visu_filename)),
            (CRYPTO['sha'], self.get_sha_sum(excursion_set_file)),
            ))
        self.provBundle.wasGeneratedBy(NIIRI['excursion_set_id_'+str(stat_num)], NIIRI['inference_id_'+str(stat_num)])
        

    def save_prov_to_files(self, showattributes=False):
        suffixName = ''
        if showattributes is False:
            suffixName = '_without_attributes'

        jsondata = self.provBundle.get_provjson(indent=4)
        JSONfile = open('./FSL_example.json', 'w');
        JSONfile.write(jsondata)
        PROVNfile = open('./FSL_example.provn', 'w');
        PROVNfile.write(self.provBundle.get_provn(4))

        dot = graph.prov_to_dot(self.provBundle, use_labels=True, show_element_attributes=showattributes)
        dot.set_dpi(200)
        dot.write_png('./FSL_example'+suffixName+'.png')

