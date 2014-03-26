'''Python implementation of NI-DM (for statistical results)

@author: Camille Maumet <c.m.j.maumet@warwick.ac.uk>
@copyright: University of Warwick 2013-2014
'''

from prov.model import ProvBundle, ProvRecord, ProvExceptionCannotUnifyAttribute, graph, ProvEntity
import prov.model.graph
import os
import numpy as np
import nibabel as nib

class NIDMStat():

    def __init__(self, *args, **kwargs):
        # Keep track of the number of coordinateSpace entity already generated
        # FIXME: Merge coordinateSpace entities if identical?
        # FIXME: Use actual URIs instead
        self.coordinateSpaceId = 1

        g = ProvBundle()
        g.add_namespace("neurolex", "http://neurolex.org/wiki/")
        g.add_namespace("fsl", "http://fsl.fmrib.ox.ac.uk/")
        g.add_namespace("nidm", "http://nidm.nidash.org/")
        g.add_namespace("niiri", "http://iri.nidash.org/")
        g.add_namespace("crypto", "http://www.w3.org/2000/10/swap/crypto#")

        

        
       
        # FIXME: Check one-tailed or two-tailed test and get test type from data
        # FIXME: We want to be able to add for than one inference activity for on graph -> create a function for that
        g.activity('niiri:inference_id', other_attributes=( ('prov:type', 'fsl:inference'), ('prov:label' , "Inference"), ('prov:statisticalTest' , 'nidm:oneTailedTtest')))
        
        g.wasGeneratedBy('niiri:search_space_id', 'niiri:inference_id')
        g.used('niiri:inference_id', 'niiri:height_threshold_id')
        g.used('niiri:inference_id', 'niiri:extent_threshold_id')
        
        # FIXME: is this really empty? If so, should be deleted
        g.entity('niiri:stat_image_properties_id', other_attributes=( 
            ('prov:type', 'fsl:statisticImageProperties'), 
            ('prov:label', 'Statistical image properties')))
        g.wasGeneratedBy('niiri:stat_image_properties_id', 'niiri:inference_id')
        
        self.provBundle = g
        
    def create_thresholds(self, *args, **kwargs):
        voxelThreshold = kwargs.pop('voxelThreshold')
        voxelPUncorr = kwargs.pop('voxelPUncorr')
        voxelPCorr = kwargs.pop('voxelPCorr')
        threshDesc = ""
        if not voxelThreshold is None:
            threshDesc = ">"+str(voxelThreshold)
        elif not voxelPUncorr is None:
            threshDesc = "p<"+str(voxelPUncorr)+" uncorr."
        elif not voxelPCorr is None:
            threshDesc = "p<"+str(voxelPCorr)+" (GRF)"
        heightThreshAllFields = {
            'prov:type': 'nidm:heightThreshold', 'prov:label': "Height Threshold: "+threshDesc,
            'prov:userSpecifiedThresholdType': threshDesc , 'prov:value': voxelThreshold,
            'nidm:pValueUncorrected': voxelPUncorr, 'fsl:pValueGRF': voxelPCorr
            }
        self.provBundle.entity('niiri:height_threshold_id', other_attributes=dict((k,v) for k,v in heightThreshAllFields.iteritems() if v is not None))

        extent = kwargs.pop('extent')
        extentPUncorr = kwargs.pop('extentPUncorr')
        extentPCorr = kwargs.pop('extentPCorr')
        threshDesc = ""
        if not extent is None:
            threshDesc = "k>"+str(extent)
        elif not extentPUncorr is None:
            threshDesc = "p<"+str(extentPUncorr)+" uncorr."
        elif not extentPCorr is None:
            threshDesc = "p<"+str(extentPCorr)+" corr."
        exentThreshAllFields = {
            'prov:type': 'nidm:extentThreshold', 'prov:label': "Extent Threshold: "+threshDesc, 'nidm:clusterSizeInVoxels': extent,
            'nidm:pValueUncorrected': extentPUncorr, 'fsl:pValueGRF': extentPCorr
        }
        self.provBundle.entity('niiri:extent_threshold_id', other_attributes=dict((k,v) for k,v in exentThreshAllFields.iteritems() if v is not None))


        self.provBundle.agent('niiri:software_id', other_attributes=( 
            ('prov:type', 'nidm:fsl'), 
            ('prov:type','prov:SoftwareAgent'),
            ('prov:label','FSL'),
            ('nidm:softwareVersion','FSL') ))

        self.provBundle.wasAssociatedWith('niiri:inference_id', 'niiri:software_id')

    def create_cluster(self, *args, **kwargs):
        clusterId = kwargs.pop('id')
        self.provBundle.entity('niiri:cluster_'+str(clusterId), other_attributes=( 
                             ('prov:type' , 'fsl:clusterLevelStatistic'), 
                             ('prov:label', "Cluster Level Statistic: "+str(clusterId)),
                             ('fsl:clusterSizeInVoxels', str(kwargs.pop('size'))),
                             ('fsl:pGRF', str(kwargs.pop('pGRF')) )))
        self.provBundle.wasDerivedFrom('niiri:cluster_'+str(clusterId), 'niiri:excursion_set_id')

        self.provBundle.entity('niiri:coordinate_'+str(clusterId)+"00", other_attributes=( 
            ('prov:type' , 'prov:location'), 
            ('prov:type' , 'nidm:coordinate'),
            ('nidm:coordinateSystem' , 'nidm:mni'),
            ('nidm:coordinate1' , str(kwargs.pop('COG1'))),
            ('nidm:coordinate2' , str(kwargs.pop('COG2'))),
            ('nidm:coordinate3' , str(kwargs.pop('COG3'))),
            ('nidm:coordinateInUnits1' , str(kwargs.pop('COG1_std'))),
            ('nidm:coordinateInUnits2' , str(kwargs.pop('COG2_std'))),
            ('nidm:coordinateInUnits3' , str(kwargs.pop('COG3_std')))
            ))
        self.provBundle.entity('niiri:centerOfGravity_'+str(clusterId), other_attributes=( 
                     ('prov:type' , 'fsl:centerOfGravity'), 
                     ('prov:label', "Center of Gravity: "+str(clusterId)),
                     ('prov:location' , 'niiri:coordinate_'+str(clusterId)+"00"))   )
        self.provBundle.wasDerivedFrom('niiri:centerOfGravity_'+str(clusterId), 'niiri:cluster_'+str(clusterId))

    def create_peak(self, *args, **kwargs):
        peakIndex = kwargs.pop('id')
        clusterIndex = kwargs.pop('clusterId');

        # FIXME: Currently assumes less than 100 peaks per cluster
        peakUniqueId = clusterIndex*100+peakIndex

        self.provBundle.entity('niiri:coordinate_'+str(peakUniqueId), other_attributes=( 
                    ('prov:type' , 'prov:location'), 
                    ('prov:type' , 'nidm:coordinate'),
                    ('prov:label' , "Coordinate "+str(peakIndex)),
                    ('nidm:coordinate1' , kwargs.pop('x')),
                    ('nidm:coordinate2' , kwargs.pop('y')),
                    ('nidm:coordinate3' , kwargs.pop('z')),
                    ('nidm:coordinateInUnits1' , str(kwargs.pop('std_x'))),
                    ('nidm:coordinateInUnits2' , str(kwargs.pop('std_y'))),
                    ('nidm:coordinateInUnits3' , str(kwargs.pop('std_z')))
                    ))
        self.provBundle.entity('niiri:peak_'+str(peakUniqueId), other_attributes=( 
            ('prov:type' , 'fsl:peakLevelStatistic'), 
            ('nidm:equivalentZStatistic', str(kwargs.pop('equivZ'))),
            ('prov:location' , 'niiri:coordinate_'+str(peakUniqueId)))         )
        self.provBundle.wasDerivedFrom('niiri:peak_'+str(peakUniqueId), 'niiri:cluster_'+str(clusterIndex))

    def create_model_fitting(self, residualsFile):
        # FIXME: Add crypto sha

        path, filename = os.path.split(residualsFile)
        self.provBundle.entity('niiri:residual_mean_squares_map_id', 
            other_attributes=( ('prov:type','nidm:residualMeanSquaresMap',), 
                               ('prov:location',"file:///path/to/"+filename ),
                               ('prov:label',"Residual Mean Squares Map" ),
                               ('prov:fileName',filename ),
                               ('crypto:sha',"TODO" ),
                               ('nidm:coordinateSpace', 'coordinate_space_id_'+str(self.coordinateSpaceId))))
        self.create_coordinate_space(residualsFile)
        self.provBundle.entity('niiri:design_matrix_id', 
            other_attributes=( ('prov:type','nidm:designMatrix',), 
                               ('prov:location', "file:///path/to/design_matrix.csv")))


    # Generate prov for contrast map
    def create_contrast_map(self, copeFile, varCopeFile, statFile, contrastName, dof):
        # Contrast id entity
        # FIXME: Get contrast weights
        self.provBundle.entity('niiri:contrast_id_'+contrastName, 
            other_attributes=( ('prov:type', 'nidm:contrast'), 
                               ('nidm:contrastName', contrastName),
                               ('nidm:contrastWeights', "TODO")))

        # Create related activities
        self.provBundle.activity('niiri:contrast_estimation_id'+contrastName, other_attributes=( 
            ('prov:type', 'fsl:contrast'),
            ('prov:label', "Contrast estimation: "+contrastName)))

        self.provBundle.activity('niiri:model_fitting_id', other_attributes=( 
            ('prov:type', 'fsl:estimation'),('prov:label', "Model fitting")))
        self.provBundle.wasAssociatedWith('niiri:model_fitting_id', 'niiri:software_id')
        self.provBundle.used('niiri:model_fitting_id', 'niiri:design_matrix_id')
        self.provBundle.wasGeneratedBy('niiri:residual_mean_squares_map_id', 'niiri:model_fitting_id')

        path, filename = os.path.split(copeFile)

        self.provBundle.entity('niiri:'+'contrast_map_id_'+contrastName, other_attributes=( 
            ('prov:type', 'nidm:contrastMap'), 
            ('nidm:coordinateSpace', 'coordinate_space_id_'+str(self.coordinateSpaceId)),
            ('prov:location', "file://./"+filename),
            ('nidm:fileName', filename),
            ('nidm:contrastName', contrastName),
            ('prov:label', "Contrast map: "+contrastName)))
        self.create_coordinate_space(copeFile)
        self.provBundle.wasGeneratedBy('niiri:contrast_map_id_'+contrastName, 'niiri:contrast_estimation_id'+contrastName)
        self.provBundle.wasAssociatedWith('niiri:contrast_estimation_id'+contrastName, 'niiri:software_id')

        path, filename = os.path.split(varCopeFile)
        # FIXME: Standard error or contrast variance...
        self.provBundle.entity('niiri:'+'contrast_standard_error_map_id_'+contrastName, other_attributes=( 
            ('prov:type', 'nidm:contrastStandardErrorMap'), 
            ('nidm:coordinateSpace', 'coordinate_space_id_'+str(self.coordinateSpaceId)),
            ('prov:location', "file://./"+filename),
            ('nidm:fileName', filename),
            ('prov:label', "Contrast variance map")))
        self.create_coordinate_space(varCopeFile)
        self.provBundle.wasGeneratedBy('niiri:contrast_standard_error_map_id_'+contrastName, 'niiri:contrast_estimation_id'+contrastName)

        path, filename = os.path.split(statFile)
        # FIXME: Remove TODOs

        # Specify a different URI for each contrast (for now just by adding contrastName in URI)
        self.provBundle.entity('niiri:statistical_map_id_'+contrastName ,
            other_attributes=(  ('prov:type', 'nidm:statisticalMap'), 
                                ('prov:label', "Statistical Map: "+contrastName) ,
                                ('prov:location', "file://./"+filename),
                                ('prov:fileName', filename),
                                ('prov:statisticType', 'nidm:tStatisticTODO'),
                                ('prov:errorDegreesOfFreedom', str(dof)),
                                ('prov:effectDegreesOfFreedom', 'TODO'),
                                ('nidm:coordinateSpace', 'coordinate_space_id_'+str(self.coordinateSpaceId)),
                                ) )
        self.create_coordinate_space(statFile)
        self.provBundle.wasGeneratedBy('niiri:statistical_map_id_'+contrastName, 'niiri:contrast_estimation_id'+contrastName)
        self.provBundle.used('niiri:inference_id', 'niiri:statistical_map_id_'+contrastName)
               
        self.provBundle.used('niiri:contrast_estimation_id'+contrastName, 'niiri:residual_mean_squares_map_id')
        self.provBundle.used('niiri:contrast_estimation_id'+contrastName, 'niiri:design_matrix_id')
        self.provBundle.used('niiri:contrast_estimation_id'+contrastName, 'niiri:contrast_id_'+contrastName)

        self.provBundle.used('niiri:inference_id', 'niiri:statistical_map_id_'+contrastName)


    # Generate prov for a coordinate space entity 
    def create_coordinate_space(self, niftiFile):
        thresImg = nib.load(niftiFile)
        thresImgHdr = thresImg.get_header()

        numDim = len(thresImg.shape)

        self.provBundle.entity('niiri:'+'coordinate_space_id_'+str(self.coordinateSpaceId), other_attributes=( 
            ('prov:type', 'nidm:coordinateSpace'), 
            ('nidm:dimensions', str(thresImg.shape).replace('(', '[').replace(')', ']')),
            ('nidm:numberOfDimensions', str(numDim)),
            ('nidm:voxelToWorldMapping', str(thresImg.get_qform()).replace('(', '[').replace(')', ']')),  
            # FIXME: How to get the coordinate system? default for FSL?
            ('nidm:coordinateSystem', "TO DO"),           
            ('nidm:voxelUnits', str(thresImgHdr.get_xyzt_units()).replace('(', '[').replace(')', ']')),
            ('nidm:voxelSize', str(thresImgHdr['pixdim'][1:(numDim+1)])),
            ('prov:label', "Coordinate space "+str(self.coordinateSpaceId))))
        self.coordinateSpaceId = self.coordinateSpaceId + 1

    # Generate prov for search space entity generated by the inference activity
    def create_search_space(self, searchSpaceFile, searchVolume, reselSizeInVoxels):
        path, filename = os.path.split(searchSpaceFile)

        self.provBundle.entity('niiri:search_space_id', other_attributes=( 
                ('prov:label', "Search space"), 
                ('prov:type', 'nidm:mask'), 
                ('prov:location', "./"+filename),
                ('nidm:coordinateSpace', 'coordinate_space_id_'+str(self.coordinateSpaceId)),
                ('nidm:searchVolumeInVoxels', str(searchVolume)),
                ('fsl:reselSizeInVoxels', str(reselSizeInVoxels))))
        self.create_coordinate_space(searchSpaceFile)
        

    def create_excursion_set(self, *args, **kwargs):
        zFileImg = kwargs.pop('zFileImg')
        path, filename = os.path.split(zFileImg)

        self.provBundle.entity('niiri:excursion_set_id', other_attributes=( 
            ('prov:type', 'fsl:excursionSet'), 
            ('prov:location', "file:./"+filename),
            ('nidm:fileName', filename),
            ('nidm:coordinateSpace', 'coordinate_space_id_'+str(self.coordinateSpaceId)),
            ('prov:label', "Excursion set"),
            ))
        self.provBundle.wasGeneratedBy('niiri:excursion_set_id', 'niiri:inference_id')
        self.create_coordinate_space(zFileImg)

    def save_prov_to_files(self):
        jsondata = self.provBundle.get_provjson(indent=4)
        JSONfile = open('./FSL_example.json', 'w');
        JSONfile.write(jsondata)
        PROVNfile = open('./FSL_example.provn', 'w');
        PROVNfile.write(self.provBundle.get_provn(4))

        dot = graph.prov_to_dot(self.provBundle, use_labels=True)
        dot.set_dpi(120)
        dot.write_png('./FSL_example.png')

