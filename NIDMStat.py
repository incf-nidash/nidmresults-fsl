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
    	g = ProvBundle()
        g.add_namespace("neurolex", "http://neurolex.org/wiki/")
        g.add_namespace("spm", "http://www.fil.ion.ucl.ac.uk/spm/ns/")
        g.add_namespace("nidm", "http://nidm.nidash.org/")
        g.add_namespace("niiri", "http://iri.nidash.org/")
        g.add_namespace("crypto", "http://www.w3.org/2000/10/swap/crypto#")

        g.entity('niiri:statistical_map_id' ,
            other_attributes=(  ('prov:type', 'nidm:statisticalMap'), 
                                ('prov:label', "Statistical Map: TODO") ,
                                ('prov:contrastName', "TODO") ,
                                ('prov:location', "file:///path/to/TODO.img"),
                                ('prov:fileName', "TODO.img"),
                                ('prov:statisticType', 'nidm:tStatisticTODO'),
                                ('prov:errorDegreesOfFreedom', 'TODO'),
                                ('prov:effectDegreesOfFreedom', 'TODO')
                                ) )
        g.entity('niiri:contrast_map_id', other_attributes=(  ('prov:type', 'nidm:contrastMap'), ('prov:type', 'nidm:contrastMap2')))
        g.wasDerivedFrom('niiri:statistical_map_id', 'niiri:contrast_map_id')

        g.entity('niiri:contrast_standard_error_map_id', other_attributes=( ('prov:type' , 'nidm:contrastStandardErrorMap'), ('prov:type' , 'nidm:contrastStandardErrorMap2')))

        g.wasDerivedFrom('niiri:statistical_map_id', 'niiri:contrast_standard_error_map_id')
        g.entity('niiri:residual_mean_squares_map_id', other_attributes=( ('prov:type','nidm:residualMeanSquaresMap',), ('prov:location',"file:///path/to/ResMS.img" )))
        g.wasDerivedFrom('niiri:contrast_standard_error_map_id', 'niiri:residual_mean_squares_map_id')
        g.entity('niiri:design_matrix_id', other_attributes=( ('prov:type','nidm:designMatrix',), ('prov:location', "file:///path/to/design_matrix.csv")))
        g.entity('niiri:contrast_id', other_attributes=( ('prov:type', 'nidm:contrast'), ('nidm:contrastName',"listening &gt; rest")))
        g.wasDerivedFrom('niiri:contrast_map_id', 'niiri:design_matrix_id')
        g.wasDerivedFrom('niiri:contrast_map_id', 'niiri:contrast_id')
        g.wasDerivedFrom('niiri:contrast_standard_error_map_id', 'niiri:design_matrix_id')
        g.wasDerivedFrom('niiri:contrast_standard_error_map_id', 'niiri:contrast_id')
        g.wasDerivedFrom('niiri:residual_mean_squares_map_id', 'niiri:design_matrix_id')

       
        # FIXME: Check one-tailed or two-tailed test and get test type from data
        g.activity('niiri:inference_id', other_attributes=( ('prov:type', 'fsl:inference'), ('prov:label' , "Inference"), ('prov:statisticalTest' , 'nidm:oneTailedTtest')))
        g.used('niiri:inference_id', 'niiri:statistical_map_id')
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
    	# FIXME: update TODOs
        heightThreshAllFields = {
            'prov:type': 'nidm:heightThreshold', 'prov:label': "Height Threshold: TODO",
            'prov:userSpecifiedThresholdType': "TODO" , 'prov:value': kwargs.pop('voxelThreshold'),
            'nidm:pValueUncorrected': kwargs.pop('voxelPUncorr'), 'fsl:pValueGRF': kwargs.pop('voxelPCorr')
            }
        self.provBundle.entity('niiri:height_threshold_id', other_attributes=dict((k,v) for k,v in heightThreshAllFields.iteritems() if v is not None))
        exentThreshAllFields = {
            'prov:type': 'nidm:extentThreshold', 'prov:label': "Extent Threshold: TODO", 'nidm:clusterSizeInVoxels': kwargs.pop('extent'),
            'nidm:pValueUncorrected': kwargs.pop('extentPUncorr'), 'fsl:pValueGRF': kwargs.pop('extentPCorr')
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
                             ('prov:type' , 'spm:clusterLevelStatistic'), 
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
                     ('prov:atLocation' , 'niiri:coordinate_'+str(clusterId)+"00"))   )
        self.provBundle.wasDerivedFrom('niiri:centerOfGravity_'+str(clusterId), 'niiri:cluster_'+str(clusterId))

    def create_peak(self, *args, **kwargs):
    	peakIndex = kwargs.pop('id')
    	self.provBundle.entity('niiri:coordinate_'+str(peakIndex), other_attributes=( 
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
        self.provBundle.entity('niiri:peak_'+str(peakIndex), other_attributes=( 
            ('prov:type' , 'spm:peakLevelStatistic'), 
            ('nidm:equivalentZStatistic', str(kwargs.pop('equivZ'))),
            ('prov:atLocation' , 'niiri:coordinate_'+str(peakIndex)))         )
        self.provBundle.wasDerivedFrom('niiri:peak_'+str(peakIndex), 'niiri:cluster_'+str(kwargs.pop('clusterId')))

    def create_mask_info(self, *args, **kwargs):
    	self.provBundle.entity('niiri:search_space_id', other_attributes=( 
                ('prov:label', "Search space"), 
                ('prov:type', 'nidm:mask'), 
                ('prov:location', "file:///path/to/mask.nii.gz")))

    def create_excursion_set(self, *args, **kwargs):
    	zFileImg = kwargs.pop('zFileImg')
    	# path, zFilename = os.path.split(kwargs.pop('zFileImg'));
    	thresImg = nib.load(zFileImg)
        thresImgHdr = thresImg.get_header()

        numDim = len(thresImg.shape)

    	self.provBundle.entity('niiri:excursion_set_id', other_attributes=( 
            ('prov:type', 'fsl:excursionSet'), 
            ('prov:location', "file:./"+zFileImg),
            ('nidm:fileName', zFileImg),
            ('nidm:voxelToWorldMapping', str(thresImg.get_qform()).replace('(', '[').replace(')', ']')),
            ('nidm:numberOfDimensions', str(numDim)),
            ('nidm:dimensions', str(thresImg.shape).replace('(', '[').replace(')', ']')),
            ('nidm:voxelUnits', str(thresImgHdr.get_xyzt_units()).replace('(', '[').replace(')', ']')),
            ('nidm:voxelSize', str(thresImgHdr['pixdim'][1:(numDim+1)])),
            ('prov:label', "Excursion set"),
            ))
        self.provBundle.wasGeneratedBy('niiri:excursion_set_id', 'niiri:inference_id')

    def save_prov_to_files(self):
    	jsondata = self.provBundle.get_provjson(indent=4)
        JSONfile = open('./FSL_example.json', 'w');
        JSONfile.write(jsondata)
        PROVNfile = open('./FSL_example.provn', 'w');
        PROVNfile.write(self.provBundle.get_provn(4))

        dot = graph.prov_to_dot(self.provBundle, use_labels=True)
        dot.set_dpi(120)
        dot.write_png('./FSL_example.png')

