from prov.model import Identifier
import os
from constants import *
import nibabel as nib
import shutil
from generic import *

class Inference(NIDMObject):
    def __init__(self, inference, height_thresh, extent_thresh, peak_criteria, cluster_criteria,
        display_mask, excursion_set, clusters, search_space, software_id):
        super(Inference, self).__init__()
        self.excursion_set = excursion_set
        self.inference_act = inference
        self.height_thresh = height_thresh
        self.extent_thresh = extent_thresh
        self.clusters = clusters
        self.software_id = software_id
        self.peak_criteria = peak_criteria
        self.cluster_criteria = cluster_criteria
        self.display_mask = display_mask
        self.search_space = search_space

    def export(self):
        # Excursion set
        self.p.update(self.excursion_set.export())

        # Height threshold
        self.p.update(self.height_thresh.export())

        # Extent threshold
        self.p.update(self.extent_thresh.export())

        # Inference activity
        self.p.update(self.inference_act.export())

        # Peak Definition
        self.p.update(self.peak_criteria.export())
        self.p.used(self.inference_act.id, self.peak_criteria.id)

        # Display Mask
        self.p.update(self.display_mask.export())
        self.p.used(self.inference_act.id, self.display_mask.id)

        # Search Space
        self.p.update(self.search_space.export())
        self.p.wasGeneratedBy(self.search_space.id, self.inference_act.id)

        # Cluster Definition
        self.p.update(self.cluster_criteria.export())
        self.p.used(self.inference_act.id, self.cluster_criteria.id)

        # Clusters and peaks
        for cluster in self.clusters:
            self.p.update(cluster.export())
            self.p.wasDerivedFrom(cluster.id, self.excursion_set.id)

        self.p.wasGeneratedBy(self.excursion_set.id, self.inference_act.id)

        self.p.wasAssociatedWith(self.inference_act.id, self.software_id)
        # self.p.wasGeneratedBy(NIIRI['search_space_id'], self.inference_act.id)
        self.p.used(self.inference_act.id, self.height_thresh.id)
        self.p.used(self.inference_act.id, self.extent_thresh.id)
        # self.p.used(self.inference_act.id, NIIRI['z_statistic_map_id_'+contrast_num])
        # self.p.used(self.inference_act.id, NIIRI['mask_id_1'])  

        return self.p

class InferenceActivity(NIDMObject):

    def __init__(self, contrast_num, contrast_name):
        super(InferenceActivity, self).__init__()
        self.id = NIIRI['inference_id_'+contrast_num]
        self.contrast_name = contrast_name

    def export(self):
        # In FSL we have a single thresholding (extent, height) applied to all contrasts 
        # FIXME: Deal with two-tailed inference?
        self.p.activity(self.id, 
            other_attributes=( (PROV['type'], NIDM['Inference']), 
                               (PROV['label'] , "Inference: "+self.contrast_name),
                               (NIDM['hasAlternativeHypothesis'] , NIDM['OneTailedTest'])))
        return self.p


class ExcursionSet(NIDMObject):
    def __init__(self, filename, stat_num, visualisation, coordinate_system, coordinate_space_id, export_dir):
        super(ExcursionSet, self).__init__(export_dir)
        self.num = stat_num
        self.file = filename
        self.id = NIIRI['contrast_estimation_id_'+self.num]
        self.visualisation = visualisation
        self.id = NIIRI['excursion_set_id_'+str(stat_num)]
        self.visu = Visualisation(visualisation, stat_num, export_dir)

        self.coord_space = CoordinateSpace(coordinate_system, coordinate_space_id, filename)

    def export(self):
        self.p.update(self.coord_space.export())

        self.p.update(self.visu.export())

        # Copy "Excursion set map" in export directory
        excursion_set_original_file = self.file
        excursion_set_file = os.path.join(self.export_dir, 'ExcursionSet.nii.gz')
        excursion_set_original_filename, excursion_set_filename = self.copy_nifti(excursion_set_original_file, excursion_set_file)

        # Create "Excursion set" entity
        self.p.entity(self.id, other_attributes=( 
            (PROV['type'], NIDM['ExcursionSet']), 
            (DCT['format'], "image/nifti"), 
            (PROV['location'], Identifier("file://./"+excursion_set_filename)),
            (NIDM['filename'], excursion_set_original_filename),
            (NIDM['filename'], excursion_set_filename),
            (NIDM['inCoordinateSpace'], self.coord_space.id),
            (PROV['label'], "Excursion Set"),
            (NIDM['visualisation'], self.visu.id),
            (CRYPTO['sha512'], self.get_sha_sum(excursion_set_file)),
            ))

        return self.p

class Visualisation(NIDMObject):
    def __init__(self, visu_filename, stat_num, export_dir):
        super(Visualisation, self).__init__(export_dir)
        self.file = visu_filename
        self.id = NIIRI['excursion_set_png_id_'+str(stat_num)]


    def export(self):
        # Copy visualisation of excursion set in export directory
        shutil.copy(self.file, self.export_dir)
        path, visu_filename = os.path.split(self.file)

        # Create "png visualisation of Excursion set" entity
        self.p.entity(self.id, other_attributes=( 
            (PROV['type'], NIDM['Image']), 
            (PROV['location'], Identifier("file://./"+visu_filename)),
            (DCT['format'], "image/png"),
            ))

        return self.p

class HeightThreshold(NIDMObject):
    def __init__(self, stat_threshold=None, p_corr_threshold=None, p_uncorr_threshold=None):
        super(HeightThreshold, self).__init__()
        if not stat_threshold and not p_corr_threshold and not p_uncorr_threshold:
            raise Exception('No threshold defined')

        self.stat_threshold = stat_threshold
        self.p_corr_threshold = p_corr_threshold
        self.p_uncorr_threshold = p_uncorr_threshold

    def export(self):
        thresh_desc = ""
        if self.stat_threshold is not None:
            thresh_desc = "Z>"+str(self.stat_threshold)
            user_threshold_type = "Z-Statistic"
        elif self.p_uncorr_threshold is not None:
            thresh_desc = "p<"+str(self.p_uncorr_threshold)+" uncorr."
            user_threshold_type = "p-value uncorrected"
        elif self.p_corr_threshold is not None:
            thresh_desc = "p<"+str(self.p_corr_threshold)+" (GRF)"
            user_threshold_type = "p-value FWE"

        # FIXME: Do we want to calculate an uncorrected p equivalent to the Z thresh? 
        # FIXME: Do we want/Can we find a corrected p equivalent to the Z thresh? 
        heightThreshAllFields = {
            PROV['type']: NIDM['HeightThreshold'], PROV['label']: "Height Threshold: "+thresh_desc,
            NIDM['userSpecifiedThresholdType']: user_threshold_type , PROV['value']: self.stat_threshold,
            NIDM['pValueUncorrected']: self.p_uncorr_threshold, NIDM['pValueFWER']: self.p_corr_threshold
            }
        self.p.entity(NIIRI['height_threshold_id'], other_attributes=dict((k,v) for k,v in heightThreshAllFields.iteritems() if v is not None))

        return self.p

class ExtentThreshold(NIDMObject):
    def __init__(self, extent=None, p_corr=None, p_uncorr=None):
        super(ExtentThreshold, self).__init__()
        self.extent = extent
        self.p_corr = p_corr
        self.p_uncorr = p_uncorr

    def export(self):
        thresh_desc = ""
        if self.extent is not None:
            thresh_desc = "k>"+str(self.extent)
            user_threshold_type = "Cluster-size in voxels"
        elif not self.p_uncorr is None:
            thresh_desc = "p<"+str(self.p_uncorr)+" uncorr."
            user_threshold_type = "p-value uncorrected"
        elif not self.p_corr is None:
            thresh_desc = "p<"+str(self.p_corr)+" corr."
            user_threshold_type = "p-value FWE"
        extent_thresh_all_fields = {
            PROV['type']: NIDM['ExtentThreshold'], PROV['label']: "Extent Threshold: "+thresh_desc, NIDM['clusterSizeInVoxels']: self.extent,
            NIDM['userSpecifiedThresholdType']: user_threshold_type, NIDM['pValueUncorrected']: self.p_uncorr, NIDM['pValueFWER']: self.p_corr
        }
        self.p.entity(NIIRI['extent_threshold_id'], other_attributes=dict((k,v) for k,v in extent_thresh_all_fields.iteritems() if v is not None))

        return self.p

class Cluster(NIDMObject):
    def __init__(self, cluster_id, size, pFWER, peaks,
        x=None,y=None,z=None,x_std=None,y_std=None,z_std=None):
        super(Cluster, self).__init__()
        self.num = cluster_id
        self.id = NIIRI['cluster_000'+str(cluster_id)]
        self.center_of_gravity = CenterOfGravity(cluster_id, x,y,z,x_std,y_std,z_std)
        self.peaks = peaks
        self.size = size
        self.pFWER = pFWER

    def export(self):
        for peak in self.peaks:
            self.p.update(peak.export())
            self.p.wasDerivedFrom(peak.id, self.id)

        self.p.update(self.center_of_gravity.export())
        self.p.wasDerivedFrom(self.center_of_gravity.id, self.id)

        # FIXME deal with multiple contrasts
        self.p.entity(self.id, other_attributes=( 
                             (PROV['type'] , NIDM['Cluster']), 
                             (PROV['label'], "Cluster 000"+str(self.num)),
                             (NIDM['clusterSizeInVoxels'], self.size),
                             (NIDM['pValueFWER'], self.pFWER )))
        return self.p

class DisplayMaskMap(NIDMObject):
    def __init__(self, contrast_num, filename, coordinate_system, coordinate_space_id, export_dir):
        super(DisplayMaskMap, self).__init__(export_dir)
        self.id = NIIRI['display_map_id_'+contrast_num]
        self.filename = filename
        self.coord_space = CoordinateSpace(coordinate_system, coordinate_space_id, filename)

    def export(self):
        # Create coordinate space entity
        self.p.update(self.coord_space.export())

        # Create "Display Mask Map" entity
        display_mask_file = os.path.join(self.export_dir, 'DisplayMask.nii.gz')
        display_mask_original_filename, display_mask_filename = self.copy_nifti(self.filename, display_mask_file)     

        self.p.entity(self.id,
            other_attributes=(  (PROV['type'], NIDM['DisplayMaskMap']), 
                                (PROV['label'] , "Display Mask Map"),
                                (DCT['format'] , "image/nifti"),
                                (NIDM['inCoordinateSpace'], self.coord_space.id),
                                (NIDM['filename'], display_mask_original_filename),
                                (NIDM['filename'], display_mask_filename),
                                (PROV['location'], Identifier("file://./"+display_mask_filename)),
                                (CRYPTO['sha512'], self.get_sha_sum(display_mask_file))))
        return self.p
        
class PeakCriteria(NIDMObject):
    def __init__(self, contrast_num, num_peak, peak_dist):
        super(PeakCriteria, self).__init__()
        self.id = NIIRI['peak_definition_criteria_id_'+contrast_num]
        self.num_peak = num_peak
        self.peak_dist = peak_dist

    def export(self):
        # Create "Peak definition criteria" entity
        self.p.entity(self.id,
            other_attributes=(  (PROV['type'], NIDM['PeakDefinitionCriteria']), 
                                (PROV['label'] , "Peak Definition Criteria"),
                                (NIDM['maxNumberOfPeaksPerCluster'], self.num_peak),
                                (NIDM['minDistanceBetweenPeaks'], self.peak_dist)))

        return self.p


class ClusterCriteria(NIDMObject):
    def __init__(self, contrast_num, connectivity):
        super(ClusterCriteria, self).__init__()
        self.id = NIIRI['cluster_definition_criteria_id_'+contrast_num]
        self.connectivity = connectivity

    def export(self):
        # Create "Cluster definition criteria" entity
        voxel_connectivity = NIDM['connected'+str(self.connectivity)+'In3D']
        
        self.p.entity(self.id,
            other_attributes=(  (PROV['type'], NIDM['ClusterDefinitionCriteria']), 
                                (PROV['label'] , "Cluster Connectivity Criterion: "+str(self.connectivity)),
                                (NIDM['hasConnectivityCriterion'], voxel_connectivity)))

        return self.p
        

class CenterOfGravity(NIDMObject):
    def __init__(self, cluster_id, x=None,y=None,z=None,x_std=None,y_std=None,z_std=None):
        super(CenterOfGravity, self).__init__()
        self.cluster_id = cluster_id
        self.id = NIIRI['center_of_gravity_'+str(cluster_id)]
        self.coordinate = Coordinate(NIIRI['COG_coordinate_000'+str(cluster_id)], '000'+str(cluster_id),x,y,z,x_std,y_std,z_std)

    def export(self):
        self.p.update(self.coordinate.export())

        self.p.entity(self.id, other_attributes=( 
                     (PROV['type'] , FSL['CenterOfGravity']), 
                     (PROV['label'], "Center of gravity "+str(self.cluster_id)),
                     (PROV['location'] , NIIRI['COG_coordinate_000'+str(self.cluster_id)]))   )

        return self.p

class SearchSpace(NIDMObject):

    def __init__(self, search_space_file, search_volume, resel_size_in_voxels, dlh, 
        coordinate_system, coordinate_space_id, export_dir):
        super(SearchSpace, self).__init__(export_dir)
        self.file = search_space_file
        self.coord_space = CoordinateSpace(coordinate_system, coordinate_space_id, self.file) 
        self.resel_size_in_voxels = resel_size_in_voxels
        self.dlh = dlh
        self.search_volume = search_volume
        self.id = NIIRI['search_space_id']

    # Generate prov for search space entity generated by the inference activity
    def export(self):
        self.p.update(self.coord_space.export())

        # Copy "Mask map" in export directory
        search_space_original_file = self.file
        search_space_file = os.path.join(self.export_dir, 'SearchSpace.nii.gz')
        search_space_original_filename, search_space_filename = self.copy_nifti(search_space_original_file, search_space_file)
        
        # Crate "Mask map" entity
        self.p.entity(self.id, other_attributes=( 
                (PROV['label'], "Search Space Map"), 
                (DCT['format'], "image/nifti"), 
                (PROV['type'], NIDM['SearchSpaceMap']), 
                (PROV['location'], Identifier("file://./"+search_space_filename)),
                (NIDM['filename'], search_space_original_filename),
                (NIDM['filename'], search_space_filename),
                (NIDM['inCoordinateSpace'], self.coord_space.id),
                (FSL['searchVolumeInVoxels'], self.search_volume),
                (CRYPTO['sha512'], self.get_sha_sum(search_space_file)),
                (FSL['reselSizeInVoxels'], self.resel_size_in_voxels),
                (FSL['dlh'], self.dlh)))

        return self.p

class Coordinate(NIDMObject):

    def __init__(self, coordinate_id, label_id, x=None, y=None, z=None, x_std=None, y_std=None, z_std=None):
        super(Coordinate, self).__init__()   
        # FIXME: coordiinate_id should not be determined externally
        self.id = coordinate_id
        self.label_id = label_id
        self.x = x
        self.y = y
        self.z = z
        self.x_std = x_std
        self.y_std = y_std
        self.z_std = z_std

    def export(self):
        # We can not have this in the dictionnary because we want to keep the duplicate prov:type attribute
        typeAndLabelAttributes = [(PROV['type'],PROV['Location']),
            (PROV['type'], NIDM['Coordinate']),
            (PROV['label'], "Coordinate "+self.label_id)]

        coordinateAttributes = {
            FSL['coordinate1InVoxels'] : self.x,
            FSL['coordinate2InVoxels'] : self.y,
            FSL['coordinate3InVoxels'] : self.z,
            NIDM['coordinate1'] : self.x_std,
            NIDM['coordinate2'] : self.y_std,
            NIDM['coordinate3'] : self.z_std
            };

        self.p.entity(self.id, 
            other_attributes=typeAndLabelAttributes+\
                list(dict((k,v) for k,v in coordinateAttributes.iteritems() if not v is None).items()))

        return self.p

class Peak(NIDMObject):
    def __init__(self, cluster_index, peak_index, equiv_z, stat_num, max_peak, *args, **kwargs):
        super(Peak, self).__init__()
        # FIXME: Currently assumes less than 10 clusters per contrast
        cluster_id = cluster_index
        # FIXME: Currently assumes less than 100 peaks 
        peak_unique_id = '000'+str(cluster_index)+'_'+str(peak_index)
        self.id = NIIRI['peak_'+str(peak_unique_id)]
        self.num = peak_unique_id
        self.equiv_z = equiv_z
        self.max_peak = max_peak
        self.coordinate = Coordinate(NIIRI['coordinate_'+str(peak_unique_id)], str(peak_unique_id), **kwargs)

    def export(self):
        self.p.update(self.coordinate.export())

        other_attributes = [ 
            (PROV['type'] , NIDM['Peak']), 
            (PROV['label'] , "Peak "+str(self.num)), 
            (NIDM['equivalentZStatistic'], self.equiv_z), 
            (PROV['location'] , self.coordinate.id)]

        if self.max_peak:
            other_attributes.insert(0, (PROV['type'], FSL['ClusterMaximumStatistic']))

        self.p.entity(self.id, other_attributes=other_attributes)

        return self.p
