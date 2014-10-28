from prov.model import Identifier
import numpy as np
import os
from constants import *
import nibabel as nib
import hashlib
from generic import *

class Contrast(NIDMObject):
    def __init__(self, contrast_num, contrast_name, weights, estimation, contrast_map, stderr_map, 
        stat_map, z_stat_map=None):
        super(Contrast, self).__init__()
        self.contrast_num = contrast_num
        self.contrast_name = contrast_name
        self.weights = weights
        self.estimation = estimation
        self.contrast_map = contrast_map
        self.stderr_map = stderr_map
        self.stat_map = stat_map
        self.z_stat_map = z_stat_map

    def export(self):
        # Estimation activity
        self.p.update(self.estimation.export())

        # Contrast weights
        self.p.update(self.weights.export())

        # Contrast Map
        self.p.update(self.contrast_map.export())
        self.p.wasGeneratedBy(self.contrast_map.id, self.estimation.id)

        # Standard Error Map
        self.p.update(self.stderr_map.export())
        self.p.wasGeneratedBy(self.stderr_map.id, self.estimation.id)

        # Statistic Map
        self.p.update(self.stat_map.export())
        self.p.wasGeneratedBy(self.stat_map.id, self.estimation.id)

        # Z Statistic Map
        if self.z_stat_map:
            self.p.update(self.z_stat_map.export())
            self.p.wasGeneratedBy(self.z_stat_map.id, self.estimation.id)

        return self.p

class ContrastWeights(NIDMObject):
    def __init__(self, contrast_num, contrast_name, contrast_weights):
        super(ContrastWeights, self).__init__()
        self.contrast_name = contrast_name
        self.contrast_weights = contrast_weights
        self.contrast_num = contrast_num
        self.id = NIIRI['contrast_id_'+self.contrast_num]

    def export(self):
        # Contrast id entity
        # FIXME: Get contrast weights
        # FIXME: Deal with F weights
        self.p.entity(self.id, 
            other_attributes=( (PROV['type'], NIDM['ContrastWeights']), 
                               (NIDM['statisticType'], NIDM['TStatistic']),
                               (PROV['label'], "Contrast Weights: "+self.contrast_name), 
                               (NIDM['contrastName'], self.contrast_name),
                               (PROV['value'], self.contrast_weights)))
        return self.p

class ContrastMap(NIDMObject):
    def __init__(self, contrast_file, contrast_num, contrast_name, coordinate_system, coordinate_space_id, export_dir):
        super(ContrastMap, self).__init__(export_dir)
        self.file = contrast_file
        self.num = contrast_num
        self.name = contrast_name
        self.id = NIIRI['contrast_map_id_'+contrast_num]
        self.coord_space = CoordinateSpace(coordinate_system, coordinate_space_id, self.file)

    def export(self):
        self.p.update(self.coord_space.export())

        # Copy contrast map in export directory
        cope_file = os.path.join(self.export_dir, 'Contrast.nii.gz')
        cope_original_filename, cope_filename = self.copy_nifti(self.file, cope_file)

        # Contrast Map entity
        path, filename = os.path.split(cope_file)
        self.p.entity(self.id, other_attributes=( 
            (PROV['type'], NIDM['ContrastMap']), 
            (DCT['format'], "image/nifti"), 
            (NIDM['inCoordinateSpace'], self.coord_space.id),
            (PROV['location'], Identifier("file://./"+cope_filename)),
            (NIDM['filename'], cope_original_filename),
            (NIDM['filename'], cope_filename),
            (NIDM['contrastName'], self.name),
            (CRYPTO['sha512'], self.get_sha_sum(cope_file)),
            (PROV['label'], "Contrast Map: "+self.name)))        
        return self.p


class ContrastStdErrMap(NIDMObject):
    def __init__(self, contrast_num, filename, is_variance, coordinate_system, coordinate_space_id, export_dir):
        super(ContrastStdErrMap, self).__init__(export_dir)
        self.file = filename
        self.id = NIIRI['contrast_standard_error_map_id_'+contrast_num]
        self.is_variance = is_variance
        self.coord_space = CoordinateSpace(coordinate_system, coordinate_space_id, filename)
        if is_variance:
            self.var_coord_space = CoordinateSpace(coordinate_system, coordinate_space_id+1, filename)

    def export(self):
        self.p.update(self.coord_space.export())

        standard_error_file = os.path.join(self.export_dir, "ContrastStandardError.nii.gz")
        if self.is_variance:
            self.p.update(self.var_coord_space.export())

            # Copy contrast variance map in export directory
            path, var_cope_filename = os.path.split(self.file)
            # FIXME: Use ContrastVariance.nii.gz?
            # var_cope_file = os.path.join(self.export_dir, var_cope_filename)
            # var_cope_original_filename, var_cope_filename = self.copy_nifti(var_cope_original_file, var_cope_file)

            # Contrast Variance Map entity
            # self.provBundle.entity('niiri:'+'contrast_variance_map_id_'+contrast_num, other_attributes=( 
            contrast_var_id = NIIRI[hashlib.md5(self.get_sha_sum(self.file)).hexdigest()]
            
            self.p.entity(contrast_var_id, other_attributes=( 
                (PROV['type'], FSL['ContrastVarianceMap']), 
                (NIDM['inCoordinateSpace'], self.var_coord_space.id),
                (CRYPTO['sha512'], self.get_sha_sum(self.file)),
                (NIDM['filename'], var_cope_filename)))
            
            # Create standard error map from contrast variance map
            var_cope_img = nib.load(self.file)
            contrast_variance = var_cope_img.get_data()

            standard_error_img = nib.Nifti1Image(np.sqrt(contrast_variance), var_cope_img.get_qform())
            nib.save(standard_error_img, standard_error_file)

        else:
            standard_error_original_file, standard_error_file = self.copy_nifti(self.file, standard_error_file)

        path, filename = os.path.split(standard_error_file)
        self.p.entity(self.id, other_attributes=( 
            (PROV['type'], NIDM['ContrastStandardErrorMap']), 
            (DCT['format'], "image/nifti"), 
            (NIDM['inCoordinateSpace'], self.coord_space.id),
            (PROV['location'], Identifier("file://./"+filename)),
            (CRYPTO['sha512'], self.get_sha_sum(standard_error_file)),
            (NIDM['filename'], filename),
            (PROV['label'], "Contrast Standard Error Map")))
        
        if self.is_variance:
            self.p.wasDerivedFrom(self.id, contrast_var_id)

        return self.p

class StatisticMap(NIDMObject):
    def __init__(self, stat_file, stat_type, contrast_num, contrast_name, dof,
                coordinate_system, coordinate_space_id, export_dir):
        super(StatisticMap, self).__init__(export_dir)
        self.num = contrast_num
        self.name = contrast_name
        self.file = stat_file
        if stat_type == 'Z':
            self.id = NIIRI['z_statistic_map_id_'+contrast_num ]
        else:
            self.id = NIIRI['statistic_map_id_'+contrast_num ]
        self.coord_space = CoordinateSpace(coordinate_system, coordinate_space_id, stat_file)
        self.stat_type = stat_type
        self.dof = dof

    def export(self):
        self.p.update(self.coord_space.export())

        # Copy Statistical map in export directory
        stat_file = os.path.join(self.export_dir, self.stat_type+'Statistic.nii.gz')
        stat_original_filename, stat_filename = self.copy_nifti(self.file, stat_file)       

        label = "Statistic Map: "+self.name
        if self.stat_type == 'Z':
            label = self.stat_type+'-'+label


        attributes = [(PROV['type'], NIDM['StatisticMap']), 
                                        (DCT['format'], "image/nifti"), 
                                        (PROV['label'], label) ,
                                        (PROV['location'], Identifier("file://./"+stat_filename)),
                                        (NIDM['statisticType'], NIDM[self.stat_type+'Statistic']), 
                                        (NIDM['filename'], stat_filename),
                                        (NIDM['filename'], stat_original_filename),
                                        (NIDM['contrastName'], self.name),
                                        (CRYPTO['sha512'], self.get_sha_sum(stat_file)),
                                        (NIDM['inCoordinateSpace'], self.coord_space.id)]

        if not self.stat_type == 'Z':
            attributes.insert(0, (NIDM['errorDegreesOfFreedom'], self.dof))
            # Check if we should have effectDegreesOfFreedom with z-stat?
            attributes.insert(0, (NIDM['effectDegreesOfFreedom'], 1.0))

        # Create "Statistic Map" entity
        # FIXME: Deal with other than t-contrast maps: dof + statisticType
        self.p.entity(self.id,
            other_attributes=attributes )
        return self.p

# class ZStatisticMap(StatisticMap)
#     # Copy Z-Statistical map in export directory
#     z_stat_file = os.path.join(self.export_dir, 'ZStatistic.nii.gz')
#     z_stat_original_filename, z_stat_filename = self.copy_nifti(z_stat_original_file, z_stat_file)

#     # Create "Z-Statistic Map" entity
#     self.provBundle.entity(NIIRI['z_statistic_map_id_'+contrast_num ],
#         other_attributes=(  (PROV['type'], NIDM['StatisticMap']), 
#                             (DCT['format'], "image/nifti"), 
#                             (PROV['label'], "Z-Statistic Map: "+contrast_name) ,
#                             (PROV['location'], Identifier("file://./"+z_stat_filename)),
#                             (NIDM['statisticType'], NIDM['ZStatistic']), 
#                             (NIDM['contrastName'], contrast_name),
#                             (NIDM['filename'], z_stat_original_filename),
#                             (NIDM['filename'], z_stat_filename),
#                             (CRYPTO['sha512'], self.get_sha_sum(z_stat_file)),
#                             (NIDM['inCoordinateSpace'], self.create_coordinate_space(z_stat_file)),
#                             ) )

class ContrastEstimation(NIDMObject):
    def __init__(self, contrast_num, contrast_name):
        super(ContrastEstimation, self).__init__()
        self.num = contrast_num
        self.name = contrast_name
        self.id = NIIRI['contrast_estimation_id_'+self.num]

    def export(self):
        self.p.activity(self.id, other_attributes=( 
            (PROV['type'], NIDM['ContrastEstimation']),
            (PROV['label'], "Contrast estimation: "+self.name)))

        return self.p

