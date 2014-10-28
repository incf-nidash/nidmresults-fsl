from prov.model import ProvBundle
import numpy as np
import os
from constants import *
import nibabel as nib
import shutil
import hashlib

class NIDMObject(object):
    def __init__(self, export_dir=None, coordinate_space_id=None):
        self.export_dir = export_dir
        self.coordinate_space_id = coordinate_space_id
        self.p = ProvBundle()
        self.id = None

    def copy_nifti(self, original_file, new_file):
        shutil.copy(original_file, new_file)
        path, new_filename = os.path.split(new_file)
        path, original_filename = os.path.split(original_file)

        return original_filename, new_filename 

    def get_sha_sum(self, nifti_file):
        nifti_img = nib.load(nifti_file)
        data = nifti_img.get_data()
        # Fix needed as in https://github.com/pymc-devs/pymc/issues/327
        if not data.flags["C_CONTIGUOUS"]:
          data = np.ascontiguousarray(data)
        return hashlib.sha512(data).hexdigest()


# Generate prov for a coordinate space entity 
class CoordinateSpace(NIDMObject):
    def __init__(self, coordinate_system, coordinate_space_id, nifti_file):
        super(CoordinateSpace, self).__init__()
        self.coordinate_system = coordinate_system
        self.nifti_file = nifti_file
        self.id_num = coordinate_space_id
        self.id = NIIRI['coordinate_space_id_'+str(coordinate_space_id)]

    def export(self):
        thresImg = nib.load(self.nifti_file)
        thresImgHdr = thresImg.get_header()

        numDim = len(thresImg.shape)

        mydict = { 
            PROV['type']: NIDM['CoordinateSpace'], 
            NIDM['dimensionsInVoxels']: str(thresImg.shape).replace('(', '[').replace(')', ']'),
            NIDM['numberOfDimensions']: numDim,
            NIDM['voxelToWorldMapping']: '%s'%', '.join(str(thresImg.get_qform()).strip('()').replace('. ', '').split()).replace('[,', '[').replace('\n', ''),
            NIDM['inWorldCoordinateSystem']: self.coordinate_system,           
            # FIXME: this gives mm, sec => what is wrong: FSL file, nibabel, other?
            # NIDM['voxelUnits']: '[%s]'%str(thresImgHdr.get_xyzt_units()).strip('()'),
            NIDM['voxelUnits']: "['mm', 'mm', 'mm']",
            NIDM['voxelSize']: '[%s]'%', '.join(map(str, thresImgHdr['pixdim'][1:(numDim+1)])),
            PROV['label']: "Coordinate space "+str(self.id_num)}
        self.p.entity(self.id, other_attributes=mydict)
        return self.p


