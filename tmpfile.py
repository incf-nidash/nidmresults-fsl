import nibabel as nib
import numpy as np
import numpy.linalg as npla

tmpfile = '/home/tommaullin/Documents/nidmresults-fsl/test/data/nidmresults-examples/fsl_thr_clustfwep05/cluster_zstat1.txt'
exc_set = '/home/tommaullin/Documents/nidmresults-fsl/test/data/nidmresults-examples/fsl_thr_clustfwep05/thresh_zstat1.nii.gz'

numberarray = np.loadtxt(tmpfile, skiprows=1)

print(numberarray)

# Read in dimensions
dim = numberarray.shape

# Create arrays of ones
zmax_coords = np.ones((dim[0], 4))
zcog_coords = np.ones((dim[0], 4))
copemax_coords = np.ones((dim[0], 4))

# Replace first 3 columns with coordinates
zmax_coords[:, :3] = numberarray[:,5:8]
zcog_coords[:, :3] = numberarray[:,8:11]
copemax_coords[:, :3] = numberarray[:,12:15]

exc_set_img = nib.load(exc_set)

# Transformation matrix from voxels to mm
voxToWorld = exc_set_img.affine

# Transformation matrix from mm back to voxels
worldToVox = npla.inv(voxToWorld)



print(zmax_coords)

print(zcog_coords)
print(copemax_coords)

zmax_mm = np.dot(zmax_coords, voxToWorld)

zcog_mm = np.dot(zcog_coords, voxToWorld)

copemax_mm = np.dot(copemax_coords, voxToWorld)



print(zmax_mm)
print(zcog_mm)
print(copemax_mm)
