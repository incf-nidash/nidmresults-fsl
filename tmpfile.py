import nibabel as nib
import numpy as np
import numpy.linalg as npla

tmpfile = '/home/tommaullin/Documents/nidmresults-fsl/test/data/nidmresults-examples/fsl_thr_clustfwep05/cluster_zstat1.txt'
exc_set = '/home/tommaullin/Documents/nidmresults-fsl/test/data/nidmresults-examples/fsl_thr_clustfwep05/thresh_zstat1.nii.gz'

numberarray = np.loadtxt(tmpfile, skiprows=1)
tab_hdr = 'Cluster Index	Voxels	P	-log10(P)	Z-MAX	Z-MAX X (vox)	Z-MAX Y (vox)	Z-MAX Z (vox)	Z-COG X (vox)	Z-COG Y (vox)	Z-COG Z (vox)	COPE-MAX	COPE-MAX X (vox)	COPE-MAX Y (vox)	COPE-MAX Z (vox)	COPE-MEAN'

print(numberarray.shape)

# Replace first 3 columns with coordinates
zmax_coords = np.insert(numberarray[:,5:8], 3, 1, axis=1)
zcog_coords = np.insert(numberarray[:,8:11], 3, 1, axis=1)
copemax_coords = np.insert(numberarray[:,12:15], 3, 1, axis=1)

exc_set_img = nib.load(exc_set)

# Transformation matrix from voxels to mm
voxToWorld = exc_set_img.affine

# Intercept.
intrcp = voxToWorld[:3, 3]

print(intrcp)

zmax_mm = np.dot(zmax_coords, voxToWorld)

zcog_mm = np.dot(zcog_coords, voxToWorld)

copemax_mm = np.dot(copemax_coords, voxToWorld)

print(zmax_mm[:, :3])
print(zmax_mm[:, :3] + intrcp)

numberarray[:,5:8] = zmax_mm[:, :3]
numberarray[:,8:11] = zcog_mm[:, :3]
numberarray[:,12:15] = copemax_mm[:, :3]

print(numberarray.shape)

cluster_mm_file ='/home/tommaullin/Documents/nidmresults-fsl/test/data/nidmresults-examples/fsl_thr_clustfwep05/cluster_zstat1_sub.txt'
 # cluster_mm_file = os.path.join(analysis_dir, 'cluster_' + prefix + str(stat_num) + '_std.txt')

np.savetxt(cluster_mm_file, numberarray, header=tab_hdr, comments='', fmt='%i %i %.2e %3g %3g %i %i %i %3g %3g %3g %i %i %i %i %i')
