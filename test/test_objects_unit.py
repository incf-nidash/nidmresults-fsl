import unittest
import os
import sys
import numpy as np

from prov.model import ProvDocument

import logging
logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.DEBUG)

RELPATH = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# Add FSL NIDM export to python path
sys.path.append(RELPATH)

# Add nidm common testing code folder to python path
NIDM_DIR = os.path.join(RELPATH, "nidm")
# In TravisCI the nidm repository will be created as a subtree, however locally the nidm
# directory will be accessed directly
logging.debug(NIDM_DIR)
if not os.path.isdir(NIDM_DIR):
    NIDM_DIR = os.path.join(os.path.dirname(RELPATH), "nidm")


NIDM_RESULTS_DIR = os.path.join(NIDM_DIR, "nidm", "nidm-results")
TERM_RESULTS_DIR = os.path.join(NIDM_RESULTS_DIR, "terms")
TEST_FOLDER = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'example001')
DATA_DIR = os.path.join(RELPATH, 'test', 'data', 'fmri.feat')

path = os.path.join(NIDM_RESULTS_DIR, "test")
sys.path.append(path)

NIDM_SCRIPT_DIR = os.path.join(NIDM_RESULTS_DIR, "scripts")
sys.path.append(NIDM_SCRIPT_DIR)

from TestCommons import *
from OwlReader import OwlReader
from CheckConsistency import *

from exporter.objects.modelfitting import *

from subprocess import call
from rdflib.graph import Graph
from Constants import *



class NIDMObjectsUnitTesting(unittest.TestCase):
    """
    Unit testing of NIDM objects (compared to examples provided in 
    nidm-results.owl)
    """

    def setUp(self):
        self.export_dir = os.path.join(TEST_FOLDER, 'nidm')
        if not os.path.isdir(self.export_dir):
            os.mkdir(self.export_dir)

        # Retreive owl file for NIDM-Results
        owl_file = os.path.join(TERM_RESULTS_DIR, 'nidm-results.owl')
        assert owl_file
        self.owl = OwlReader(owl_file)

        self.doc = ProvDocument()
        # self.bundle = ProvBundle(identifier=NIIRI[software_lc+'_results_id'])

        self.provn_file = os.path.join(self.export_dir, 'unit_test.provn')

        namespaces_file = os.path.join(TERM_RESULTS_DIR, "templates", \
            "Namespaces.txt")
        namespaces_fid = open(namespaces_file)
        self.prefixes = namespaces_fid.read()
        namespaces_fid.close()

        self.to_delete_files = [self.provn_file]
        self.gt_ttl_files = list()

    def test_design_matrix(self):
        mat = np.matrix('1 2; 3 4')

        mat_image = os.path.join(os.path.dirname(TEST_FOLDER), "data", \
            "fmri.feat", "design.png")

        design_matrix = DesignMatrix(mat, mat_image, self.export_dir)
        self.doc.update(design_matrix.export())

        # In the FSL export the design matrix contains both the Design Matrix
        # entity and the Image entity representing the design matrix 
        # visualisation.
        self.to_delete_files.append(os.path.join(self.export_dir, \
            "DesignMatrix.csv"))
        self.to_delete_files.append(os.path.join(self.export_dir, \
            "DesignMatrix.png")) 

        gt_file = self.owl.get_example(NIDM['DesignMatrix'])
        self.gt_ttl_files = [os.path.join(TERM_RESULTS_DIR, \
            gt_file.replace("file://./", "")), 
            os.path.join(TERM_RESULTS_DIR, "examples", "Image-DesignMatrix.txt")]

        self._create_gt_and_compare("Design Matrix")

    def test_data(self):
        data = Data(grand_mean_scaling=True, target=100.0)
        self.doc.update(data.export())

        gt_file = self.owl.get_example(NIDM['Data'])
        self.gt_ttl_files.append(os.path.join(TERM_RESULTS_DIR, \
            gt_file.replace("file://./", "")))

        self._create_gt_and_compare("Data")

# INDEPEDENT_CORR = NIDM['IndependentError']
# SERIALLY_CORR = NIDM['SeriallyCorrelatedError']
# COMPOUND_SYMMETRY_CORR = NIDM['CompoundSymmetricError']
# ARBITRARILY_CORR = NIDM['ArbitriralyCorrelatedError']


    # def test_error_model_indepdt_global(self):
    #     error_distribution = GAUSSIAN_DISTRIBUTION
    #     variance_homo = True
    #     variance_spatial = SPATIALLY_GLOBAL
    #     dependance = INDEPEDENT_CORR
    #     dependance_spatial = SPATIALLY_GLOBAL

    #     error_model = ErrorModel(error_distribution, variance_homo, 
    #         variance_spatial, dependance, dependance_spatial)
    #     self.doc.update(error_model.export())

    #     nidm_classes = {
    #         "ErrorModel": dict(
    #             error_model_id="niiri:error_model_id",
    #             noise_distribution="nidm:GaussianDistribution",
    #             variance_homo="true",
    #             variance_spatial="nidm:SpatiallyGlobal",
    #             dependence="nidm:IndependentError",
    #             dependence_spatial="nidm:SpatiallyLocal"
    #         )
    #         }
    #     self._create_gt_and_compare(nidm_classes, "Data")

    def _create_gt_and_compare(self, class_name):
        # Write-out current example in a provn file and convert to turtle
        provn_fid = open(self.provn_file, 'w')
        provn_fid.write(self.doc.get_provn())
        provn_fid.close()

        ttl_file = self.provn_file.replace(".provn", ".ttl")

        call("provconvert -infile "+self.provn_file+" -outfile "+ttl_file, \
            shell=True)

        self.to_delete_files.append(ttl_file)

        # Load current example graph
        ex_graph = Graph()
        ex_graph.parse(source=ttl_file, format='turtle')
       
        # Read and concatenate ground truth files
        gt = ""
        for gt_ttl_file in self.gt_ttl_files:
            gt_fid = open(gt_ttl_file)
            # What is described in the examples to be at any path is relative 
            # in export
            gt = gt.replace("/path/to/", "./")
            gt = gt+gt_fid.read()
            gt_fid.close()

        gt_graph = Graph()
        gt = self.prefixes+gt
        gt_graph.parse(data=gt, format='turtle')

        # Compare graphs
        found_diff = compare_graphs(ex_graph, gt_graph)

        if found_diff:
            raise Exception("Difference in "+class_name+".")

    def tearDown(self):
        # Delete files created for testing
        for to_delete_file in self.to_delete_files:
            if os.path.isfile(to_delete_file):
                os.remove(to_delete_file)

        os.rmdir(self.export_dir)

if __name__ == '__main__':
    unittest.main()