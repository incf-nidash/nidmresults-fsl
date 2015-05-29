#!/usr/bin/env python
"""
Test of NIDM FSL export tool


@author: Camille Maumet <c.m.j.maumet@warwick.ac.uk>
@copyright: University of Warwick 2013-2014
"""
import unittest
import os
from rdflib.graph import Graph
import shutil
import sys
import glob
import json

import logging
logger = logging.getLogger(__name__)
# Display log messages in console
logging.basicConfig(level=logging.INFO, format='%(levelname)s - %(message)s')

RELPATH = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# Add FSL NIDM export to python path
sys.path.append(RELPATH)

# Add nidm common testing code folder to python path
NIDM_DIR = os.path.join(RELPATH, "nidm")
# In TravisCI the nidm repository will be created as a subtree, however locally
# the nidm directory will be accessed directly
logging.debug(NIDM_DIR)
if not os.path.isdir(NIDM_DIR):
    NIDM_DIR = os.path.join(os.path.dirname(RELPATH), "nidm")
    # The FSL export to NIDM will only be run locally (for now)
    from nidmfsl.fsl_exporter.fsl_exporter import FSLtoNIDMExporter

NIDM_RESULTS_DIR = os.path.join(NIDM_DIR, "nidm", "nidm-results")
TERM_RESULTS_DIRNAME = "terms"
TEST_DIR = os.path.dirname(os.path.abspath(__file__))
TEST_DATA_DIR = os.path.join(TEST_DIR, "data")

path = os.path.join(NIDM_RESULTS_DIR, "test")
sys.path.append(path)


from TestResultDataModel import TestResultDataModel, ExampleGraph
from TestCommons import *
from CheckConsistency import *

from ddt import ddt, data

# Find all test examples to be compared with ground truth
test_files = glob.glob(os.path.join(TEST_DIR, 'ex*', '*.ttl'))
logging.info("Test files:\n\t" + "\n\t".join(test_files))


@ddt
class TestFSLResultDataModel(unittest.TestCase, TestResultDataModel):

    @classmethod
    def setUpClass(cls):
        # *** Once for all, run the export
        for ttl in test_files:
            test_dir = os.path.dirname(ttl)

            # If test data is available (usually if the test is run locally)
            # then compute a fresh export
            with open(os.path.join(test_dir, 'config.json')) as data_file:
                metadata = json.load(data_file)
            data_dir = os.path.join(TEST_DATA_DIR, metadata["data_dir"])

            #  Turtle file obtained with FSL NI-DM export tool
            provn = ttl.replace(".ttl", ".provn")

            if os.path.isdir(data_dir):
                logging.debug("Computing NIDM FSL export")

                # Export to NIDM using FSL export tool
                # fslnidm = FSL_NIDM(feat_dir=DATA_DIR_001);
                fslnidm = FSLtoNIDMExporter(feat_dir=data_dir, version="1.0.0")
                fslnidm.parse()
                export_dir = fslnidm.export()
                # Copy provn export to test directory
                shutil.copy(os.path.join(export_dir, 'nidm.provn'),
                            os.path.join(provn))
                shutil.copy(os.path.join(export_dir, 'nidm.ttl'),
                            os.path.join(ttl))

    def setUp(self):
        # Retreive owl file for NIDM-Results
        owl_file = os.path.join(NIDM_RESULTS_DIR, TERM_RESULTS_DIRNAME,
                                'nidm-results.owl')
        import_files = glob.glob(
            os.path.join(os.path.dirname(owl_file),
                         os.pardir, os.pardir, "imports", '*.ttl'))

        TestResultDataModel.setUp(self, owl_file, import_files)

        self.ex_graphs = dict()

        for ttl in test_files:
            test_dir = os.path.dirname(ttl)
            with open(os.path.join(test_dir, 'config.json')) as data_file:
                metadata = json.load(data_file)
            gt_file = [os.path.join(NIDM_RESULTS_DIR, x)
                       for x in metadata["ground_truth"]]
            inclusive = metadata["inclusive"]

            self.ex_graphs[ttl] = ExampleGraph(
                owl_file, ttl, gt_file, inclusive)

    @data(*test_files)
    def test_class_consistency_with_owl(self, ttl):
        """
        Test: Check that the classes used in the ttl file are defined in the
        owl file.
        """

        ex = self.ex_graphs[ttl]

        # for ex in self.ex_graphs:
        # FIXME: change example name depending on graph
        my_exception = ex.owl.check_class_names(
            ex.graph, "FSL example00")

        # FIXME (error message display should be simplified as only one
        # example...)
        if my_exception:
            error_msg = ""
            for unrecognised_class_name, examples in my_exception.items():
                error_msg += unrecognised_class_name + \
                    " (from " + ', '.join(examples) + ")"
            raise Exception(error_msg)

    @data(*test_files)
    def test_attributes_consistency_with_owl(self, ttl):
        """
        Test: Check that the attributes used in the ttl file comply with their
        definition (range, domain) specified in the owl file.
        """

        ex = self.ex_graphs[ttl]

        my_exception = ex.owl.check_attributes(
            ex.graph, "FSL example001")

        # FIXME (error message display should be simplified as only one
        # example...)
        error_msg = ""
        if my_exception[0]:
            for unrecognised_attribute, example_names \
                    in my_exception[0].items():
                error_msg += unrecognised_attribute + \
                    " (from " + ', '.join(example_names) + ")"
        if my_exception[1]:
            for unrecognised_range, example_names \
                    in my_exception[1].items():
                error_msg += unrecognised_range + \
                    " (from " + ', '.join(example_names) + ")"

        if error_msg:
            raise Exception(error_msg)

    # FIXME: If terms PR is accepted then these tests should be moved to
    # TestResultDataModel.py
    @data(*test_files)
    def test_examples_match_ground_truth(self, ttl):
        """
        Test03: Comparing that the ttl file generated by FSL and the expected
        ttl file (generated manually) are identical
        """

        ex = self.ex_graphs[ttl]

        for gt_file in ex.gt_ttl_files:
            logging.info("Ground truth ttl: " + gt_file)

            # RDF obtained by the ground truth export
            gt = Graph()
            gt.parse(gt_file, format='turtle')

            self.compare_full_graphs(gt, ex.graph, ex.exact_comparison)

            if self.my_execption:
                raise Exception(self.my_execption)

if __name__ == '__main__':
    unittest.main()
