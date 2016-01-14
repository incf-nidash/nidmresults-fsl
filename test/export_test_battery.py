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
logging.basicConfig(filename='debug.log', level=logging.DEBUG, filemode='w',
                    format='%(levelname)s - %(message)s')

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


from TestResultDataModel import TestResultDataModel
from TestCommons import *
from CheckConsistency import *


if __name__ == '__main__':
    # Read config json file to find nidmresults-examples repository
    with open(os.path.join(TEST_DIR, 'config.json')) as data_file:
        metadata = json.load(data_file)
    data_dir = os.path.join(TEST_DATA_DIR, metadata["data"])
    print data_dir

    # Find all test examples to be compared with ground truth
    test_data = glob.glob(os.path.join(data_dir, 'fsl_*'))
    print test_data

    # # For test name readability remove path to test file
    # test_files = [x.replace(TEST_DIR, "") for x in test_files]
    # logging.info("Test files:\n\t" + "\n\t".join(test_files))

    #     version = metadata["version"]

    # *** Once for all, run the export
    for data_dir in test_data:
        test_name = os.path.basename(data_dir)

        version = "1.2.0"

        if os.path.isdir(data_dir):
            logging.debug("Computing NIDM FSL export")

            # Export to NIDM using FSL export tool
            # fslnidm = FSL_NIDM(feat_dir=DATA_DIR_001);
            fslnidm = FSLtoNIDMExporter(feat_dir=data_dir, version=version)
            fslnidm.parse()
            export_dir = fslnidm.export()

            # Copy provn export to test directory
            test_export_dir = os.path.join(TEST_DIR, 'ex_' + test_name)
            if not os.path.exists(test_export_dir):
                os.makedirs(test_export_dir)
            shutil.copy(os.path.join(export_dir, 'nidm.provn'),
                        os.path.join(test_export_dir, 'nidm.provn'))
            shutil.copy(os.path.join(export_dir, 'nidm.ttl'),
                        os.path.join(test_export_dir, 'nidm.ttl'))

            config_file = os.path.join(data_dir, 'config.json')
            shutil.copy(os.path.join(data_dir, 'config.json'),
                        os.path.join(test_export_dir, 'config.json'))

            gt_dir = os.path.join(TEST_DIR, 'ground_truth')
            if not os.path.exists(gt_dir):
                os.makedirs(gt_dir)

            with open(config_file) as config:
                metadata = json.load(config)

            for gt in metadata["ground_truth"]:
                gt_file = os.path.join(data_dir, "..", "ground_truth", gt)
                sub_gt_dir = os.path.join(gt_dir, os.path.dirname(gt))
                if not os.path.exists(sub_gt_dir):
                    os.makedirs(sub_gt_dir)
                shutil.copy(gt_file, os.path.join(
                    sub_gt_dir, os.path.basename(gt)))