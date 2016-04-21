#!/usr/bin/env python
"""
Test of NIDM FSL export tool


@author: Camille Maumet <c.m.j.maumet@warwick.ac.uk>
@copyright: University of Warwick 2013-2014
"""
import os
import shutil
import sys
import glob
import json
import copy
import zipfile
import git
import subprocess

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


from TestCommons import *
from CheckConsistency import *


if __name__ == '__main__':
    config_file = os.path.join(TEST_DIR, 'config.json')
    if os.path.isfile(config_file):
        # Read config json file to find nidmresults-examples repository
        with open(config_file) as data_file:
            metadata = json.load(data_file)
        test_data_dir = metadata["data"]
    else:
        # Pull nidmresults-examples repository
        test_data_dir = os.path.join(TEST_DATA_DIR, "nidmresults-examples")

    if not os.path.isdir(os.path.join(test_data_dir, ".git")):
        logging.debug("Cloning to " + test_data_dir)
        # Cloning test data repository
        data_repo = git.Repo.clone_from(
            "https://github.com/incf-nidash/nidmresults-examples.git",
            test_data_dir)
    else:
        # Updating test data repository
        logging.debug("Updating repository at " + test_data_dir)
        subprocess.call(
            ["cd " + test_data_dir + "; git checkout new_ground_truth"],
            shell=True)
        subprocess.call(
            ["cd " + test_data_dir + "; git pull origin new_ground_truth"],
            shell=True)
        # "git stash" gives the repo one more chance to checkout the files
        # if something failed (e.g. git lfs error)
        subprocess.call(
            ["cd " + test_data_dir + "; git stash"],
            shell=True)
        # Just to check that everything went fine
        subprocess.call(
            ["cd " + test_data_dir + "; git status"],
            shell=True)

    # Find all test examples to be compared with ground truth
    test_data_cfg = glob.glob(os.path.join(test_data_dir, '*/config.json'))

    # # For test name readability remove path to test file
    # test_files = [x.replace(TEST_DIR, "") for x in test_files]
    # logging.info("Test files:\n\t" + "\n\t".join(test_files))

    #     version = metadata["version"]

    # *** Once for all, run the export
    EXPORTED_TEST_DIR = os.path.join(TEST_DIR, 'exported')
    if os.path.isdir(EXPORTED_TEST_DIR):
        shutil.rmtree(EXPORTED_TEST_DIR)
        os.mkdir(EXPORTED_TEST_DIR)

    for cfg in test_data_cfg:
        with open(cfg) as data_file:
            logging.debug(data_file)
            metadata = json.load(data_file)

        data_dir = os.path.dirname(cfg)

        if metadata["software"].lower() == "fsl":
            test_name = os.path.basename(data_dir)
            logging.debug("Computing NIDM FSL export for " + test_name)

            versions = metadata["versions"]

            for version in versions:
                version_str = version.replace(".", "")

                if os.path.isdir(data_dir):
                    # Remove existent NIDM exports (if an export already exist
                    # the program migth be stopped waiting on user output)
                    for nidmpack in glob.glob(os.path.join(
                            data_dir, "*.nidm.zip")):
                        os.remove(nidmpack)

                    # Export to NIDM using FSL export tool
                    # fslnidm = FSL_NIDM(feat_dir=DATA_DIR_001);
                    fslnidm = FSLtoNIDMExporter(
                        feat_dir=data_dir, version=version, zipped=True)
                    fslnidm.parse()
                    zipped_dir = fslnidm.export()
                    print 'NIDM export available at '+zipped_dir

                    # Copy provn export to test directory
                    test_export_dir = os.path.join(
                        EXPORTED_TEST_DIR,
                        'ex_' + test_name + '_' + version_str)

                    if not os.path.exists(test_export_dir):
                        os.makedirs(test_export_dir)

                    with zipfile.ZipFile(zipped_dir) as z:
                        z.extract('nidm.ttl', test_export_dir)
                        z.extract('nidm.provn', test_export_dir)

                    cfg_file = os.path.join(test_export_dir, 'config.json')

                    test_metadata = copy.copy(metadata)
                    del test_metadata['versions']
                    del test_metadata['software']
                    test_metadata['version'] = version

                    with open(cfg_file, 'w') as outfile:
                        json.dump(test_metadata,
                                  outfile,
                                  sort_keys=True,
                                  indent=4,
                                  separators=(',', ': '))

                    gt_dir = os.path.join(EXPORTED_TEST_DIR, '_ground_truth')
                    if not os.path.exists(gt_dir):
                        os.makedirs(gt_dir)

                    # with open(config_file) as config:
                    #     metadata = json.load(config)

                    for gt in metadata["ground_truth"]:
                        gt_file = os.path.join(
                            data_dir, "..", "_ground_truth", version, gt)
                        version_dir = os.path.join(gt_dir, version)
                        if not os.path.exists(version_dir):
                            os.makedirs(version_dir)

                        sub_gt_dir = os.path.join(
                            version_dir, os.path.dirname(gt))
                        if not os.path.exists(sub_gt_dir):
                            os.makedirs(sub_gt_dir)
                        shutil.copy(gt_file, os.path.join(
                            sub_gt_dir, os.path.basename(gt)))

                    # delete nidm export folder
                    os.remove(zipped_dir)
