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
import subprocess
import argparse

import logging
logger = logging.getLogger(__name__)
# Display log messages in console
logging.basicConfig(filename='debug.log', level=logging.DEBUG, filemode='w',
                    format='%(levelname)s - %(message)s')

RELPATH = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# Add FSL NIDM export to python path
sys.path.append(RELPATH)

TEST_DIR = os.path.dirname(os.path.abspath(__file__))
TEST_DATA_DIR = os.path.join(TEST_DIR, "data")

if __name__ == '__main__':
    # Arguments and description
    parser = argparse.ArgumentParser(
        description='Export FSL NIDM-Results test data as NIDM packs.')
    parser.add_argument(
        "-o", "--out_dir",
        help='Path to output directory.', default=None)
    args = parser.parse_args()

    test_only = True
    if args.out_dir is None:
        EXPORTED_TEST_DIR = os.path.join(TEST_DIR, 'exported')
    else:
        test_only = False
        EXPORTED_TEST_DIR = args.out_dir

    # # Parse feat dir and export to NIDM
    # fslnidm = FSLtoNIDMExporter(
    #     out_dirname=args.output_name, zipped=(not args.directory_output),
    #     version=args.nidm_version, feat_dir=args.feat_dir, groups=args.group)
    # fslnidm.parse()
    # output_path = fslnidm.export()

    ###############

    config_file = os.path.join(TEST_DIR, 'config.json')
    if os.path.isfile(config_file):
        # Read config json file to find nidmresults-examples repository
        with open(config_file) as data_file:
            metadata = json.load(data_file)
        test_data_dir = metadata["data"]
    else:
        # Pull nidmresults-examples repository
        test_data_dir = os.path.join(TEST_DATA_DIR, "nidmresults-examples")

    # Find all test examples to be compared with ground truth
    test_data_cfg = glob.glob(os.path.join(test_data_dir, '*/config.json'))

    # # For test name readability remove path to test file
    # test_files = [x.replace(TEST_DIR, "") for x in test_files]
    # logging.info("Test files:\n\t" + "\n\t".join(test_files))

    #     version = metadata["version"]

    # *** Once for all, run the export
    if os.path.isdir(EXPORTED_TEST_DIR):
        if not test_only:
            msg = EXPORTED_TEST_DIR+" already exists, overwrite?"
            if not input("%s (y/N) " % msg).lower() == 'y':
                quit("Bye.")
            if os.path.isdir(EXPORTED_TEST_DIR):
                shutil.rmtree(EXPORTED_TEST_DIR)
            else:
                os.remove(EXPORTED_TEST_DIR)
        else:
            shutil.rmtree(EXPORTED_TEST_DIR)

    os.mkdir(EXPORTED_TEST_DIR)

    for cfg in test_data_cfg:
        with open(cfg) as data_file:
            logging.debug(data_file)
            metadata = json.load(data_file)

        data_dir = os.path.dirname(cfg)
        if data_dir.endswith("/"):
            data_dir = data_dir[:-1]

        test_name = os.path.basename(data_dir)
        if metadata["software"].lower() == "fsl":
            logging.debug("Computing NIDM FSL export for " + test_name)

            versions = metadata["versions"]

            if "group_names" in metadata:
                group_names = metadata["group_names"]
                num_subjects = metadata["num_subjects"]
            else:
                group_names = None
                num_subjects = None

            for version in versions:
                version_str = version.replace(".", "")

                if os.path.isdir(data_dir):
                    # Remove existent NIDM exports (if an export already exist
                    # the program migth be stopped waiting on user output)
                    for nidmpack in glob.glob(os.path.join(
                            data_dir, "*.nidm.zip")):
                        os.remove(nidmpack)

                    if test_name == "fsl_full_examples001":

                        # For our test case full_examples001 we need to change
                        # path to the onset files (stored in the feat folder)
                        # so that the type of model (mixed, event, block)
                        # can be retreived
                        fsf_file = os.path.join(data_dir, "design.fsf")
                        fsf_cp = os.path.join(data_dir, "design_cp.fsf")
                        shutil.copy(fsf_file, fsf_cp)

                        with open(fsf_file, 'r') as fsf:

                            design = fsf.read()

                        with open(fsf_file, 'w') as fsf:

                            fsf.write(design.replace(
                                "/storage/wmsmfe/fsl_course_data/fmri_fluency",
                                data_dir))

                    # Export to NIDM using FSL export tool
                    # fslnidm = FSL_NIDM(feat_dir=DATA_DIR_001);
                    featdir_arg = str(data_dir)
                    group_arg = ""
                    if num_subjects and \
                            version not in ["1.0.0", "1.1.0", "1.2.0"]:
                        for label, numsub in \
                                list(zip(group_names, num_subjects)):
                            group_arg += " -g " + label + " " + str(numsub)
                    if version:
                        version_arg = " -n " + version

                    nidmfsl_cmd = [
                        "nidmfsl " + featdir_arg + group_arg + version_arg]
                    print("\nRunning " + str(nidmfsl_cmd))
                    subprocess.check_call(nidmfsl_cmd, shell=True)

                    if test_name == "fsl_full_examples001":
                        shutil.move(fsf_cp, fsf_file)

                    zipped_dir = os.path.join(
                        data_dir, os.path.basename(data_dir) + ".nidm.zip")

                    # Copy ttl export to test directory
                    test_export_dir = os.path.join(
                        EXPORTED_TEST_DIR,
                        'ex_' + test_name + '_' + version_str)

                    if test_only:
                        if not os.path.exists(test_export_dir):
                            os.makedirs(test_export_dir)
                        with zipfile.ZipFile(zipped_dir) as z:
                            z.extract('nidm.ttl', test_export_dir)
                            z.extract('nidm.json', test_export_dir)
                    else:
                        zipped_ver_dir = os.path.join(
                            EXPORTED_TEST_DIR,
                            os.path.basename(data_dir) +
                            "_" + version.replace(".", "") + ".nidm.zip")
                        shutil.copyfile(zipped_dir, zipped_ver_dir)

                    cfg_file = os.path.join(test_export_dir, 'config.json')

                    test_metadata = copy.copy(metadata)
                    del test_metadata['versions']
                    del test_metadata['software']
                    test_metadata['version'] = version

                    if test_only:
                        with open(cfg_file, 'w') as outfile:
                            json.dump(test_metadata,
                                      outfile,
                                      sort_keys=True,
                                      indent=4,
                                      separators=(',', ': '))

                        gt_dir = os.path.join(
                            EXPORTED_TEST_DIR, '_ground_truth')
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
