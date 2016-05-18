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
    def retry_lfs_download(cmd):
        MAX_ITER = 200
        it = 0
        while (it < MAX_ITER):
            # "git stash" gives the repo one more chance to checkout
            # the git-lfs files
            try:
                print cmd
                out = subprocess.check_output(cmd, shell=True)
                print out
                break
            except subprocess.CalledProcessError as e:
                if e.returncode == 128:
                    it = it + 1
                    print 'Retry #'+str(it)

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
        try:
            logging.debug("Cloning to " + test_data_dir)
            repo_https = \
                "https://github.com/incf-nidash/nidmresults-examples.git"
            clone_cmd = ["cd " + test_data_dir + "; git clone " + repo_https]
            print clone_cmd
            subprocess.check_call(clone_cmd, shell=True)
        except subprocess.CalledProcessError as e:
            # 128 -> git-lfs download error: "Error downloading object"
            if e.returncode == 128:
                # "git stash" gives the repo one more chance to checkout the
                # git-lfs files if the download failed
                stash_cmd = ["cd " + test_data_dir + "; git stash"]
                retry_lfs_download(stash_cmd)
    else:
        # Updating test data repository
        logging.debug("Updating repository at " + test_data_dir)

        # Check current branch and status
        subprocess.call(["cd " + test_data_dir + "; git branch"], shell=True)
        subprocess.call(["cd " + test_data_dir + "; git status"], shell=True)

        # If we are in a different branch, checkout (this test is useful so
        # that we only stash untracked files if in a different bramch)
        branch_name = "new_ground_truth"
        try:
            # Start from a clean state: stash local changes including
            # untracked and then checkout
            checkout_cmd = ["cd " + test_data_dir +
                            "; git stash --include-untracked" +
                            "; git checkout " + branch_name]
            print checkout_cmd
            subprocess.check_call(checkout_cmd, shell=True)
        except subprocess.CalledProcessError as e:
            # 128 -> git-lfs download error: "Error downloading object"
            if e.returncode == 128:
                retry_lfs_download(checkout_cmd)

        # Pull latest updates
        pull_cmd = ["cd " + test_data_dir +
                    "; git pull origin " + branch_name]
        try:
            print pull_cmd
            subprocess.check_call(pull_cmd, shell=True)
        except subprocess.CalledProcessError as e:
            # 128 -> git-lfs download error: "Error downloading object"
            if e.returncode == 128:
                # "git stash" gives the repo one more chance to checkout the
                # git-lfs files if the download failed
                stash_cmd = ["cd " + test_data_dir + "; git stash"]
                retry_lfs_download(stash_cmd)

        # Check current branch and status
        subprocess.call(["cd " + test_data_dir + "; git branch"], shell=True)
        subprocess.call(["cd " + test_data_dir + "; git status"], shell=True)

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
        if data_dir.endswith("/"):
            data_dir = data_dir[:-1]

        if metadata["software"].lower() == "fsl":
            test_name = os.path.basename(data_dir)
            logging.debug("Computing NIDM FSL export for " + test_name)

            versions = metadata["versions"]

            if "num_subjects" in metadata:
                num_subjects = metadata["num_subjects"]
            else:
                num_subjects = None
            if "group_names" in metadata:
                group_names = metadata["group_names"]
            else:
                group_names = None

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
                    featdir_arg = str(data_dir)
                    numsubs_arg = ""
                    groupnmes_arg = ""
                    if num_subjects:
                        numsubs_arg = " " + " ".join(map(str, num_subjects))
                        if group_names:
                            groupnmes_arg = \
                                " --group_names " + " ".join(group_names)
                    if version:
                        version_arg = " -v " + version

                    nidmfsl_cmd = [
                        "nidmfsl " + featdir_arg + numsubs_arg +
                        groupnmes_arg + version_arg]
                    print "Running " + str(nidmfsl_cmd)
                    subprocess.check_call(nidmfsl_cmd, shell=True)

                    zipped_dir = os.path.join(
                        data_dir, os.path.basename(data_dir) + ".nidm.zip")

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
