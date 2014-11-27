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

import logging
logger = logging.getLogger(__name__)
# Display log messages in console
logging.basicConfig(level=logging.INFO, format='%(levelname)s - %(message)s')

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
    # The FSL export to NIDM will only be run locally (for now)
    from nidmfsl.fsl_exporter.fsl_exporter import FSLtoNIDMExporter

NIDM_RESULTS_DIR = os.path.join(NIDM_DIR, "nidm", "nidm-results")
TERM_RESULTS_DIR = os.path.join(NIDM_RESULTS_DIR, "terms")
TEST_DIR = os.path.dirname(os.path.abspath(__file__))
TEST_DIR_001 = os.path.join(TEST_DIR, 'example001')
TEST_DIR_002 = os.path.join(TEST_DIR, 'example002')
DATA_DIR_001 = os.path.join(RELPATH, 'test', 'data', 'fmri_one_contrast.feat')
DATA_DIR_002 = os.path.join(RELPATH, 'test', 'data', 'fmri_one_contrast_voxelwise.feat')

path = os.path.join(NIDM_RESULTS_DIR, "test")
sys.path.append(path)


from TestResultDataModel import TestResultDataModel
from TestCommons import *
from CheckConsistency import *

class TestFSLResultDataModel(unittest.TestCase, TestResultDataModel):
    """
    Tests based on the analysis of single-subject fmri fluency data as 
    described at 
    http://fsl.fmrib.ox.ac.uk/fslcourse/lectures/practicals/feat1/index.html 
    but with only *1 contrast: Generation*
    """
    
    @classmethod
    def setUpClass(cls):
        # *** Once for all, run the export and convert provn to ttl
        for test_dir, data_dir in [ (TEST_DIR_001, DATA_DIR_001),
                                    (TEST_DIR_002,DATA_DIR_002)]:

            #  Turtle file obtained with FSL NI-DM export tool
            provn = os.path.join(test_dir, 'FSL_example.provn');
            ttl = os.path.join(test_dir, 'FSL_example.ttl');

            # If test data is available (usually if the test is run locally) then 
            # compute a fresh export
            if os.path.isdir(data_dir):
                logging.debug("Computing NIDM FSL export")

                # Export to NIDM using FSL export tool
                # fslnidm = FSL_NIDM(feat_dir=DATA_DIR_001);
                fslnidm = FSLtoNIDMExporter(feat_dir=data_dir, version="0.2.0")
                fslnidm.parse()
                export_dir = fslnidm.export()
                # Copy provn export to test directory
                shutil.copy(os.path.join(export_dir, 'nidm.provn'), 
                            os.path.join(provn))
                shutil.copy(os.path.join(export_dir, 'nidm.ttl'), 
                            os.path.join(ttl))

            # # Equivalent turtle file converted using the ProvStore API
            # # Local file to save the turtle export (and avoid multiple calls to ProvStore)
            # ttl = os.path.join(test_dir, 'FSL_example.ttl');
            # ttl_url = get_turtle(provn)
            # ttl_fid = open(ttl, 'w')
            # ttl_fid.write(urllib2.urlopen(ttl_url).read())
            # ttl_fid.close()

    def setUp(self):
        TestResultDataModel.setUp(self) 
        self.ttl_001 = os.path.join(TEST_DIR_001, 'FSL_example.ttl');
        self.ttl_002 = os.path.join(TEST_DIR_002, 'FSL_example.ttl');

        # RDF obtained by the FSL export 
        self.graph_001 = Graph()
        # self.ttl_001 = os.path.join(self.test_dir, 'fsl', 'export', 'test01', 'fsl_nidm.ttl');
        self.graph_001.parse(self.ttl_001, format='turtle')
        self.graph_002 = Graph()
        self.graph_002.parse(self.ttl_002, format='turtle')

        # Retreive owl file for NIDM-Results
        self.owl_file = os.path.join(TERM_RESULTS_DIR, 'nidm-results.owl')

        # Move in test dir (storage of prov file)
        # fsl_test_dir = os.path.join(RELPATH, 'test')

    def test01_class_consistency_with_owl(self):
        for graph in [self.graph_001, self.graph_002]:
            my_exception = check_class_names(graph, "FSL example00", owl_file=self.owl_file)

            # FIXME (error message display should be simplified when only one example...)
            if my_exception:
                error_msg = ""
                for unrecognised_class_name, examples in my_exception.items():
                    error_msg += unrecognised_class_name+" (from "+', '.join(examples)+")"
                raise Exception(error_msg)


    def test02_attributes_consistency_with_owl(self):
        for graph in [self.graph_001, self.graph_002]:
            my_exception = check_attributes(graph, "FSL example001", owl_file=self.owl_file)

            # FIXME (error message display should be simplified when only one example...)
            error_msg = ""
            if my_exception[0]:
                for unrecognised_attribute, example_names in my_exception[0].items():
                    error_msg += unrecognised_attribute+" (from "+', '.join(example_names)+")"
            if my_exception[1]:
                for unrecognised_range, example_names in my_exception[1].items():
                    error_msg += unrecognised_range+" (from "+', '.join(example_names)+")"
            if error_msg:
                raise Exception(error_msg)

    # FIXME: If terms PR is accepted then these tests should be moved to TestResultDataModel.py
    def test03_ex1_auditory_singlesub_full_graph(self):
        """
        Test03: Comparing that the ttl file generated by FSL and the expected 
        ttl file (generated manually) are identical
        """
        ground_truth_dir = os.path.join(NIDM_RESULTS_DIR,'fsl', "example001")
        ground_truth_ttl = os.path.join(ground_truth_dir, 'fsl_nidm.ttl');
        logging.info("Ground truth ttl: "+ground_truth_ttl)

        # RDF obtained by the ground truth export
        gt = Graph()
        gt.parse(ground_truth_ttl, format='turtle')

        self.compare_full_graphs(gt, self.graph_001)

        if self.my_execption:
            raise Exception(self.my_execption)

    def test_voxelwise_threshold_fwe05(self):
        """
        Check that minimal set of relations needed to describe FWE p<0.05 
        voxel-wise thresholding is present.
        """
        ground_truth_dir = os.path.join(NIDM_RESULTS_DIR,'fsl', "example002")
        ground_truth_ttl = os.path.join(ground_truth_dir, 'fsl_nidm.ttl');
        logging.info("Ground truth ttl: "+ground_truth_ttl)

        # RDF obtained by the ground truth export
        gt = Graph()
        gt.parse(ground_truth_ttl, format='turtle')

        self.compare_full_graphs(gt, self.graph_002, True)

        if self.my_execption:
            raise Exception(self.my_execption)


    @classmethod
    def tearDownClass(cls):
        # Delete temporarily written out ttl file
        os.path.join(TEST_DIR_001, 'FSL_example.ttl');
        # os.remove(ttl_001)

    # '''Test02: Test availability of attributes needed to perform a meta-analysis as specified in use-case *1* at: http://wiki.incf.org/mediawiki/index.php/Queries'''
    # def test02_metaanalysis_usecase1(self):
    #     prefixInfo = """
    #     prefix prov: <http://www.w3.org/ns/prov#>
    #     prefix fsl: <http://www.fil.ion.ucl.ac.uk/fsl/ns/#>
    #     prefix nidm: <http://nidm.nidash.org/>

    #     """
    #     # Look for:
    #     # - "location" of "Contrast map",
    #     # - "location" of "Contrast variance map",
    #     # - "prov:type" in "nidm" namespace of the analysis software.
    #     query = prefixInfo+"""
    #     SELECT ?cfile ?efile ?stype WHERE {
    #      ?aid a fsl:contrast ;
    #           prov:wasAssociatedWith ?sid.
    #      ?sid a prov:Agent;
    #           a prov:SoftwareAgent;
    #           a ?stype . 
    #      FILTER regex(str(?stype), "nidm") 
    #      ?cid a nidm:contrastMap ;
    #           prov:wasGeneratedBy ?aid ;
    #           prov:atLocation ?cfile .
    #      ?eid a nidm:contrastStandardErrorMap ;
    #           prov:wasGeneratedBy ?aid ;
    #           prov:atLocation ?efile .
    #     }
    #     """

    #     if not self.successful_retreive(self.graph_001.query(query), 'ContrastMap and ContrastStandardErrorMap'):
    #         raise Exception(self.my_execption)

    # '''Test03: Test availability of attributes needed to perform a meta-analysis as specified in use-case *2* at: http://wiki.incf.org/mediawiki/index.php/Queries'''
    # def test03_metaanalysis_usecase2(self):
    #     prefixInfo = """
    #     prefix prov: <http://www.w3.org/ns/prov#>
    #     prefix fsl: <http://www.fil.ion.ucl.ac.uk/fsl/ns/#>
    #     prefix nidm: <http://nidm.nidash.org/>

    #     """

    #     # Look for:
    #     # - "location" of "Contrast map",
    #     # - "prov:type" in "nidm" namespace of the analysis software.
    #     query = prefixInfo+"""
    #     SELECT ?cfile ?efile ?stype WHERE {
    #      ?aid a fsl:contrast ;
    #           prov:wasAssociatedWith ?sid.
    #      ?sid a prov:Agent;
    #           a prov:SoftwareAgent;
    #           a ?stype . 
    #      FILTER regex(str(?stype), "nidm") 
    #      ?cid a nidm:contrastMap ;
    #           prov:wasGeneratedBy ?aid ;
    #           prov:atLocation ?cfile .
    #     }
    #     """

    #     if not self.successful_retreive(self.graph_001.query(query), 'ContrastMap and ContrastStandardErrorMap'):
    #         raise Exception(self.my_execption)

    # '''Test04: Test availability of attributes needed to perform a meta-analysis as specified in use-case *3* at: http://wiki.incf.org/mediawiki/index.php/Queries'''
    # def test04_metaanalysis_usecase3(self):
    #     prefixInfo = """
    #     prefix prov: <http://www.w3.org/ns/prov#>
    #     prefix fsl: <http://www.fil.ion.ucl.ac.uk/fsl/ns/#>
    #     prefix nidm: <http://nidm.nidash.org/>

    #     """

    #     # Look for:
    #     # - "location" of "Statistical map",
    #     # - "nidm:errorDegreesOfFreedom" in "Statistical map".
    #     query = prefixInfo+"""
    #     SELECT ?sfile ?dof WHERE {
    #      ?sid a nidm:statisticalMap ;
    #           prov:atLocation ?sfile ;
    #           nidm:errorDegreesOfFreedom ?dof .
    #     }
    #     """

    #     if not self.successful_retreive(self.graph_001.query(query), 'ContrastMap and ContrastStandardErrorMap'):
    #         raise Exception(self.my_execption)

    # '''Test05: Test availability of attributes needed to perform a meta-analysis as specified in use-case *4* at: http://wiki.incf.org/mediawiki/index.php/Queries'''
    # def test05_metaanalysis_usecase4(self):
    #     prefixInfo = """
    #     prefix prov: <http://www.w3.org/ns/prov#>
    #     prefix fsl: <http://www.fil.ion.ucl.ac.uk/fsl/ns/#>
    #     prefix nidm: <http://nidm.nidash.org/>

    #     """

    #     # Look for:
    #     # - For each "Peak" "equivZStat" and"coordinate1" (and optionally "coordinate2" and "coordinate3"),
    #     # - "clusterSizeInVoxels" of "height threshold"
    #     # - "value" of "extent threshold"
    #     query = prefixInfo+"""
    #     SELECT ?equivz ?coord1 ?coord2 ?coord3 ?ethresh ?hthresh WHERE {
    #      ?pid a fsl:peakLevelStatistic ;
    #         prov:atLocation ?cid ;
    #         nidm:equivalentZStatistic ?equivz ;
    #         prov:wasDerivedFrom ?clid .
    #      ?cid a nidm:coordinate;
    #         nidm:coordinate1 ?coord1 .
    #         OPTIONAL { ?cid nidm:coordinate2 ?coord2 }
    #         OPTIONAL { ?cid nidm:coordinate3 ?coord3 }
    #      ?iid a nidm:inference .
    #      ?esid a fsl:excursionSet;
    #         prov:wasGeneratedBy ?iid .
    #      ?setid a fsl:setLevelStatistic;
    #         prov:wasDerivedFrom ?esid .
    #      ?clid a fsl:clusterLevelStatistic;
    #         prov:wasDerivedFrom ?setid .
    #      ?tid a nidm:extentThreshold ;
    #         nidm:clusterSizeInVoxels ?ethresh .
    #      ?htid a nidm:heightThreshold ;
    #         prov:value ?hthresh .
    #     }
    #     """

    #     if not self.successful_retreive(self.graph_001.query(query), 'ContrastMap and ContrastStandardErrorMap'):
    #         raise Exception(self.my_execption)

if __name__ == '__main__':
    unittest.main()

