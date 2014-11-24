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
    from fsl_exporter.fsl_exporter import FSLtoNIDMExporter

NIDM_RESULTS_DIR = os.path.join(NIDM_DIR, "nidm", "nidm-results")
TERM_RESULTS_DIR = os.path.join(NIDM_RESULTS_DIR, "terms")
TEST_FOLDER = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'example001')
DATA_DIR = os.path.join(RELPATH, 'test', 'data', 'fmri.feat')

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
        
        #  Turtle file obtained with FSL NI-DM export tool
        fsl_export_provn = os.path.join(TEST_FOLDER, 'FSL_example.provn');

        # If test data is available (usually if the test is run locally) then 
        # compute a fresh export
        if os.path.isdir(DATA_DIR):
            logging.debug("Computing NIDM FSL export")

            # Export to NIDM using FSL export tool
            # fslnidm = FSL_NIDM(feat_dir=DATA_DIR);
            fslnidm = FSLtoNIDMExporter(feat_dir=DATA_DIR, version="0.2.0")
            fslnidm.parse()
            export_dir = fslnidm.export()

            # Copy provn export to test directory
            shutil.copy(os.path.join(export_dir, 'nidm.provn'), 
                        os.path.join(fsl_export_provn))

        # Equivalent turtle file converted using the ProvStore API
        # Local file to save the turtle export (and avoid multiple calls to ProvStore)
        fsl_export_ttl = os.path.join(TEST_FOLDER, 'FSL_example.ttl');
        fsl_export_ttl_url = get_turtle(fsl_export_provn)
        ttl_fid = open(fsl_export_ttl, 'w')
        ttl_fid.write(urllib2.urlopen(fsl_export_ttl_url).read())
        ttl_fid.close()

    def setUp(self):
        TestResultDataModel.setUp(self) 
        self.ground_truth_dir = os.path.join(NIDM_RESULTS_DIR,'fsl', 'example001')
        self.fsl_export_ttl = os.path.join(TEST_FOLDER, 'FSL_example.ttl');

        # RDF obtained by the FSL export 
        self.fslexport = Graph()
        # self.fsl_export_ttl = os.path.join(self.test_dir, 'fsl', 'export', 'test01', 'fsl_nidm.ttl');
        self.fslexport.parse(self.fsl_export_ttl, format='turtle')

        # Retreive owl file for NIDM-Results
        self.owl_file = os.path.join(TERM_RESULTS_DIR, 'nidm-results.owl')

        # Move in test dir (storage of prov file)
        # fsl_test_dir = os.path.join(RELPATH, 'test')

    def test01_class_consistency_with_owl(self):
        my_exception = check_class_names(self.fslexport, "FSL example001", owl_file=self.owl_file)

        # FIXME (error message display should be simplified when only one example...)
        if my_exception:
            error_msg = ""
            for unrecognised_class_name, examples in my_exception.items():
                error_msg += unrecognised_class_name+" (from "+', '.join(examples)+")"
            raise Exception(error_msg)


    def test02_attributes_consistency_with_owl(self):
        my_exception = check_attributes(self.fslexport, "FSL example001", owl_file=self.owl_file)

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
        ground_truth_ttl = os.path.join(self.ground_truth_dir, 'fsl_nidm.ttl');
        logging.info("Ground truth ttl: "+ground_truth_ttl)

        # RDF obtained by the ground truth export
        gt = Graph()
        gt.parse(ground_truth_ttl, format='turtle')

        self.compare_full_graphs(gt, self.fslexport)

        if self.my_execption:
            raise Exception(self.my_execption)

    @classmethod
    def tearDownClass(cls):
        # Delete temporarily written out ttl file
        os.path.join(TEST_FOLDER, 'FSL_example.ttl');
        # os.remove(fsl_export_ttl)

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

    #     if not self.successful_retreive(self.fslexport.query(query), 'ContrastMap and ContrastStandardErrorMap'):
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

    #     if not self.successful_retreive(self.fslexport.query(query), 'ContrastMap and ContrastStandardErrorMap'):
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

    #     if not self.successful_retreive(self.fslexport.query(query), 'ContrastMap and ContrastStandardErrorMap'):
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

    #     if not self.successful_retreive(self.fslexport.query(query), 'ContrastMap and ContrastStandardErrorMap'):
    #         raise Exception(self.my_execption)

if __name__ == '__main__':
    unittest.main()

