"""
FSL-specific classes and classes overloaded to add FSL-specific attributes.

@author: Camille Maumet <c.m.j.maumet@warwick.ac.uk>
@copyright: University of Warwick 2013-2014
"""
from nidmfsl.exporter.objects.generic import NIDMObject
from nidmfsl.exporter.objects.constants import *
import logging
import uuid

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.DEBUG)

class Software(NIDMObject):
    # FIXME software should be generic and then overloaded

    """
    Class representing a Software entity.
    """

    def __init__(self, feat_version):
        super(Software, self).__init__()
        self.feat_version = feat_version
        self.name = "FSL"
        self.id = NIIRI[str(uuid.uuid4())]
        # Retreive FSL version from feat version 
        # (cf. https://github.com/incf-nidash/nidm-results_fsl/issues/3)
        if feat_version == "6.00":
            self.version = "fsl-5_0_x"
        elif feat_version == "5.98":
            self.version = "fsl-4_1_x"
        elif feat_version == "5.92":
            self.version = "fsl-4_0_x"
        elif feat_version == "5.91":
            self.version = "fsl-4_0_1"
        elif feat_version == "5.90":
            self.version = "fsl-4_0"
        elif feat_version == "5.61":
            self.version = "fsl-3_3"
        elif feat_version == "5.4":
            self.version = "fsl-3_2_1"
        elif feat_version == "5.1":
            self.version = "fsl-3_1_1_x"
        else:
            logging.debug("FSL version unknow for feat version: \""+\
                feat_version+"\"")
            self.version = "fsl-unknown"

    def export(self):
        """
        Create prov entities and activities.
        """
        self.p.agent(self.id, 
            other_attributes=((PROV['type'], NIDM[self.name]), 
                            (PROV['type'], PROV['SoftwareAgent']),
                            (PROV['label'],self.name),
                            (NIDM['softwareVersion'],self.version),
                            (FSL['featVersion'], self.feat_version)))

        return self.p