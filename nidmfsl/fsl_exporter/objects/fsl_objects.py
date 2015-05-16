"""
FSL-specific classes and classes overloaded to add FSL-specific attributes.

@author: Camille Maumet <c.m.j.maumet@warwick.ac.uk>
@copyright: University of Warwick 2013-2014
"""
from nidmresults.objects.generic import NIDMObject
from nidmresults.objects.constants import *
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
        self.type = NLX_FSL
        self.prov_type = PROV['Agent']
        # Retreive FSL version from feat version
        # (cf. https://github.com/incf-nidash/nidm-results_fsl/issues/3)
        if feat_version == "6.00":
            self.version = "5.0.x"
        elif feat_version == "5.98":
            self.version = "4.1.x"
        elif feat_version == "5.92":
            self.version = "4.0.x"
        elif feat_version == "5.91":
            self.version = "4.0.1"
        elif feat_version == "5.90":
            self.version = "4.0"
        elif feat_version == "5.61":
            self.version = "3.3"
        elif feat_version == "5.4":
            self.version = "3.2.1"
        elif feat_version == "5.1":
            self.version = "3.1.1.x"
        else:
            logging.debug("FSL version unknow for feat version: \"" +
                          feat_version + "\"")
            self.version = "unknown"

    def export(self):
        """
        Create prov entities and activities.
        """
        self.add_attributes((
            (PROV['type'], NLX_FSL),
            (PROV['type'], PROV['SoftwareAgent']),
            (PROV['label'], self.name),
            (NIDM_SOFTWARE_VERSION, self.version),
            (FSL_FEAT_VERSION, self.feat_version)))

        return self.p
