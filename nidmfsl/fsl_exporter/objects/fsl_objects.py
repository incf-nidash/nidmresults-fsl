"""
FSL-specific classes and classes overloaded to add FSL-specific attributes.

@author: Camille Maumet <c.m.j.maumet@warwick.ac.uk>
@copyright: University of Warwick 2013-2014
"""
from nidmresults.objects.generic import ExporterSoftware, NeuroimagingSoftware
from nidmresults.objects.constants import *
import nidmfsl
import logging
import uuid

logger = logging.getLogger(__name__)


class FSLNeuroimagingSoftware(NeuroimagingSoftware):
    """
    Class representing a Software entity.
    """

    def __init__(self, feat_version):
        self.id = NIIRI[str(uuid.uuid4())]
        self.feat_version = feat_version
        # Retreive FSL version from feat version
        # (cf. https://github.com/incf-nidash/nidm-results_fsl/issues/3)
        feat_to_fsl_versions = {
            "6.00":   "5.0.x",
            "5.98":   "4.1.x",
            "5.92":   "4.0.x",
            "5.91":   "4.0.1",
            "5.90":   "4.0",
            "5.61":   "3.3",
            "5.4":    "3.2.1",
            "5.1":    "3.1.1.x"
        }

        if feat_version in feat_to_fsl_versions:
            version = feat_to_fsl_versions[feat_version]
        else:
            logging.debug("FSL version unknow for feat version: \"" +
                          feat_version + "\"")
            version = "unknown"

        super(FSLNeuroimagingSoftware, self).__init__("fsl", version)

    def export(self, nidm_version):
        """
        Create prov entities and activities.
        """
        super(FSLNeuroimagingSoftware, self).export(nidm_version)
        self.add_attributes([(FSL_FEAT_VERSION, self.feat_version)])

        return self.p


class FSLExporterSoftware(ExporterSoftware):
    """
    Class representing a Software entity.
    """

    def __init__(self):
        self.id = NIIRI[str(uuid.uuid4())]
        super(FSLExporterSoftware, self).__init__(
            NIDM_FSL, nidmfsl.__version__)

    def export(self, nidm_version):
        """
        Create prov entities and activities.
        """
        super(FSLExporterSoftware, self).export(nidm_version)

        return self.p
