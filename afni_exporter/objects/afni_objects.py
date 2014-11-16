"""
AFNI-specific classes and classes overloaded to add AFNI-specific attributes.

@author: Rick Reynolds/Camille Maumet
@copyright:
"""
from exporter.objects.generic import NIDMObject
from exporter.objects.constants import *
import logging
import uuid

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.DEBUG)

class Software(NIDMObject):
    # FIXME software should be generic and then overloaded

    """
    Class representing a Software entity.
    """

    def __init__(self, version):
        super(Software, self).__init__()
        self.version = version
        self.name = "AFNI"
        self.id = NIIRI[str(uuid.uuid4())]
        # Retreive AFNI version from afni -ver

    def export(self):
        """
        Create prov entities and activities.
        """
        self.p.agent(self.id, 
            other_attributes=((PROV['type'], NIDM[self.name]), 
                            (PROV['type'], PROV['SoftwareAgent']),
                            (PROV['label'],self.name),
                            (NIDM['softwareVersion'],self.version)))

        return self.p
