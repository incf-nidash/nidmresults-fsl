"""
FSL-specific classes and classes overloaded to add FSL-specific attributes.

@author: Camille Maumet <c.m.j.maumet@warwick.ac.uk>
@copyright: University of Warwick 2013-2014
"""
from exporter.objects.generic import NIDMObject
from exporter.objects.constants import *

class Software(NIDMObject):
    # FIXME software should be generic and then overloaded

    """
    Class representing a Software entity.
    """

    def __init__(self, feat_version):
        super(Software, self).__init__()
        self.feat_version = feat_version
        self.name = "FSL"
        self.id = NIIRI['software_id']

    def export(self):
        """
        Create prov entities and activities.
        """
        self.p.agent(self.id, 
            other_attributes=((PROV['type'], NIDM[self.name]), 
                            (PROV['type'], PROV['SoftwareAgent']),
                            (PROV['label'],self.name),
                            (FSL['featVersion'], self.feat_version)))


        return self.p