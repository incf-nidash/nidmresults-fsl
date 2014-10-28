from exporter.objects.generic import NIDMObject
from exporter.objects.constants import *

# NIDM objects in FSL namespace or having FSL-specific attributes

class Software(NIDMObject):
    def __init__(self, feat_version):
        super(Software, self).__init__()
        self.feat_version = feat_version
        self.name = "FSL"
        self.id = NIIRI['software_id']

    def export(self):
        self.p.agent(self.id, 
            other_attributes=((PROV['type'], NIDM[self.name]), 
                            (PROV['type'], PROV['SoftwareAgent']),
                            (PROV['label'],self.name),
                            (FSL['featVersion'], self.feat_version)))


        return self.p