'''Generic class to export software packages (FSL, AFNI, ...) results into NI-DM

@author: Camille Maumet <c.m.j.maumet@warwick.ac.uk>
@copyright: University of Warwick 2013-2014
'''

from prov.model import ProvBundle, ProvDocument
import os
import datetime
from objects.constants import *
from objects.modelfitting import *
from objects.contrast import *
from objects.inference import *

''' Parse an FSL result directory to extract the pieces information stored in NIDM-Results
'''
class NIDMExporter():
    
    def parse(self):
        # Software agent
        self.software = self.find_software()

        # Model Parameters Estimation activity and entities      
        self.model_fittings = self.find_model_fitting()
        # Contrast Estimation activity and entities
        self.contrasts = self.find_contrasts()
        # Inference activity and entities
        self.inferences = self.find_inferences()
                
        # Create namespaces
        self.provDocument = ProvDocument();
        self.add_namespaces()


    def get_model_fitting(self, mf_id):
        for model_fitting in self.model_fittings:
            if model_fitting.activity.id == mf_id:
                return model_fitting
        raise Exception("Model fitting activity with id: "+str(mf_id)+" not found.")

    def get_contrast(self, con_id):
        for contrasts in self.contrasts.values():
            for contrast in contrasts:
                if contrast.estimation.id == con_id:
                    return contrast
        raise Exception("Contrast activity with id: "+str(con_id)+" not found.")

    def export(self):
        self.create_bundle(self.version)

        self.provBundle.update(self.software.export())

        for model_fitting in self.model_fittings:
            self.provBundle.update(model_fitting.export())
            self.provBundle.wasAssociatedWith(model_fitting.activity.id, self.software.id)

        for (model_fitting_id, pe_ids), contrasts in self.contrasts.items():
            model_fitting = self.get_model_fitting(model_fitting_id)
            for contrast in contrasts:
                self.provBundle.update(contrast.export())
                self.provBundle.used(contrast.estimation.id, model_fitting.rms_map.id)
                self.provBundle.used(contrast.estimation.id, model_fitting.mask_map.id)
                self.provBundle.wasAssociatedWith(contrast.estimation.id, self.software.id)
                for pe_id in pe_ids:
                    self.provBundle.used(contrast.estimation.id, pe_id)

        for contrast_id, inferences in self.inferences.items():
            contrast = self.get_contrast(contrast_id)
            for inference in inferences:
                self.provBundle.update(inference.export())
                if contrast.z_stat_map:
                    used_id = contrast.z_stat_map.id
                else:
                    used_id = contrast.stat_map.id
                self.provBundle.used(inference.id, used_id)
                self.provBundle.wasAssociatedWith(inference.id, self.software.id)

        # for element_to_export in self.elements_to_export:
        #     other = element_to_export.export()
        #     self.provBundle.update(other)

        self.save_prov_to_files()

    def add_namespaces(self):
        self.provDocument.add_namespace("neurolex", "http://neurolex.org/wiki/")
        self.provDocument.add_namespace(FSL)
        self.provDocument.add_namespace(NIDM)
        self.provDocument.add_namespace(NIIRI)
        self.provDocument.add_namespace(CRYPTO)
        self.provDocument.add_namespace(DCT)

    def create_bundle(self, version):
        software_lc = self.software.name.lower()
        software_uc = self.software.name.upper()

        self.provBundle = ProvBundle(identifier=NIIRI[software_lc+'_results_id'])

        self.provDocument.entity(NIIRI[software_lc+'_results_id'], 
            other_attributes=( (PROV['type'], PROV['Bundle'],), 
                               (PROV['label'],software_uc+" Results" ),
                               (NIDM['objectModel'],NIDM[software_uc+'Results']),
                               (NIDM['version'], version))
            )

        self.provDocument.wasGeneratedBy(NIIRI[software_lc+'_results_id'], 
            time=str(datetime.datetime.now().time()))

    # Method used for model estimation is directly be inferred from the error term
    def _get_model_parameters_estimations(self, error_model):
        model_parameters_estimations = list()
        if error_model.dependance == INDEPEDENT_CORR:
            if error_model.variance_homo:
                estimation_method = ESTIMATION_OLS
            else:
                estimation_method = ESTIMATION_WLS
        else:
            estimation_method = ESTIMATION_GLS

        mpe = ModelParametersEstimation(estimation_method, self.software.id)

        return mpe

    def save_prov_to_files(self, showattributes=False):
        self.provDocument.add_bundle(self.provBundle)

        suffixName = ''
        if showattributes is False:
            suffixName = '_without_attributes'

        # jsondata = self.provBundle.get_provjson(indent=4)
        # JSONfile = open(os.path.join(self.export_dir, 'nidm.json'), 'w');
        # JSONfile.write(jsondata)
        PROVNfile = open(os.path.join(self.export_dir, 'nidm.provn'), 'w');
        PROVNfile.write(self.provDocument.get_provn(4))