#!/usr/bin/python
"""
Export neuroimaging results created with AFNI.

@author: R. Reynolds/C. Maumet
originally: @author: Camille Maumet <c.m.j.maumet@warwick.ac.uk>
"""

import sys
import os
from afni_exporter.afni_exporter import AFNItoNIDMExporter

if __name__ == "__main__":
    # Remove first argument (script name)
    num_args = len(sys.argv)-1
    sys.argv.pop(0)
    args = sys.argv

    usage = "Usage: python nidm-results_afni.py path/to/dataset clustsim"

    if num_args != 2:
        raise Exception(usage)

    dset = args[0]
    cset = args[1]
    # p_unc = args[2]	rcr - add option processing
    # p_cor = args[3]

    # check for existance given various extensions
    #if not os.path.isdir(feat_dir):
    #    raise Exception("Unknown directory: "+str(feat_dir))

    afninidm = AFNItoNIDMExporter(dset=dset, csim_dset=cset,
		 		  # p_uncor=p_unc, p_cor=p_cor,
				  nidm_ver="0.2.0")
    afninidm.parse()
    export_dir = afninidm.export()

    print 'NIDM export available at: '+str(export_dir)
