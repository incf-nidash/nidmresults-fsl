#!/usr/bin/env python
from FSLparser import FSL_NIDM

# pathTofslResDir = '/Users/cmaumet/Data/fsl_practicals/fsl_course_data/fmri/fmri_fluency/fmri.feat'
pathTofslResDir = '/Users/cmaumet/Data/fsl_practicals/fsl_course_data/fmri/av/fmri.feat'

fslnidm = FSL_NIDM(featDir=pathTofslResDir);
# fslnidm.add_report_file('./data/report_poststats.html')
# fslnidm.add_zstat_file('./data/cluster_zstat1-short.html','./data/cluster_zstat1_std-short.html')
fslnidm.save_prov_to_files()
