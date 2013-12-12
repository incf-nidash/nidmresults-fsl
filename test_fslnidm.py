#!/usr/bin/env python
from FSLparser import FSL_NIDM

fslnidm = FSL_NIDM();
fslnidm.add_report_file('./data/report_poststats.html')
fslnidm.add_zstat_file('./data/cluster_zstat1-short.html','./data/cluster_zstat1_std-short.html')
fslnidm.save_prov_to_files()
