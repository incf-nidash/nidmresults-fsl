NIDM-Results for FSL
====================

Export mass-univariate neuroimaging results computed in FSL (using FEAT)
as NIDM-Results packs.

A *NIDM-Results pack* is a compressed file containing a NIDM-Results
serialization and some or all of the referenced image data files in
compliance with `NIDM-Results specification`_.

Usage
'''''

::

usage: nidmfsl [-h] [-g [GROUP_NAMES [GROUP_NAMES ...]]] [-o OUTPUT_NAME] [-d]
               [-v [NIDM_VERSION]] [--version]
               feat_dir [numsubjects [numsubjects ...]]

NIDM-Results exporter for FSL Feat.

positional arguments:
  feat_dir              Path to feat directory.
  numsubjects           Number of subjects per group.

optional arguments:
  -h, --help            show this help message and exit
  -g [GROUP_NAMES [GROUP_NAMES ...]], --group_names [GROUP_NAMES [GROUP_NAMES ...]]
                        Label for each group.
  -o OUTPUT_NAME, --output_name OUTPUT_NAME
                        Name of the output. A ".nidm.zip" or ".nidm" (when -d
                        is used) suffix will be appended.
  -d, --directory-output
                        Produces a .nidm directory rather than a .nidm.zip
                        file.
  -v [NIDM_VERSION], --nidm_version [NIDM_VERSION]
                        NIDM-Results version to use (default: latest).
  --version             show program's version number and exit

Requirements
''''''''''''

-  `nidmresults`_

Installation
''''''''''''

::

        pip install nidmfsl

.. _NIDM-Results specification: http://nidm.nidash.org/specs/nidm-results.html
.. _nidmresults: http://pypi.python.org/pypi/nidmresults