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

       nidmfsl [-h] [-g [GROUP_NAMES [GROUP_NAMES ...]]] [-o OUTPUT_NAME] [-d]
                   [-v [NIDM_VERSION]] [--version]
                   feat_dir [numsubjects [numsubjects ...]]

Requirements
''''''''''''

-  `nidmresults`_

Installation
''''''''''''

::

        pip install nidmfsl

.. _NIDM-Results specification: http://nidm.nidash.org/specs/nidm-results.html
.. _nidmresults: http://pypi.python.org/pypi/nidmresults