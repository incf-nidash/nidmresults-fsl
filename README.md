NIDM-Results export for FSL
===========================

Export of FSL FEAT statistical results using the NeuroImaging Data Model ([NIDM Results]).

Install
-------

    $ pip install nidmfsl

Usage
-----

    $ nidmfsl [-h] [-g [GROUP_NAMES [GROUP_NAMES ...]]] [-o OUTPUT_NAME] [-d] [-v [VERSION]] feat_dir [numsubjects [numsubjects ...]]

Requirements
------------

-   [rdflib]
-   [prov]
-   [nibabel]
-   [numpy]
-   [nidmresults]

  [NIDM Results]: http://nidm.nidash.org/specs/nidm-results.html
  [rdflib]: http://rdflib.readthedocs.org/en/latest/
  [prov]: https://github.com/trungdong/prov
  [nibabel]: http://nipy.org/nibabel/
  [numpy]: http://www.numpy.org/
  [nidmresults]: https://github.com/incf-nidash/nidmresults/