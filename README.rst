NIDM-Results export for FSL
===========================

Export of FSL FEAT statistical results using the NeuroImaging Data Model
(`NIDM Results`_).

Install
-------

::

    $ pip install nidmfsl

Usage
-----

::

$ nidmfsl [-h] [-g [GROUP_NAMES [GROUP_NAMES ...]]] [-o OUTPUT_NAME] [-d] [-v [VERSION]] feat_dir [numsubjects [numsubjects ...]]


Requirements
------------

-  `rdflib`_
-  `prov`_
-  `nibabel`_
-  `numpy`_
-  `nidmresults`_

.. _NIDM Results: http://nidm.nidash.org/specs/nidm-results.html
.. _prov: https://github.com/trungdong/prov
.. _nibabel: http://nipy.org/nibabel/
.. _numpy: http://www.numpy.org/
.. _nidmresults: https://github.com/incf-nidash/nidmresults/
.. _rdflib: http://rdflib.readthedocs.org/en/latest/