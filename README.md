
# NIDM-Results for FSL

[![Build Status](https://travis-ci.org/incf-nidash/nidmresults-fsl.svg?branch=master)](https://travis-ci.org/incf-nidash/nidmresults-fsl)

Export mass-univariate neuroimaging results computed in FSL (using FEAT) as NIDM-Results packs.

A *NIDM-Results pack* is a compressed file containing a NIDM-Results serialization and some or all of the referenced image data files in compliance with [NIDM-Results specification](http://nidm.nidash.org/specs/nidm-results.html).

##### Usage
```
usage: nidmfsl [-h] [-g GROUP_NAME NUM_SUBJECTS] [-o OUTPUT_NAME] [-d]
               [-n NIDM_VERSION] [--version]
               feat_dir

NIDM-Results exporter for FSL Feat.

positional arguments:
  feat_dir              Path to feat directory.

optional arguments:
  -h, --help            show this help message and exit
  -g GROUP_NAME NUM_SUBJECTS, --group GROUP_NAME NUM_SUBJECTS
                        Group label followed by number of subjects
  -o OUTPUT_NAME, --output_name OUTPUT_NAME
                        Name of the output. A ".nidm.zip" or ".nidm" (when -d
                        is used) suffix will be appended.
  -d, --directory-output
                        Produces a .nidm directory rather than a .nidm.zip
                        file.
  -n NIDM_VERSION, --nidm_version NIDM_VERSION
                        NIDM-Results version to use (default: latest).
  --version             show program's version number and exit
```


##### Installation

To install, run the below command in the bash terminal.
```
    pip install nidmfsl
```

##### Compatible with 
FSL version 5.0.9

##### Testing

###### Requirements

To run the tests for this repository, the following must be installed

- Git LFS. Installation instructions for Git LFS can be found [here](https://git-lfs.github.com/).

- The python packages `vcr` and `ddt`. These can be installed using the below commands in the bash terminal:

```
pip install vcrpy
pip install ddt
```

In addition, the test data must also be downloaded from the `nidmresults-examples` repository to `<path_to_this_repository>/test/data/nidmresults-examples`.

```
git lfs clone https://github.com/incf-nidash/nidmresults-examples.git <path_to_this_repository>/test/data/nidmresults-examples
```

###### Running the tests

The below command can be used to generate the test cases.
```
python <path_to_this_repository>/test/export_test_battery.py
```

Folowing this, the test cases can be verified against the ground truth provided in the `nidmresults-examples` repository using the below command.
```
cd <path_to_this_repository>/test/
python -m unittest discover
```
