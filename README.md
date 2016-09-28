
# NIDM-Results for FSL

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
```
    pip install nidmfsl
```

##### Compatible with 
FSL version 5.0.9
