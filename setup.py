from setuptools import setup, find_packages

readme = open('README.md').read()

requirements = [
    'prov',
    'nibabel',
    'numpy'    
]

setup(
    name = "nidmfsl",
    version = "0.1.0",
    author = "Camille Maumet",
    author_email = "c.m.j.maumet@warwick.ac.uk",
    description = ("Export of FSL statistical results using NIDM"
        " as specified at http://nidm.nidash.org/specs/nidm-results.html."),
    license = "BSD",
    keywords = "Prov, NIDM, Provenance",
    scripts=['bin/nidm-results_fsl'],
    # packages=['nidmfsl', 'test'],
    packages=find_packages(),
    package_dir={
        'prov': 'prov'
    },
    long_description=readme,
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Topic :: Utilities",
        "License :: OSI Approved :: MIT License",
    ],
)