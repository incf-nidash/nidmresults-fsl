from setuptools import setup, find_packages

readme = open('README.rst').read()

reqs = [line.strip() for line in open('requirements.txt').readlines()]
requirements = list(filter(None, reqs))

setup(
    name="nidmfsl",
    version="2.1.0",
    author="Camille Maumet",
    author_email="c.m.j.maumet@warwick.ac.uk",
    description=("Export of FSL statistical results using NIDM\
 as specified at http://nidm.nidash.org/specs/nidm-results.html."),
    license = "BSD",
    keywords = "Prov, NIDM, Provenance",
    scripts=['bin/nidmfsl'],
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
        "Programming Language :: Python :: 3",
    ],
    install_requires=requirements
)
