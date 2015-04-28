#!/usr/bin/env python
try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

try:
    from setuptools import find_packages
except:
    pass

setup(
    name='pyms',
    version='2.1',
    author='Sean O\'Callaghan',
    author_email='spoc@unimelb.edu.au',
    packages=find_packages(),
    url='https://github.com/ma-bio21/pyms',
    download_url = 'https://github.com/ma-bio21/pyms/tarball/0.1',
    license='LICENSE.txt',
    description='Python Toolkit for Mass Spectrometry',
    long_description=('PyMS is a collection of PyMS libraries '
                      'which can be used to perform peak picking ',
                      'alignment by dynamic programming, and '
                      'gap-filling in Gas Chromatography Mass '
                      'Spectrometry experiments'),
    install_requires=[
        # PyMS also requires both scipy and numpy. Installation
        # of scipy and numpy should be performed before PyMS 
        # installation, as these packages do not play
        # well with automated package installers
        # In addition, matplotlib is required if graphical output is 
        # desired. 
        # A good solution which provides all of these dependencies is
        # Anaconda
        #"numpy >= 1.7.1",
        #"scipy >= 0.12.0",
        "pymzml >= 0.7.5",
        "pycluster >= 1.49",
        ],
)
