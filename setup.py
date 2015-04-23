#!/usr/bin/env python
from distutils.core import setup
setup(
    name='pyms',
    version='2.0',
    author='Sean O\'Callaghan',
    author_email='spoc@unimelb.edu.au',
    packages=['pyms'],
    url='http://bioinformatics.bio21.unimelb.edu.au/pyms',
    license='LICENSE.txt',
    description='Python Toolkit for Mass Spectrometry',
    long_description=('PyMS is a collection of PyMS libraries'
                      'which can be used to perform peak picking',
                      'alignment by dynamic programming, and'
                      'gap-filling in Gas Chromatography Mass'
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
        "pymzml >= 0.7.5"
        "pycluster >= 1.49"
],
)
