# pyms
The original source is a python package called pyms that was developed and published in 2012.
See here: https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-13-115

The code was stored on google code, you can see the archive here: http://code.google.com/p/pyms/

The upstream repository here on github from @ma-bio21 seems more up to date than the code found in the archive, 
most notably the implementation of using mpi4py in loading ANDI files. There are other differences though, 
some bug fixes and added functionality. 

Also a note: in the setup.py there is a requirement for the Pycluster library and it is used in the 
files found in the ./Peak/List/DPA directory, but installing the library from pip failed on my machine 
(Linux Mint 18 Sarah with python 2.7.12 installed).

A diff of the two directories list the following files as different:
- ./Deconvolution/BillerBiemann/Function.py 
- ./Display/Class.py 
- ./Display/Class.pyc 
- ./Experiment/IO.py 
- ./Gapfill/Class.py 
- ./Gapfill/Function.py 
- ./Gapfill/__init__.py 
- ./GCMS/Class.py 
- ./GCMS/IO/ANDI/Function.py 
- ./GCMS/IO/MZML/Function.py 
- ./GCMS/IO/MZML/__init__.py 
- ./Peak/Class.py 
- ./Peak/List/DPA/Class_old.py 
- ./Peak/List/DPA/Class.py 
- ./Peak/List/DPA/Function_new.py 
- ./Peak/List/DPA/Function.py 
- ./Utils/DP.py 
