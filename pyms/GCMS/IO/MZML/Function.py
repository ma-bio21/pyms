"""
Functions for reading manufacturer specific ANDI-MS data files
"""

 #############################################################################
 #                                                                           #
 #    PyMS software for processing of metabolomic mass-spectrometry data     #
 #    Copyright (C) 2005-2012 Vladimir Likic                                 #
 #                                                                           #
 #    This program is free software; you can redistribute it and/or modify   #
 #    it under the terms of the GNU General Public License version 2 as      #
 #    published by the Free Software Foundation.                             #
 #                                                                           #
 #    This program is distributed in the hope that it will be useful,        #
 #    but WITHOUT ANY WARRANTY; without even the implied warranty of         #
 #    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the          #
 #    GNU General Public License for more details.                           #
 #                                                                           #
 #    You should have received a copy of the GNU General Public License      #
 #    along with this program; if not, write to the Free Software            #
 #    Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.              #
 #                                                                           #
 #############################################################################

import pymzml

import numpy

from pyms.GCMS.Class import GCMS_data
from pyms.GCMS.Class import Scan
from pyms.Utils.Utils import is_str
from pyms.Utils.Error import error, stop

try:
    from mpi4py import MPI
except:
    pass

def mzML_reader(file_name):

    """
    @summary: A reader for mzML files, returns
        a GC-MS data object

    @param file_name: The name of the mzML file
    @type file_name: StringType

    @author: Sean O'Callaghan
    """

    if not is_str(file_name):
        error("'file_name' must be a string")
    try:
        mzml_file = pymzml.run.Reader(file_name)
    except:
        error("Cannot open file '%s'" % file_name)

    try:# avoid printing from each rank
        comm = MPI.COMM_WORLD
        rank = comm.Get_rank()
        size = comm.Get_size()
    
        if rank ==0:
            file_names = []
        
            for i in range(1,size):
                recv_buffer = ""
                file_n = comm.recv(recv_buffer, i)
                file_names.append(file_n)

            print " -> Reading mzML files:"
            print file_name
            for file_n in file_names:
                print file_n
        else:
            comm.send(file_name, dest=0)
    except:
        print " -> Reading mzML file '%s'" % (file_name)

    scan_list = []
    time_list = []

    for spectrum in mzml_file:
        mass_list = []
        intensity_list = []

        for mz,i in spectrum.peaks:
            mass_list.append(mz)
            intensity_list.append(i)

        #scan_list.append(Scan(mass_list, intensity_list))
        for element in spectrum.xmlTree:
            # For some reason there are spectra with no time value, 
            # Ignore these????????????
            if element.get('accession') == "MS:1000016": #time value
                # We need time in seconds not minutes
                time_list.append(60*float(element.get('value')))
                scan_list.append(Scan(mass_list, intensity_list))

    #print "time:", len(time_list)
    #print "scan:", len(scan_list)

    data = GCMS_data(time_list, scan_list)

    return data

