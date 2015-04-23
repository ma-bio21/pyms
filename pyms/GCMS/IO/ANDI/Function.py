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

import math, copy

import numpy

from pyms.GCMS.Class import GCMS_data
from pyms.GCMS.Class import Scan
from pyms.Utils.Error import error, stop
from pyms.Utils.Utils import is_str, is_int, is_float, is_number, is_list
from pyms.Utils.Time import time_str_secs
from pycdf import *

try:
    from mpi4py import MPI
except:
    pass

def ANDI_reader(file_name):

    """
    @summary: A reader for ANDI-MS NetCDF files, returns
        a GC-MS data object

    @param file_name: The name of the ANDI-MS file
    @type file_name: StringType

    @author: Qiao Wang
    @author: Andrew Isaac
    @author: Vladimir Likic
    """

    ## TODO: use 'point_count' and allow for zero len scans

    # the keys used to retrieve certain data from the NetCDF file
    __MASS_STRING = "mass_values"
    __INTENSITY_STRING = "intensity_values"
    __TIME_STRING = "scan_acquisition_time"

    if not is_str(file_name):
        error("'file_name' must be a string")
    try:
        file = CDF(file_name)
    except CDFError:
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

            print " -> Reading netCDF files:"
            print file_name
            for file_n in file_names:
                print file_n
        else:
            comm.send(file_name, dest=0)
    except:
        print " -> Reading netCDF file '%s'" % (file_name)





    scan_list = []
    mass = file.var(__MASS_STRING)
    intensity = file.var(__INTENSITY_STRING)
    mass_values = mass.get().tolist()
    mass_list = []
    mass_previous = mass_values[0]
    mass_list.append(mass_previous)
    intensity_values = intensity.get().tolist()
    intensity_list = []
    intensity_previous = intensity_values[0]
    intensity_list.append(intensity_previous)
    if not len(mass_values) == len(intensity_values):
        error("length of mass_list is not equal to length of intensity_list !")
    for i in range(len(mass_values) - 1):
        # assume masses in ascending order until new scan
        if mass_previous <= mass_values[i + 1]:
            #print mass_values[i+1]
            mass_list.append(mass_values[i + 1])
            mass_previous = mass_values[i + 1]
            intensity_list.append(intensity_values[i + 1])
            intensity_previous = intensity_values[i + 1]
        # new scan
        else:
            scan_list.append(Scan(mass_list, intensity_list))
            #print "Added scan"
            mass_previous = mass_values[i + 1]
            intensity_previous = intensity_values[i + 1]
            mass_list = []
            intensity_list = []
            mass_list.append(mass_previous)
            intensity_list.append(intensity_previous)
    # store final scan
    scan_list.append(Scan(mass_list, intensity_list))
    time = file.var(__TIME_STRING)
    time_list = time.get().tolist()

    # sanity check
    if not len(time_list) == len(scan_list):
        error("number of time points (%d) does not equal the number of scans (%d)"%(len(time_list), len(scan_list)))

    data = GCMS_data(time_list, scan_list)

    return data

def ANDI_writer(file_name, im):

    """
    @summary: A reader for ANDI-MS NetCDF files, returns
        a GC-MS data object

    @param file_name: The name of the ANDI-MS file
    @type file_name: StringType
    @param im: The IntensityMatrix
    @type file_name: pyms.GCMS.Class.IntensityMatrix

    @author: Andrew Isaac
    """

    # netCDF header info for compatability
    # attributes
  #dataset_completeness   0 CHAR     6 C1+C2
  #dataset_origin         4 CHAR    16 Santa Clara, CA
  #experiment_date_time_stamp   7 CHAR    20 20081218044500+1100
  #experiment_title       6 CHAR     7 mix ma
  #experiment_type       10 CHAR    25 Centroided Mass Spectrum
  #external_file_ref_0    9 CHAR     8 MA_5C.M
  #languages              3 CHAR     8 English
  #ms_template_revision   1 CHAR     6 1.0.1
  #netcdf_file_date_time_stamp   5 CHAR    20 20090114001531+1100
  #netcdf_revision        2 CHAR     6 2.3.2
  #number_of_times_calibrated  12 INT      1 0
  #number_of_times_processed  11 INT      1 1
  #operator_name          8 CHAR    12 Dave and Su
  #raw_data_intensity_format  25 CHAR     6 Float
  #raw_data_mass_format  23 CHAR     6 Float
  #raw_data_time_format  24 CHAR     6 Short
  #sample_state          13 CHAR    12 Other State
  #test_detector_type    18 CHAR    20 Electron Multiplier
  #test_ionization_mode  16 CHAR    16 Electron Impact
  #test_ionization_polarity  17 CHAR    18 Positive Polarity
  #test_ms_inlet         15 CHAR    17 Capillary Direct
  #test_resolution_type  19 CHAR    20 Constant Resolution
  #test_scan_direction   21 CHAR     3 Up
  #test_scan_function    20 CHAR    10 Mass Scan
  #test_scan_law         22 CHAR     7 Linear
  #test_separation_type  14 CHAR    18 No Chromatography

    # dimensions
  #_128_byte_string       6    128
  #_16_byte_string        3     16
  #_255_byte_string       7    255
  #_2_byte_string         0      2
  #_32_byte_string        4     32
  #_4_byte_string         1      4
  #_64_byte_string        5     64
  #_8_byte_string         2      8
  #error_number          10      1
  #instrument_number     12      1
  #point_number           9 554826   X
  #range                  8      2
  #scan_number           11   9865

    # variables
  #a_d_coaddition_factor   2 SHORT      0 scan_number(9865)
  #a_d_sampling_rate      1 DOUBLE     0 scan_number(9865)
  #actual_scan_number     7 INT        0 scan_number(9865)
  #error_log              0 CHAR       0 error_number(1), _64_byte_string(64)
  #flag_count            15 INT        0 scan_number(9865)
  #instrument_app_version  27 CHAR       0 instrument_number(1),
#_32_byte_string(32)
  #instrument_comments   28 CHAR       0 instrument_number(1),
#_32_byte_string(32)
  #instrument_fw_version  25 CHAR       0 instrument_number(1),
#_32_byte_string(32)
  #instrument_id         20 CHAR       0 instrument_number(1),
#_32_byte_string(32)
  #instrument_mfr        21 CHAR       0 instrument_number(1),
#_32_byte_string(32)
  #instrument_model      22 CHAR       0 instrument_number(1),
#_32_byte_string(32)
  #instrument_name       19 CHAR       0 instrument_number(1),
#_32_byte_string(32)
  #instrument_os_version  26 CHAR       0 instrument_number(1),
#_32_byte_string(32)
  #instrument_serial_no  23 CHAR       0 instrument_number(1),
#_32_byte_string(32)
  #instrument_sw_version  24 CHAR       0 instrument_number(1),
#_32_byte_string(32)
  #intensity_values      18 FLOAT      3 point_number(554826)
  #inter_scan_time        5 DOUBLE     0 scan_number(9865)
  #mass_range_max        10 DOUBLE     0 scan_number(9865)
  #mass_range_min         9 DOUBLE     0 scan_number(9865)
  #mass_values           16 FLOAT      2 point_number(554826)
  #point_count           14 INT        0 scan_number(9865)
  #resolution             6 DOUBLE     0 scan_number(9865)
  #scan_acquisition_time   3 DOUBLE     0 scan_number(9865)
  #scan_duration          4 DOUBLE     0 scan_number(9865)
  #scan_index            13 INT        0 scan_number(9865)
  #time_range_max        12 DOUBLE     0 scan_number(9865)
  #time_range_min        11 DOUBLE     0 scan_number(9865)
  #time_values           17 FLOAT      2 point_number(554826)
  #total_intensity        8 DOUBLE     1 scan_number(9865)

    # variable information
#intensity_values attributes

  #name                 idx type   len value
  #-------------------- --- ----   --- -----
  #add_offset             1 DOUBLE   1 0.0
  #scale_factor           2 DOUBLE   1 1.0
  #units                  0 CHAR    26 Arbitrary Intensity Units

#mass_values attributes

  #name                 idx type   len value
  #-------------------- --- ----   --- -----
  #scale_factor           1 DOUBLE   1 1.0
  #units                  0 CHAR     4 M/Z

#time_values attributes

  #name                 idx type   len value
  #-------------------- --- ----   --- -----
  #scale_factor           1 DOUBLE   1 1.0
  #units                  0 CHAR     8 Seconds

#total_intensity attributes

  #name                 idx type   len value
  #-------------------- --- ----   --- -----
  #units                  0 CHAR    26 Arbitrary Intensity Units

    # netCDF dimension names
    __POINT_NUMBER = "point_number"
    __SCAN_NUMBER = "scan_number"

    # the keys used to create certain data from the NetCDF file
    __MASS_STRING = "mass_values"
    __INTENSITY_STRING = "intensity_values"
    __TIME_STRING = "scan_acquisition_time"
    __POINT_COUNT = "point_count"

    if not is_str(file_name):
        error("'file_name' must be a string")
    try:
        # Open netCDF file in overwrite mode, creating it if inexistent.
        nc = CDF(file_name, NC.WRITE|NC.TRUNC|NC.CREATE)
        # Automatically set define and data modes.
        nc.automode()
    except CDFError:
        error("Cannot create file '%s'" % file_name)

    mass_list = im.get_mass_list()
    time_list = im.get_time_list()

    # direct access, don't modify
    intensity_matrix = im.intensity_matrix

    # compress by ignoring zero intensities
    # included for consistency with imported netCDF format
    mass_values = []
    intensity_values = []
    point_count_values = []
    for row in xrange(len(intensity_matrix)):
        pc = 0  # point count
        for col in xrange(len(intensity_matrix[0])):  # all rows same len
            if (intensity_matrix[row][col] > 0):
                mass_values.append(mass_list[col])
                intensity_values.append(intensity_matrix[row][col])
                pc += 1
        point_count_values.append(pc)

    # sanity checks
    if not len(time_list) == len(point_count_values):
        error("number of time points does not equal the number of scans")

    # create dimensions
    # total number of data points
    dim_point_number = nc.def_dim(__POINT_NUMBER, len(mass_values))
    # number of scans
    dim_scan_number = nc.def_dim(__SCAN_NUMBER, len(point_count_values))

    # create variables
    # points
    var_mass_values = nc.def_var(__MASS_STRING, NC.FLOAT, dim_point_number)
    var_intensity_values = nc.def_var(__INTENSITY_STRING, NC.FLOAT,
        dim_point_number)
    # scans
    var_time_list = nc.def_var(__TIME_STRING, NC.DOUBLE, dim_scan_number)
    var_point_count_values = nc.def_var(__POINT_COUNT, NC.INT,
        dim_scan_number)

    # populate variables
    # points
    var_mass_values[:] = mass_values
    var_intensity_values[:] = intensity_values
    # scans
    var_time_list[:] = time_list
    var_point_count_values[:] = point_count_values

    # close file
    nc.close()

