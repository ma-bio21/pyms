"""
Classes to model GC-MS data
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

import numpy
import math
import copy
import operator

from pyms.Utils.Error import error
from pyms.Utils.Utils import is_str, is_int, is_array, is_list, is_number
from pyms.Utils.IO import open_for_writing, close_for_writing, save_data
from pyms.Utils.Math import mean, std, median
from pyms.Utils.Time import time_str_secs

class GCMS_data(object):

    """
    @summary: Generic object for GC-MS data. Contains raw data
        as a list of scans and times

    @author: Qiao Wang
    @author: Andrew Isaac
    @author: Vladimir Likic
    """

    def __init__(self, time_list, scan_list):

        """
        @summary: Initialize the GC-MS data

        @param time_list: List of scan retention times
        @type time_list: ListType
        @param scan_list: List of Scan objects
        @type scan_list: ListType

        @author: Qiao Wang
        @author: Andrew Isaac
        @author: Vladimir Likic
        """

        if not is_list(time_list) or not is_number(time_list[0]):
            error("'time_list' must be a list of numbers")

        if not is_list(scan_list) or not isinstance(scan_list[0], Scan):
            error("'scan_list' must be a list of Scan objects")

        self.__set_time(time_list)
        self.__scan_list = scan_list
        self.__set_min_max_mass()
        self.__calc_tic()

    def __len__(self):

        """
        @summary: Returns the length of the data object,
            defined as the number of scans

        @return: Number of scans
        @rtype: IntType

        @author: Vladimir Likic
        """

        return len(self.__scan_list)

    def __set_time(self, time_list):

        """
        @summary: Sets time-related properties of the data

        @param time_list: List of retention times
        @type time_list: ListType

        @author: Vladimir Likic
        """

        # calculate the time step, its spreak, and along the way
        # check that retention times are increasing
        time_diff_list = []

        for ii in range(len(time_list)-1):
            t1 = time_list[ii]
            t2 = time_list[ii+1]
            if not t2 > t1:
                error("problem with retention times detected")
            time_diff = t2 - t1
            time_diff_list.append(time_diff)

        time_step = mean(time_diff_list)
        time_step_std = std(time_diff_list)

        self.__time_list = time_list
        self.__time_step = time_step
        self.__time_step_std = time_step_std
        self.__min_rt = min(time_list)
        self.__max_rt = max(time_list)

    def __set_min_max_mass(self):

        """
        @summary: Sets the min and max mass value

        @author: Qiao Wang
        @author: Andrew Isaac
        @author: Vladimir Likic
        """

        mini = self.__scan_list[0].get_min_mass()
        maxi = self.__scan_list[0].get_max_mass()
        for scan in self.__scan_list:
            tmp_mini = scan.get_min_mass()
            tmp_maxi = scan.get_max_mass()
            if tmp_mini < mini:
                mini = tmp_mini
            if tmp_maxi > maxi:
                maxi = tmp_maxi
        self.__min_mass = mini
        self.__max_mass = maxi

    def get_min_mass(self):

        """
        @summary: Get the min mass value over all scans

        @return: The minimum mass of all the data
        @rtype: FloatType

        @author: Qiao Wang
        @author: Andrew Isaac
        @author: Vladimir Likic
        """

        return self.__min_mass

    def get_max_mass(self):

        """
        @summary: Get the max mass value over all scans

        @return: The maximum mass of all the data
        @rtype: FloatType

        @author: Qiao Wang
        @author: Andrew Isaac
        @author: Vladimir Likic
        """

        return self.__max_mass

    def get_index_at_time(self, time):

        """
        @summary: Returns the nearest index corresponding to the given
            time

        @param time: Time in seconds
        @type time: FloatType

        @return: Nearest index corresponding to given time
        @rtype: IntType

        @author: Lewis Lee
        @author: Tim Erwin
        @author: Vladimir Likic
        """

        if not is_number(time):
            error("'time' must be a number")

        if (time < self.__min_rt) or (time > self.__max_rt):
            error("time %.2f is out of bounds (min: %.2f, max: %.2f)" %
                  (time, self.__min_rt, self.__max_rt))

        time_list = self.__time_list
        time_diff_min = self.__max_rt
        ix_match = None

        for ix in range(len(time_list)):

            time_diff = math.fabs(time-time_list[ix])

            if time_diff < time_diff_min:
                ix_match = ix
                time_diff_min = time_diff

        return ix_match

    def get_time_list(self):

        """
        @summary: Returns the list of each scan retention time

        @return: A list of each scan retention time
        @rtype: ListType

        @author: Qiao Wang
        @author: Andrew Isaac
        @author: Vladimir Likic
        """

        return copy.deepcopy(self.__time_list)

    def get_scan_list(self):

        """
        @summary: Return a list of the scan objects

        @return: A list of scan objects
        @rtype: ListType

        @author: Qiao Wang
        @author: Andrew Isaac
        @author: Vladimir Likic
        """

        return copy.deepcopy(self.__scan_list)

    def __get_scan_list(self):
        return self.__scan_list

    scan_list = property(__get_scan_list)

    def get_tic(self):

        """
        @summary: Returns the total ion chromatogram

        @return: Total ion chromatogram
        @rtype: pyms.GCMS.Class.IonChromatogram

        @author: Andrew Isaac
        """

        return self.__tic

    def __calc_tic(self):
        """
        @summary: Calculate the total ion chromatogram

        @return: Total ion chromatogram
        @rtype: pyms.GCMS.Class.IonChromatogram

        @author: Qiao Wang
        @author: Andrew Isaac
        @author: Vladimir Likic
        """

        intensities = []
        for scan in self.__scan_list:
            intensities.append(sum(scan.get_intensity_list()))
        ia = numpy.array(intensities)
        rt = copy.deepcopy(self.__time_list)
        tic = IonChromatogram(ia, rt)

        self.__tic = tic

    def trim(self, begin=None, end=None):

        """
        @summary: trims data in the time domain

        @param begin: begin parameter designating start time or
            scan number
        @type begin: IntType or StrType
        @param end: end parameter designating start time or
            scan number
        @type end: IntType or StrType

            The arguments 'begin' and 'end' can be either integers
            (in which case they are taken as the first/last scan
            number for trimming) or strings in which case they are
            treated as time strings and converted to scan numbers.

            At present both 'begin' and 'end' must be of the same
            type, either both scan numbers or time strings.

        @author: Vladimir Likic
        """

        # trim called with defaults, or silly arguments
        if begin == None and end == None:
            print "Nothing to do."
            return # exit immediately

        N = len(self.__scan_list)

        # process 'begin' and 'end'
        if begin == None:
            first_scan = 0
        elif is_int(begin):
            first_scan = begin-1
        elif is_str(begin):
            time = time_str_secs(begin)
            first_scan = self.get_index_at_time(time) + 1
        else:
            error("invalid 'begin' argument")

        if end == None:
            last_scan = N-1
        elif is_int(end):
            last_scan = end
        elif is_str(end):
            time = time_str_secs(end)
            last_scan = self.get_index_at_time(time) + 1
        else:
            error("invalid 'end' argument")

        # sanity checks
        if not last_scan > first_scan:
            error("last scan=%d, first scan=%d" % (last_scan, first_scan))
        elif first_scan < 0:
            error("scan number must be greater than one")
        elif last_scan > N-1:
            error("last scan=%d, total number of scans=%d" % (last_scan, N))

        print "Trimming data to between %d and %d scans" % \
                (first_scan+1, last_scan+1)

        scan_list_new = []
        time_list_new = []
        for ii in range(len(self.__scan_list)):
            if ii >= first_scan and ii <= last_scan:
                scan = self.__scan_list[ii]
                time = self.__time_list[ii]
                scan_list_new.append(scan)
                time_list_new.append(time)


        # update info
        self.__scan_list = scan_list_new
        self.__set_time(time_list_new)
        self.__set_min_max_mass()
        self.__calc_tic()

    def info(self, print_scan_n=False):

        """
        @summary: Prints some information about the data

        @param print_scan_n: If set to True will print the number
            of m/z values in each scan
        @type print_scan_n: BooleanType

        @author: Vladimir Likic
        """

        # print the summary of simply attributes
        print " Data retention time range: %.3f min -- %.3f min" % \
                (self.__min_rt/60.0, self.__max_rt/60)
        print " Time step: %.3f s (std=%.3f s)" % ( self.__time_step, \
                self.__time_step_std )
        print " Number of scans: %d" % ( len(self.__scan_list) )
        print " Minimum m/z measured: %.3f" % ( self.__min_mass )
        print " Maximum m/z measured: %.3f" % ( self.__max_mass )

        # calculate median number of m/z values measured per scan
        n_list = []
        for ii in range(len(self.__scan_list)):
            scan = self.__scan_list[ii]
            n = len(scan)
            n_list.append(n)
            if print_scan_n: print n
        mz_mean = mean(n_list)
        mz_median = median(n_list)
        print " Mean number of m/z values per scan: %d" % ( mz_mean )
        print " Median number of m/z values per scan: %d" % ( mz_median )

    def write(self, file_root):

        """
        @summary: Writes the entire raw data to two files, one
            'file_root'.I.csv (intensities) and 'file_root'.mz.csv
            (m/z values).

            This method writes two CSV files, containing intensities
            and corresponding m/z values. In general these are not
            two-dimensional matrices, because different scans may
            have different number of m/z values recorded.

        @param file_root: The root for the output file names
        @type file_root: StringType

        @author: Vladimir Likic
        """

        if not is_str(file_root):
            error("'file_root' must be a string")

        file_name1 = file_root + ".I.csv"
        file_name2 = file_root + ".mz.csv"

        print " -> Writing intensities to '%s'" % ( file_name1 )
        print " -> Writing m/z values to '%s'" % ( file_name2 )

        fp1 = open_for_writing(file_name1)
        fp2 = open_for_writing(file_name2)

        for ii in range(len(self.__scan_list)):

            scan = self.__scan_list[ii]

            intensity_list = scan.get_intensity_list()
            mass_list = scan.get_mass_list()

            for ii in range(len(intensity_list)):
                v = intensity_list[ii]
                if ii == 0:
                    fp1.write("%.4f" % (v))
                else:
                    fp1.write(",%.4f" % (v))
            fp1.write("\n")

            for ii in range(len(mass_list)):
                v = mass_list[ii]
                if ii == 0:
                    fp2.write("%.4f" % (v))
                else:
                    fp2.write(",%.4f" % (v))
            fp2.write("\n")

        close_for_writing(fp1)
        close_for_writing(fp2)

    def write_intensities_stream(self, file_name):

        """
        @summary: Writes all intensities to a file

        @param file_name: Output file name
        @type file_name: StringType

        This function loop over all scans, and for each scan
        writes intensities to the file, one intenisity per
        line. Intensities from different scans are joined
        without any delimiters.

        @author: Vladimir Likic
        """

        if not is_str(file_name):
            error("'file_name' must be a string")

        N = len(self.__scan_list)

        print" -> Writing scans to a file"

        fp = open_for_writing(file_name)

        for ii in range(len(self.__scan_list)):
            scan = self.__scan_list[ii]
            intensities = scan.get_intensity_list()
            for I in intensities:
                fp.write("%8.4f\n" % ( I ) )

        close_for_writing(fp)

class Scan(object):

    """
    @summary: Generic object for a single Scan's raw data

    @author: Qiao Wang
    @author: Andrew Isaac
    @author: Vladimir Likic
    """

    def __init__(self, mass_list, intensity_list):

        """
        @summary: Initialize the Scan data

        @param mass_list: mass values
        @type mass_list: ListType

        @param intensity_list: intensity values
        @type intensity_list: ListType

        @author: Qiao Wang
        @author: Andrew Isaac
        @author: Vladimir Likic
        """
        #print "mass_list[0]",mass_list[0]
        if not is_list(mass_list) or not is_number(mass_list[0]):
            error("'mass_list' must be a list of numbers")
        if not is_list(intensity_list) or \
           not is_number(intensity_list[0]):
            error("'intensity_list' must be a list of numbers")

        self.__mass_list = mass_list
        self.__intensity_list = intensity_list
        self.__min_mass = min(mass_list)
        self.__max_mass = max(mass_list)

    def __len__(self):

        """
        @summary: Returns the length of the Scan object

        @return: Length of Scan
        @rtype: IntType

        @author: Andrew Isaac
        """

        return len(self.__mass_list)

    def get_mass_list(self):

        """
        @summary: Returns the masses for the current scan

        @return: the masses
        @rtype: ListType

        @author: Qiao Wang
        @author: Andrew Isaac
        @author: Vladimir Likic
        """
 
        return copy.deepcopy(self.__mass_list)

    def __get_mass_list(self):
        return self.__mass_list

    mass_list = property(__get_mass_list)


    def get_intensity_list(self):

        """
        @summary: Returns the intensities for the current scan

        @return: the intensities
        @rtype: ListType

        @author: Qiao Wang
        @author: Andrew Isaac
        @author: Vladimir Likic
        """
       
        return copy.deepcopy(self.__intensity_list)

    def __get_intensity_list(self):
        return self.__intensity_list

    intensity_list = property(__get_intensity_list)
  

    def get_min_mass(self):

        """
        @summary: Returns the minimum m/z value in the scan

        @return: Minimum m/z
        @rtype: Float

        @author: Andrew Isaac
        """

        return self.__min_mass

    def get_max_mass(self):

        """
        @summary: Returns the maximum m/z value in the scan

        @return: Maximum m/z
        @rtype: Float

        @author: Andrew Isaac
        """

        return self.__max_mass

class IntensityMatrix(object):

    """
    @summary: Intensity matrix of binned raw data

    @author: Andrew Isaac
    """

    def __init__(self, time_list, mass_list, intensity_matrix):

        """
        @summary: Initialize the IntensityMatrix data

        @param time_list: Retention time values
        @type time_list: ListType

        @param mass_list: Binned mass values
        @type mass_list: ListType

        @param intensity_matrix: Binned intensity values per scan
        @type intensity_matrix: ListType

        @author: Andrew Isaac
        """

        # sanity check
        if not is_list(time_list) or not is_number(time_list[0]):
            error("'time_list' must be a list of numbers")
        if not is_list(mass_list) or not is_number(mass_list[0]):
            error("'mass_list' must be a list of numbers")
        if not is_list(intensity_matrix) or \
           not is_list(intensity_matrix[0]) or \
           not is_number(intensity_matrix[0][0]):
            error("'intensity_matrix' must be a list, of a list, of numbers")
        if not len(time_list) == len(intensity_matrix):
            error("'time_list' is not the same length as 'intensity_matrix'")
        if not len(mass_list) == len(intensity_matrix[0]):
            error("'mass_list' is not the same size as 'intensity_matrix'"
                " width")

        self.__time_list = time_list
        self.__mass_list = mass_list
        self.__intensity_matrix = intensity_matrix

        self.__min_mass = min(mass_list)
        self.__max_mass = max(mass_list)

        # Direct access for speed (DANGEROUS)
        self.intensity_matrix = self.__intensity_matrix

        # Try to include parallelism.
        try:
            from mpi4py import MPI
            comm = MPI.COMM_WORLD
            num_ranks = comm.Get_size()
            rank = comm.Get_rank()
            M, N = len(intensity_matrix), len(intensity_matrix[0])
            lrr = (rank*M/num_ranks, (rank + 1)*M/num_ranks)
            lcr = (rank*N/num_ranks, (rank + 1)*N/num_ranks)
            m, n = (lrr[1] - lrr[0], lcr[1] - lcr[0])
            self.comm = comm
            self.num_ranks = num_ranks
            self.rank = rank
            self.M = M
            self.N = N
            self.local_row_range = lrr
            self.local_col_range = lcr
            self.m = m
            self.n = n

        # If we can't import mpi4py then continue in serial.
        except:
            pass

    def get_local_size(self):
        """
        @summary: Gets the local size of intensity matrix.

        @return: Number of rows and cols
        @rtype: IntType

        @author: Luke Hodkinson
        """

        # Check for parallel.
        if hasattr(self, 'comm'):
            return self.m, self.n

        # If serial call the regular routine.
        return self.get_size()

    def get_size(self):

        """
        @summary: Gets the size of intensity matrix

        @return: Number of rows and cols
        @rtype: IntType

        @author: Qiao Wang
        @author: Andrew Isaac
        @author: Luke Hodkinson
        @author: Vladimir Likic
        """

        n_scan = len(self.__intensity_matrix)
        n_mz = len(self.__intensity_matrix[0])

        return n_scan, n_mz

    def iter_ms_indices(self):
        """
        @summary: Iterates over local row indices

        @return: Current row index
        @rtype: IntType

        @author: Luke Hodkinson
        """

        # Check for parallel.
        if hasattr(self, 'comm'):
            # At the moment we assume we break the matrix into contiguous
            # ranges. We've allowed for this to change by wrapping up the
            # iteration in this method.
            for i in xrange(self.local_row_range[0], self.local_row_range[1]):
                yield i

        else:
            # Iterate over global indices.
            n_scan = len(self.__intensity_matrix)
            for i in xrange(0, n_scan):
                yield i

    def iter_ic_indices(self):
        """
        @summary: Iterate over local column indices

        @return: Current column index
        @rtype: IntType

        @author: Luke Hodkinson
        """

        # Check for parallel.
        if hasattr(self, 'comm'):
            # At the moment we assume we break the matrix into contiguous
            # ranges. We've allowed for this to change by wrapping up the
            # iteration in this method.
            for i in xrange(self.local_col_range[0], self.local_col_range[1]):
                yield i

        else:
            # Iterate over global indices.
            n_mz = len(self.__intensity_matrix[0])
            for i in xrange(0, n_mz):
                yield i

    def set_ic_at_index(self, ix, ic):

        """
        @summary: Sets the ion chromatogram specified by index to a new
            value

        @param ix: Index of an ion chromatogram in the intensity data
            matrix to be set
        @type ix: IntType
        @param ic: Ion chromatogram that will be copied at position 'ix'
            in the data matrix
        @type: IonChromatogram

        The length of the ion chromatogram must match the appropriate
        dimension of the intensity matrix.

        @author: Vladimir Likic
        """

        if not is_int(ix):
            error("index not an integer")

        # this returns an numpy.array object
        ia = ic.get_intensity_array()

        # check if the dimension is ok
        if len(ia) != len(self.__intensity_matrix):
            error("ion chromatogram incompatible with the intensity matrix")
        else:
           N = len(ia)

        # Convert 'ia' to a list. By convention, the attribute
        # __intensity_matrix of the class IntensityMatrix is a list
        # of lists. This makes pickling instances of IntensityMatrix
        # practically possible, since pickling numpy.array objects
        # produces ten times larger files compared to pickling python
        #lists.
        ial = ia.tolist()
        for i in range(N):
            self.__intensity_matrix[i][ix] = ial[i]

    def get_ic_at_index(self, ix):

        """
        @summary: Returns the ion chromatogram at the specified index

        @param ix: Index of an ion chromatogram in the intensity data
            matrix
        @type ix: IntType

        @return: Ion chromatogram at given index
        @rtype: IonChromatogram

        @author: Qiao Wang
        @author: Andrew Isaac
        @author: Vladimir Likic
        """

        if not is_int(ix):
            error("index not an integer")
        try:
            ia = []
            for i in range(len(self.__intensity_matrix)):
                ia.append(self.__intensity_matrix[i][ix])
        except IndexError:
            error("index out of bounds.")

        ic_ia = numpy.array(ia)
        mass = self.get_mass_at_index(ix)
        rt = copy.deepcopy(self.__time_list)

        return IonChromatogram(ic_ia, rt, mass)

    def get_ic_at_mass(self, mass = None):

        """
        @summary: Returns the ion chromatogram for the specified mass.
            The nearest binned mass to mass is used.

            If no mass value is given, the function returns the total
            ion chromatogram.

        @param mass: Mass value of an ion chromatogram
        @type mass: IntType

        @return: Ion chromatogram for given mass
        @rtype: IonChromatogram

        @author: Andrew Isaac
        @author: Vladimir Likic
        """

        if mass == None:
            return self.get_tic()

        if mass < self.__min_mass or mass > self.__max_mass:
            print "min mass: ", self.__min_mass, "max mass:", self.__max_mass
            error("mass is out of range")

        ix = self.get_index_of_mass(mass)

        return self.get_ic_at_index(ix)

    def get_mass_list(self):

        """
        @summary: Returns a list of the binned masses

        @return: Binned mass list
        @rtype: ListType

        @author: Qiao Wang
        @author: Andrew Isaac
        @author: Vladimir Likic
        """

        return copy.deepcopy(self.__mass_list)

    def get_ms_at_index(self, ix):

        """
        @summary: Returns a mass spectrum for a given scan index

        @param ix: The index of the scan
        @type ix: IntType

        @return: Mass spectrum
        @rtype: pyms.GCMS.Class.MassSpectrum

        @author: Andrew Isaac
        """

        # TODO: should a deepcopy be returned?

        if not is_int(ix):
            error("index not an integer")

        scan = self.get_scan_at_index(ix)

        return MassSpectrum(self.__mass_list, scan)

    def get_scan_at_index(self, ix):

        """
        @summary: Returns the spectral intensities for scan index

        @param ix: The index of the scan
        @type ix: IntType

        @return: Intensity values of scan spectra
        @rtype: ListType

        @author: Andrew Isaac
        """

        if not is_int(ix):
            error("index not an integer")

        if ix < 0 or ix >= len(self.__intensity_matrix):
            error("index out of range")

        return copy.deepcopy(self.__intensity_matrix[ix])

    def get_min_mass(self):

        """
        @summary: Returns the maximum binned mass

        @return: The maximum binned mass
        @rtype: FloatType

        @author: Andrew Isaac
        """

        return self.__min_mass

    def get_max_mass(self):

        """
        @summary: Returns the maximum binned mass

        @return: The maximum binned mass
        @rtype: FloatType

        @author: Andrew Isaac
        """

        return self.__max_mass

    def get_mass_at_index(self, ix):

        """
        @summary: Returns binned mass at index

        @param ix: Index of binned mass
        @type ix: IntType

        @return: Binned mass
        @rtype: IntType

        @author: Andrew Isaac
        """

        if not is_int(ix):
            error("index not an integer")

        if ix < 0 or ix >= len(self.__mass_list):
            error("index out of range")

        return self.__mass_list[ix]

    def get_index_of_mass(self, mass):

        """
        @summary: Returns the index of mass in the list of masses

        The nearest binned mass to given mass is used.

        @param mass: Mass to lookup in list of masses
        @type mass: FloatType

        @return: Index of mass closest to given mass
        @rtype: IntType

        @author: Andrew Isaac
        """

        best = self.__max_mass
        ix = 0
        for ii in range(len(self.__mass_list)):
            tmp = abs(self.__mass_list[ii] - mass)
            if tmp < best:
                best = tmp
                ix = ii
        return ix

    def get_matrix_list(self):

        """
        @summary: Returns a copy of the intensity matrix as a
            list of lists of floats

        @return: Matrix of intensity values
        @rtype: ListType

        @author: Andrew Isaac
        """

        return copy.deepcopy(self.__intensity_matrix)

    def get_time_list(self):

        """
        @summary: Returns a copy of the time list

        @return: List of retention times
        @rtype: ListType

        @author: Andrew Isaac
        """

        return copy.deepcopy(self.__time_list)

    def get_index_at_time(self, time):

        """
        @summary: Returns the nearest index corresponding to the given time

        @param time: Time in seconds
        @type time: FloatType

        @return: Nearest index corresponding to given time
        @rtype: IntType

        @author: Lewis Lee
        @author: Tim Erwin
        @author: Vladimir Likic
        """

        if not is_number(time):
            error("'time' must be a number")

        if time < min(self.__time_list) or time > max(self.__time_list):
            error("time %.2f is out of bounds (min: %.2f, max: %.2f)" %
                  (time, self.__min_rt, self.__max_rt))

        time_list = self.__time_list
        time_diff_min = max(self.__time_list)
        ix_match = None

        for ix in range(len(time_list)):

            time_diff = math.fabs(time-time_list[ix])

            if time_diff < time_diff_min:
                ix_match = ix
                time_diff_min = time_diff

        return ix_match

    def crop_mass(self, mass_min, mass_max):

        """
        @summary: Crops mass spectrum

        @param mass_min: Minimum mass value
        @type mass_min: IntType or FloatType
        @param mass_max: Maximum mass value
        @type mass_max: IntType or FloatType

        @return: none
        @rtype: NoneType

        @author: Andrew Isaac
        """

        if not is_number(mass_min) or not is_number(mass_max):
            error("'mass_min' and 'mass_max' must be numbers")
        if mass_min >= mass_max:
            error("'mass_min' must be less than 'mass_max'")
        if mass_min < self.__min_mass:
            error("'mass_min' is less than the smallest mass: %.3f" %
                self.__min_mass)
        if mass_max > self.__max_mass:
            error("'mass_max' is greater than the largest mass: %.3f" %
                self.__max_mass)

        # pre build mass_list and list of indecies
        mass_list = self.__mass_list
        new_mass_list = []
        ii_list = []
        for ii in range(len(mass_list)):
             mass = mass_list[ii]
             if mass >= mass_min and mass <= mass_max:
                 new_mass_list.append(mass)
                 ii_list.append(ii)

        # update intensity matrix
        im = self.__intensity_matrix
        for spec_jj in range(len(im)):
            new_spec = []
            for ii in ii_list:
                new_spec.append(im[spec_jj][ii])
            im[spec_jj] = new_spec

        self.__mass_list = new_mass_list
        self.__min_mass = min(new_mass_list)
        self.__max_mass = max(new_mass_list)

    def null_mass(self, mass):

        """
        @summary: Ignore given (closest) mass in spectra

        @param mass: Mass value to remove
        @type mass: IntType or FloatType

        @author: Andrew Isaac
        """

        if not is_number(mass):
            error("'mass' must be numbers")
        if mass < self.__min_mass or mass > self.__max_mass:
            error("'mass' not in mass range: %.3f to %.3f" % (self.__min_mass, \
                self.__max_mass))

        ii = self.get_index_of_mass(mass)

        im = self.__intensity_matrix
        for spec_jj in range(len(im)):
            im[spec_jj][ii] = 0

    def reduce_mass_spectra(self, N=5):

        """
        @summary: Reduces mass spectra by retaining top N
        intensities, discarding all other intensities.

        @param N: The number of top intensities to keep
        @type N: IntType

        @author: Vladimir Likic
        """

        # loop over all mass spectral scans
        for ii in range(len(self.__intensity_matrix)):

            # get the next mass spectrum as list of intensities
            intensity_list = self.__intensity_matrix[ii]
            n = len(intensity_list)

            # get the indices of top N intensities
            top_indices = range(n)
            top_indices.sort(key=lambda i: intensity_list[i], reverse=True)
            top_indices = top_indices[:N]

            # initiate new mass spectrum, and retain only top N intensities
            intensity_list_new = []

            for jj in range(n):
                intensity_list_new.append(0.0)
                if jj in top_indices:
                    intensity_list_new[jj] = intensity_list[jj]

            self.__intensity_matrix[ii] = intensity_list_new

    def export_ascii(self, root_name, format='dat'):

        """
        @summary: Exports the intensity matrix, retention time vector, and
        m/z vector to the ascii format

        By default, export_ascii("NAME") will create NAME.im.dat, NAME.rt.dat,
        and NAME.mz.dat where these are the intensity matrix, retention
        time vector, and m/z vector in tab delimited format. If format='csv',
        the files will be in the CSV format, named NAME.im.csv, NAME.rt.csv,
        and NAME.mz.csv.

        @param root_name: Root name for the output files
        @type root_name: StringType

        @return: none
        @rtype: NoneType

        @author: Milica Ng
        @author: Andrew Isaac
        @author: Vladimir Likic
        """

        if not is_str(root_name):
            error("'root_name' is not a string")

        if format == 'dat':
            separator = " "
            extension = ".dat"
        elif format == 'csv':
            separator = ","
            extension = ".csv"
        else:
            error("unkown format '%s'. Only 'dat' or 'csv' supported" % format)

        # export 2D matrix of intensities
        vals = self.__intensity_matrix
        save_data(root_name+'.im'+extension, vals, sep=separator)

        # export 1D vector of m/z's, corresponding to rows of
        # the intensity matrix
        mass_list = self.__mass_list
        save_data(root_name+'.mz'+extension, mass_list, sep=separator)

        # export 1D vector of retention times, corresponding to
        # columns of the intensity matrix
        time_list = self.__time_list
        save_data(root_name+'.rt'+extension, time_list, sep=separator)

    def export_leco_csv(self, file_name):

        """
        @summary: Exports data in LECO CSV format

        @param file_name: File name
        @type file_name: StringType

        @return: none
        @rtype: NoneType

        @author: Andrew Isaac
        @author: Vladimir Likic
        """

        if not is_str(file_name):
            error("'file_name' is not a string")

        mass_list = self.__mass_list
        time_list = self.__time_list
        vals = self.__intensity_matrix

        fp = open_for_writing(file_name)

        # Format is text header with:
        # "Scan","Time",...
        # and the rest is "TIC" or m/z as text, i.e. "50","51"...
        # The following lines are:
        # scan_number,time,value,value,...
        # scan_number is an int, rest seem to be fixed format floats.
        # The format is 0.000000e+000

        # write header
        fp.write("\"Scan\",\"Time\"")
        for ii in mass_list:
            if is_number(ii):
                fp.write(",\"%d\"" % int(ii))
            else:
                error("mass list datum not a number")
        fp.write("\r\n")  # windows CR/LF

        # write lines
        for ii in range(len(time_list)):
            fp.write("%s,%#.6e" % (ii, time_list[ii]))
            for jj in range(len(vals[ii])):
                if is_number(vals[ii][jj]):
                    fp.write(",%#.6e" % (vals[ii][jj]))
                else:
                    error("datum not a number")
            fp.write("\r\n")

        close_for_writing(fp)

    def import_leco_csv(self, file_name):
        """
        @summary: Imports data in LECO CSV format

        @param file_name: File name
        @type file_name: StringType

        @return: Data as an IntensityMatrix
        @rtype: pyms.GCMS.Class.IntensityMatrix

        @author: Andrew Isaac
        """

        if not is_str(file_name):
            error("'file_name' not a string")

        lines_list = open(file_name,'r')
        data = []
        time_list = []
        mass_list = []

        # Format is text header with:
        # "Scan","Time",...
        # and the rest is "TIC" or m/z as text, i.e. "50","51"...
        # The following lines are:
        # scan_number,time,value,value,...
        # scan_number is an int, rest seem to be fixed format floats.
        # The format is 0.000000e+000

        num_mass = 0
        FIRST = True
        HEADER = True
        data_col = -1
        time_col = -1
        # get each line
        for line in lines_list:
            cols = -1
            data_row = []
            if len(line.strip()) > 0:
                data_list = line.strip().split(',')
                # get each value in line
                for item in data_list:
                    item = item.strip()
                    item = item.strip('\'"')  # remove quotes (in header)

                    # Get header
                    if HEADER:
                        cols += 1
                        if len(item) > 0:
                            if item.lower().find("time") > -1:
                                time_col = cols
                            try:
                                value = float(item)
                                # find 1st col with number as header
                                if FIRST and value > 1:  # assume >1 mass
                                    data_col = cols
                                    # assume time col is previous col
                                    if time_col < 0:
                                        time_col = cols -1
                                    FIRST = False
                                mass_list.append(value)
                                num_mass += 1
                            except ValueError:
                                pass
                    # Get rest
                    else:
                        cols += 1
                        if len(item) > 0:
                            try:
                                value = float(item)
                                if cols == time_col:
                                    time_list.append(value)
                                elif cols >= data_col:
                                    data_row.append(value)
                            except ValueError:
                                pass

                # check row length
                if not HEADER:
                    if len(data_row) == num_mass:
                        data.append(data_row)
                    else:
                        print ("Warning: ignoring row")

                HEADER = False

        # check col lengths
        if len(time_list) != len(data):
            print ("Warning: number of data rows and time list length differ")

        self.__mass_list = mass_list
        self.__time_list = time_list
        self.__intensity_matrix = data
        # Direct access for speed (DANGEROUS)
        self.intensity_matrix = self.__intensity_matrix

## get_ms_at_time()

class IonChromatogram(object):

    """
    @summary: Models an ion chromatogram

    An ion chromatogram is a set of intensities as a function of retention
    time. This can can be either m/z channel intensities (for example, ion
    chromatograms at m/z=65), or cumulative intensities over all measured
    m/z. In the latter case the ion chromatogram is total ion chromatogram
    (TIC).

    The nature of an IonChromatogram object can be revealed by inspecting
    the value of the attribute '__mass'. This is se to the m/z value of the
    ion chromatogram, or to None for TIC.

    @author: Lewis Lee
    @author: Vladimir Likic
    """

    def __init__(self, ia, time_list, mass=None):

        """
        @param ia: Ion chromatogram intensity values
        @type ia: numpy.array
        @param time_list: A list of ion chromatogram retention times
        @type time_list: ListType
        @param mass: Mass of ion chromatogram (Null if TIC)
        @type mass: IntType

        @author: Lewis Lee
        @author: Vladimir Likic
        @author: Vladimir Likic
        """

        if not isinstance(ia, numpy.ndarray):
            error("'ia' must be a numpy array")

        if not is_list(time_list) or not is_number(time_list[0]):
            error("'time_list' must be a list of numbers")

        if len(ia) != len(time_list):
            error("Intensity array and time list differ in length")

        self.__ia = ia
        self.__time_list = time_list
        self.__mass = mass
        self.__time_step = self.__calc_time_step(time_list)
        self.__min_rt = min(time_list)
        self.__max_rt = max(time_list)

    def __len__(self):

        """
        @summary: Returns the length of the IonChromatogram object

        @return: Length of ion chromatogram
        @rtype: IntType

        @author: Lewis Lee
        @author: Vladimir Likic
        """

        return self.__ia.size
    
    def __sub__(self, Other):
        """
        @summary: Subtracts another IC from the current one
        
        @param other: Another IC
        @type other: pyms.GCMS.IonChromatogram
        """
        
        ia_for_sub = Other.get_intensity_array()
        
        for i in range(self.__ia.size):
            self.__ia[i] = self.__ia[i] - ia_for_sub[i]
            

    def get_intensity_at_index(self, ix):

        """
        @summary: Returns intensity at given index

        @param ix: An index
        @type ix: IntType

        @return: Intensity value
        @rtype: FloatType

        @author: Lewis Lee
        @author: Vladimir Likic
        """

        if not is_int(ix):
            error("index not an integer")

        if ix < 0 or ix > self.__ia.size - 1:
            error("index out of bounds")

        return self.__ia[ix]

    def get_intensity_array(self):

        """
        @summary: Returns the entire intensity array

        @return: Intensity array
        @rtype: numpy.ndarray

        @author: Lewis Lee
        @author: Vladimir Likic
        """

        return self.__ia

    def get_time_at_index(self, ix):

        """
        @summary: Returns time at given index

        @param ix: An index
        @type ix: IntType

        @return: Time value
        @rtype: FloatType

        @author: Lewis Lee
        @author: Vladimir Likic
        """

        if not is_int(ix):
            error("index not an integer")

        if ix < 0 or ix > len(self.__time_list) - 1:
            error("index out of bounds")

        return self.__time_list[ix]

    def get_time_list(self):

        """
        @summary: Returns the time list

        @return: Time list
        @rtype: ListType

        @author: Lewis Lee
        @author: Vladimir Likic
        """

        return self.__time_list
    
    def get_mass(self):
        
        """
        @summary: Returns the m/z channel of the IC
        
        @return: m/z channel of the IC
        @rtype: intType
        
        @author: Sean O'Callaghan
        """
        if self.__mass == None:
            error("TIC has no m/z label")
        
        return self.__mass
        

    def get_time_step(self):

        """
        @summary: Returns the time step

        @return: Time step
        @rtype: FloatType

        @author: Lewis Lee
        @author: Vladimir Likic
        """

        return self.__time_step

    def __calc_time_step(self, time_list):

        """
        @summary: Calculates the time step

        @param time_list: A list of retention times
        @type time_list: ListType

        @return: Time step value
        @rtype: FloatType

        @author: Lewis Lee
        @author: Vladimir Likic
        """

        td_list = []
        for ii in range(len(time_list)-1):
            td = time_list[ii+1]-time_list[ii]
            td_list.append(td)

        td_array = numpy.array(td_list)
        time_step = td_array.mean()

        return time_step

    def get_index_at_time(self, time):

        """
        @summary: Returns the nearest index corresponding to the given time

        @param time: Time in seconds
        @type time: FloatType

        @return: Nearest index corresponding to given time
        @rtype: IntType

        @author: Lewis Lee
        @author: Tim Erwin
        @author: Milica Ng
        @author: Vladimir Likic
        """

        if not is_number(time):
            error("'time' must be a number")

        if time < self.__min_rt or time > self.__max_rt:
            error("time %.2f is out of bounds (min: %.2f, max: %.2f)" %
                  (time, self.__min_rt, self.__max_rt))

        time_list = self.__time_list
        time_diff_min = self.__max_rt
        ix_match = None

        for ix in range(len(time_list)):


            time_diff = math.fabs(time-time_list[ix])

            if time_diff < time_diff_min:
                ix_match = ix
                time_diff_min = time_diff

        return ix_match

    def is_tic(self):

        """
        @summary: Returns True if the ion chromatogram is a total ion
            chromatogram (TIC), or False otherwise

        @return: A boolean value indicating if the ion chromatogram
            is a total ion chromatogram (True) or not (False)
        @rtype: BooleanType

        @author: Lewis Lee
        @author: Vladimir Likic
        """

        if self.__mass == None:
            return True
        else:
            return False

    def set_intensity_array(self, ia):

        """
        @summary: Sets the value for the intensity array

        @param ia: An array of new intensity values
        @type ia: numpy.ndarray

        @return: none
        @rtype: NoneType

        @author: Vladimir Likic
        """

        self.__ia = ia


    def write(self, file_name, minutes=False):

        """
        @summary: Writes the ion chromatogram to the specified file

        @param file_name: Output file name
        @type file_name: StringType
        @param minutes: A boolean value indicating whether to write
            time in minutes
        @type minutes: BooleanType

        @return: none
        @rtype: NoneType

        @author: Lewis Lee
        @author: Vladimir Likic
        """

        if not is_str(file_name):
            error("'file_name' must be a string")

        fp = open_for_writing(file_name)

        time_list = copy.deepcopy(self.__time_list)

        if minutes:
            for ii in range(len(time_list)):
                time_list[ii] = time_list[ii]/60.0

        for ii in range(len(time_list)):
            fp.write("%8.4f %#.6e\n" % (time_list[ii], self.__ia[ii]))

        close_for_writing(fp)

class MassSpectrum(object):

    """
    @summary: Models a binned mass spectrum

    @author: Andrew Isaac
    @author: Qiao Wang
    @author: Vladimir Likic
    """

    def __init__(self, mass_list, intensity_list):

        """
        @summary: Initialise the MassSpectrum

        @param mass_list: List of binned masses
        @type mass_list: ListType
        @param intensity_list: List of binned intensities
        @type intensity_list: ListType

        @author: Andrew Isaac
        @author: Qiao Wang
        @author: Vladimir Likic
        """

        if not is_list(mass_list) or not is_number(mass_list[0]):
            error("'mass_list' must be a list of numbers")
        if not is_list(intensity_list) or \
           not is_number(intensity_list[0]):
            error("'intensity_list' must be a list of numbers")
        if not len(mass_list) == len(intensity_list):
            error("'mass_list' is not the same size as 'intensity_list'")

        #TODO: should these be public, or accessed through methods???
        self.mass_list = mass_list
        self.mass_spec = intensity_list

    def __len__(self):

        """
        @summary: Length of the MassSpectrum

        @return: Length of the MassSpectrum (Number of bins)
        @rtype: IntType

        @author: Andrew Isaac
        @author: Qiao Wang
        @author: Vladimir Likic
        """

        return len(self.mass_list)

