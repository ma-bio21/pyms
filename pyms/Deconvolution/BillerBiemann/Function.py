"""
BillerBiemann deconvolution algorithm
Provides a class to perform Biller and Biemann deconvolution
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
import copy

from pyms.Utils.Error import error
from pyms.Utils.Utils import is_list, is_number, is_int
from pyms.GCMS.Class import IonChromatogram, MassSpectrum
from pyms.Peak.Class import Peak

# If psyco is installed, use it to speed up running time
try:
    import psyco
    psyco.full()
except:
    pass

#######################
# structure
# 1) find local maxima per ion, store intensity and scan index
# 2) sum across N scans to compensate for scan type
# 3) sum ions belonging to each maxima scan
#######################

def BillerBiemann(im, points=3, scans=1):

    """
    @summary: BillerBiemann Deconvolution

        Deconvolution based on the algorithm of Biller and Biemann (1974)

    @param im: An IntensityMatrix object
    @type im: pyms.GCMS.Class.IntensityMatrix
    @param points: Peak if maxima over 'points' number of scans (Default 3)
    @type points: IntType
    @param scans: To compensate for spectra skewing,
        peaks from 'scans' scans are combined (Default 1).
    @type scans: IntType

    @return: List of Peak objects
    @rtype: ListType

    @author: Andrew Isaac
    """

    rt_list = im.get_time_list()
    mass_list = im.get_mass_list()
    peak_list = []
    maxima_im = get_maxima_matrix(im, points, scans)
    numrows = len(maxima_im)
    for row in range(numrows):
        if sum(maxima_im[row]) > 0:
            rt = rt_list[row]
            ms = MassSpectrum(mass_list, maxima_im[row])
            peak = Peak(rt, ms)
            peak.set_pt_bounds([0,row,0])  # store IM index for convenience
            peak_list.append(peak)

    return peak_list

def rel_threshold(pl, percent=2):

    """
    @summary: Remove ions with relative intensities less than the given
        relative percentage of the maximum intensity.

    @param pl: A list of Peak objects
    @type pl: ListType
    @param percent: Threshold for relative percentage of intensity (Default 2%)
    @type percent: FloatType

    @return: A new list of Peak objects with threshold ions
    @rtype: ListType

    @author: Andrew Isaac
    """

    if not is_number(percent) or percent <= 0:
        error("'percent' must be a number > 0")

    pl_copy = copy.deepcopy(pl)
    new_pl = []
    for p in pl_copy:
        ms = p.get_mass_spectrum()
        ia = ms.mass_spec
        # assume max(ia) big so /100 1st
        cutoff = (max(ia)/100.0)*float(percent)
        for i in range(len(ia)):
            if ia[i] < cutoff:
                ia[i] = 0
        ms.mass_spec = ia
        p.set_mass_spectrum(ms)
        new_pl.append(p)
    return new_pl

def num_ions_threshold(pl, n, cutoff):

    """
    @summary: Remove Peaks where there are less than a given number of ion
        intensities above the given threshold

    @param pl: A list of Peak objects
    @type pl: ListType
    @param n: Minimum number of ions that must have intensities above the cutoff
    @type n: IntType
    @param cutoff: The minimum intensity threshold
    @type cutoff: FloatType

    @return: A new list of Peak objects
    @rtype: ListType

    @author: Andrew Isaac
    """

    pl_copy = copy.deepcopy(pl)
    new_pl = []
    for p in pl_copy:
        ms = p.get_mass_spectrum()
        ia = ms.mass_spec
        ions = 0
        for i in range(len(ia)):
            if ia[i] >= cutoff:
                ions += 1
        if ions >= n:
            new_pl.append(p)
    return new_pl

def sum_maxima(im, points=3, scans=1):

    """
    @summary: Reconstruct the TIC as sum of maxima

    @param im: An IntensityMatrix object
    @type im: pyms.GCMS.Class.IntensityMatrix
    @param points: Peak if maxima over 'points' number of scans
    @type points: IntType
    @param scans: To compensate for spectra scewing,
        peaks from 'scans' scans are combined.
    @type scans: IntType

    @return: The reconstructed TIC
    @rtype: pyms.GCMS.Class.IonChromatogram

    @author: Andrew Isaac
    """

    maxima_im = get_maxima_matrix(im, points)
    sums = []
    numrows = len(maxima_im)
    half = int(scans/2)
    for row in range(numrows):
        val = 0
        for ii in range(scans):
            if row - half + ii >= 0 and row - half + ii < numrows:
                val += maxima_im[row - half + ii].sum()
        sums.append(val)
    tic = IonChromatogram(numpy.array(sums), im.get_time_list())

    return tic

def get_maxima_indices(ion_intensities, points=3):

    """
    @summary: Find local maxima.

    @param ion_intensities: A list of intensities for a single ion
    @type ion_intensities: ListType
    @param points: Peak if maxima over 'points' number of scans
    @type points: IntType

    @return: A list of scan indices
    @rtype: ListType

    @author: Andrew Isaac
    """

    if not is_list(ion_intensities) or not is_number(ion_intensities[0]):
        error("'ion_intensities' must be a List of numbers")

    # find peak inflection points
    # use a 'points' point window
    # for a plateau after a rise, need to check if it is the left edge of
    # a peak
    peak_point = []
    edge = -1
    points = int(points)
    half = int(points/2)
    points = 2*half+1  # ensure odd number of points
    for index in range(len(ion_intensities)-points+1):
        left = ion_intensities[index:index+half]
        mid = ion_intensities[index+half]
        right = ion_intensities[index+half+1:index+points]
        # max in middle
        if mid > max(left) and mid > max(right):
            peak_point.append(index+half)
            edge = -1  # ignore previous rising edge
        # flat from rise (left of peak?)
        if mid > max(left) and mid == max(right):
            edge = index+half  # ignore previous rising edge, update latest
        # fall from flat
        if mid == max(left) and mid > max(right):
            if edge > -1:
                centre = int((edge+index+half)/2)  # mid point
                peak_point.append(centre)
            edge = -1

    return peak_point

def get_maxima_list(ic, points=3):

    """
    @summary: List of retention time and intensity of local maxima for ion

    @param ic: An IonChromatogram object
    @type ic: pyms.GCMS.Class.IonChromatogram
    @param points: Peak if maxima over 'points' number of scans
    @type points: IntType

    @return: A list of retention time and intensity of local maxima for ion
    @rtype: ListType

    @author: Andrew Isaac
    """

    peak_point = get_maxima_indices(ic.get_intensity_array(), points)
    mlist = []
    for index in range(len(peak_point)):
        rt = ic.get_time_at_index(peak_point[index])
        intens = ic.get_intensity_at_index(peak_point[index])
        mlist.append([rt, intens])
    return mlist

def get_maxima_list_reduced(ic, mp_rt, points=13, window=3):

    """
    @summary: List of retention time and intensity of local maxima for ion
              Only peaks around a specific retention time are recorded
              created for use with gap filling algorithm
    @param ic: An IonChromatogram object
    @type ic: pyms.GCMS.Class.IonChromatogram
    @param mp_rt: The retention time of the missing peak
    @type ic: floatType
    @param points: Peak if maxima over 'points' number of scans
    @type points: IntType
    @param window: The window around the mp_rt where peaks should
                   be recorded
    @type window: intType

    @return: A list of retention time and intensity of local maxima for ion
    @rtype: ListType

    @author: Andrew Isaac
    """

    peak_point = get_maxima_indices(ic.get_intensity_array(), points)
    mlist = []
    for index in range(len(peak_point)):
        rt = ic.get_time_at_index(peak_point[index])
        if (rt > float(mp_rt) - window) and (rt < float(mp_rt) + window):
            intens = ic.get_intensity_at_index(peak_point[index])
            mlist.append([rt, intens])
        else:
            pass
    return mlist


def get_maxima_matrix(im, points=3, scans=1):

    """
    @summary: Get matrix of local maxima for each ion

    @param im: A list of intensities for a single ion
    @type im: ListType
    @param points: Peak if maxima over 'points' number of scans
    @type points: IntType
    @param scans: To compensate for spectra scewing,
        peaks from 'scans' scans are combined (Default 1).
    @type scans: IntType

    @return: A matrix of each ion and scan and intensity at ion peaks
    @rtype: ListType

    @author: Andrew Isaac
    """

    numrows, numcols = im.get_size()
    # zeroed matrix, size numrows*numcols
    maxima_im = numpy.zeros((numrows, numcols))
    raw_im = numpy.array(im.get_matrix_list())

    for col in range(numcols):  # assume all rows have same width
        # 1st, find maxima
        maxima = get_maxima_indices(raw_im[:,col], points)
        # 2nd, fill intensities
        for row in maxima:
            maxima_im[row, col] = raw_im[row, col]

    # combine spectra within 'scans' scans.
    half = int(scans/2)
    for row in range(numrows):
        tic = 0
        best = 0
        loc = 0
        # find best in scans
        for ii in range(scans):
            if row - half + ii >= 0 and row - half + ii < numrows:
                tic = maxima_im[row - half + ii].sum()
                # find largest tic of scans
                if tic > best:
                    best = tic
                    loc = ii
        # move and add others to best
        for ii in range(scans):
            if row - half + ii >= 0 and row - half + ii < numrows and ii != loc:
                for col in range(numcols):
                    maxima_im[row - half + loc, col] += \
                        maxima_im[row - half + ii, col]
                    maxima_im[row - half + ii, col] = 0

    return maxima_im
