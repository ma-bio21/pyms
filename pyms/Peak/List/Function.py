"""
Functions related to Peak modification
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

import numpy, math

#from pyms.Utils.Error import error
from pyms.Peak.Class import Peak
from pyms.Utils.Utils import is_str, is_list
from pyms.Utils.Math import median
from pyms.GCMS.Class import MassSpectrum

# If psyco is installed, use it to speed up running time
try:
    import psyco
    psyco.full()
except:
    pass

def composite_peak(peak_list, minutes=False):

    """
    @summary: Create a peak that consists of a composite spectrum from all
        spectra in the list of peaks.

    @param peak_list: A list of peak objects
    @type peak_list: ListType
    @param minutes: Return retention time as minutes
    @type minutes: BooleanType

    @return: Peak Object with combined mass spectra of 'peak_list'
    @type: pyms.Peak.Class.Peak

    @author: Andrew Isaac
    """

    first = True
    count = 0
    avg_rt = 0
    new_ms = None
    for peak in peak_list:
        if peak is not None:
            ms = peak.get_mass_spectrum()
            spec = numpy.array(ms.mass_spec, dtype='d')
            if first:
                avg_spec = numpy.zeros(len(ms.mass_spec), dtype='d')
                mass_list = ms.mass_list
                first = False
            # scale all intensities to [0,100]
            max_spec = max(spec)/100.0
            if max_spec > 0:
                spec = spec/max_spec
            else:
                spec = spec*0
            avg_rt += peak.get_rt()
            avg_spec += spec
            count += 1
    if count > 0:
        avg_rt = avg_rt/count
        #if minutes == True:
            #avg_rt = avg_rt/60.0
        avg_spec = avg_spec/count
        avg_spec = avg_spec.tolist()  # list more compact than ndarray
        new_ms = MassSpectrum(mass_list, avg_spec)
        return Peak(avg_rt, new_ms, minutes)
    else:
        return None

def fill_peaks(data, peak_list, D, minutes=False):

    """
    @summary: Gets the best matching Retention Time and spectra from 'data' for
        each peak in the peak list.

    @param data: A data IntensityMatrix that has the same mass range as the
        peaks in the peak list
    @type data: pyms.GCMS.Class.IntensityMatrix
    @param peak_list: A list of peak objects
    @type peak_list: ListType
    @param D: Peak width standard deviation in seconds.  Determines search
        window width.
    @type D: FloatType
    @param minutes: Return retention time as minutes
    @type minutes: BooleanType

    @return: List of Peak Objects
    @type: ListType

    @author: Andrew Isaac
    """

    # Test for best match in range where RT weight is greater than _TOL
    _TOL = 0.001
    cutoff = D*math.sqrt(-2.0*math.log(_TOL))

    # Penalise for neighboring peaks
    # reweight so RT weight at nearest peak is _PEN
    _PEN = 0.5

    datamat = data.get_matrix_list()
    mass_list = data.get_mass_list()
    datatimes = data.get_time_list()
    minrt = min(datatimes)
    maxrt = max(datatimes)
    rtl = 0
    rtr = 0
    new_peak_list = []
    for ii in xrange(len(peak_list)):
        spec = peak_list[ii].get_mass_spectrum().mass_spec
        spec = numpy.array(spec, dtype='d')
        rt = peak_list[ii].get_rt()
        spec_SS = numpy.sum(spec**2, axis=0)

        # get neighbour RT's
        if ii > 0:
            rtl = peak_list[ii-1].rt
        if ii < len(peak_list)-1:
            rtr = peak_list[ii+1].rt
        # adjust weighting for neighbours
        rtclose = min(abs(rt-rtl), abs(rt-rtr))
        Dclose = rtclose/math.sqrt(-2.0*math.log(_PEN))

        if Dclose > 0:
            Dclose = min(D, Dclose)
        else:
            Dclose = D

        # Get bounds
        rtlow = rt - cutoff
        if rtlow < minrt:
            rtlow = minrt
        lowii = data.get_index_at_time(rtlow)

        rtup = rt + cutoff
        if rtup > maxrt:
            rtup = maxrt
        upii = data.get_index_at_time(rtup)

        # Get sub matrix of scans in bounds
        submat = datamat[lowii:upii+1]
        submat = numpy.array(submat, dtype='d')
        subrts = datatimes[lowii:upii+1]
        subrts = numpy.array(subrts, dtype='d')

        submat_SS = numpy.sum(submat**2, axis=1)

        # transpose spec (as matrix) for dot product
        spec = numpy.transpose([spec])
        # dot product on rows

        toparr = numpy.dot(submat, spec)
        botarr = numpy.sqrt(spec_SS*submat_SS)

        # convert back to 1-D array
        toparr = toparr.ravel()

        # scaled dot product of each scan
        cosarr = toparr/botarr

        # RT weight of each scan
        rtimearr = numpy.exp(-((subrts-rt)/float(Dclose))**2 / 2.0)

        # weighted scores
        scorearr = cosarr*rtimearr

        # index of best score
        best_ii = scorearr.argmax()

        # Add new peak
        bestrt = subrts[best_ii]
        bestspec = submat[best_ii].tolist()
        ms = MassSpectrum(mass_list, bestspec)
        new_peak_list.append(Peak(bestrt, ms, minutes))

    return new_peak_list
