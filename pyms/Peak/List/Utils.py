"""
Utilities for manipulation of peak lists
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

import math

from pyms.Peak.Class import Peak
from pyms.Utils.Error import error
from pyms.Utils.Utils import is_list, is_str
from pyms.Utils.Time import time_str_secs

def is_peak_list(peaks):

    """
    @summary: Returns True if 'peaks' is a valid peak list, False
    otherwise

    @param peaks: A list of peak objects
    @type peaks: ListType

    @return: A boolean indicator
    @rtype: BooleanType

    @author: Vladimir Likic
    """

    flag = True

    if not is_list(peaks):
        flag = False
    else:
        for item in peaks:
           if not isinstance(item, Peak):
              flag = False

    return flag

def sele_peaks_by_rt(peaks, rt_range):

    """
    @summary: Selects peaks from a retention time range

    @param peaks: A list of peak objects
    @type peaks: ListType
    @param rt_range: A list of two time strings, specifying lower and
           upper retention times
    @type rt_range: ListType
    @return: A list of peak objects
    @rtype: ListType
    """

    if not is_peak_list(peaks):
        error("'peaks' not a peak list")

    if not is_list(rt_range):
        error("'rt_range' not a list")
    else:
        if len(rt_range) != 2:
            error("'rt_range' must have exactly two elements")

        if not is_str(rt_range[0]) or not is_str(rt_range[1]):
            error("lower/upper retention time limits must be strings")

    rt_lo = time_str_secs(rt_range[0])
    rt_hi = time_str_secs(rt_range[1])

    if not rt_lo < rt_hi:
        error("lower retention time limit must be less than upper")

    peaks_sele = []

    for peak in peaks:
        rt = peak.get_rt()
        if rt > rt_lo and rt < rt_hi:
            peaks_sele.append(peak)

    #print "%d peaks selected" % (len(peaks_sele))

    return peaks_sele
