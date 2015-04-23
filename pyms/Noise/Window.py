"""
Moving window noise filter
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

import copy
import numpy

from pyms.GCMS.Function import is_ionchromatogram, ic_window_points
from pyms.Utils.Math import median

__DEFAULT_WINDOW = 3

def window_smooth(ic, window=__DEFAULT_WINDOW, median=False):

    """
    @summary: Applies window smoothing on ion chromatogram

    @param ic: The input ion chromatogram
    @type ic: pyms.GCMS.Class.IonChromatogram
    @param window: The window selection parameter. This can be an integer
        or time string. If integer, taken as the number of points. If a
        string, must of the form "<NUMBER>s" or "<NUMBER>m", specifying
        a time in seconds or minutes, respectively
    @type window: IntType or StringType
    @param median: An indicator whether the mean or median window smoothing
        to be used
    @type median: Booleantype

    @return: Smoothed ion chromatogram
    @rtype: pyms.GCMS.Class.IonChromatogram

    @author: Vladimir Likic
    """

    if not is_ionchromatogram(ic):
        error("'ic' not an IonChromatogram object")

    ia = ic.get_intensity_array()

    wing_length = ic_window_points(ic, window, half_window=True)

    if median:
        ia_denoise = __median_window(ia, wing_length)
    else:
        ia_denoise = __mean_window(ia, wing_length)

    ic_denoise = copy.deepcopy(ic)
    ic_denoise.set_intensity_array(ia_denoise)

    return ic_denoise

def window_smooth_im(im, window=__DEFAULT_WINDOW, median=False):
    """
    @summary: Applies window smoothing on Intensity Matrix

              Simply wraps around the window smooth function above

    @param im: The input Intensity Matrix
    @type im: pyms.GCMS.Class.IntensityMatrix
    @param window: The window selection parameter. 
    @type window: IntType or StringType
    
    @param median: An indicator whether the mean or median window smoothing
        to be used
    @type median: Booleantype

    @return: Smoothed Intensity Matrix
    @rtype: pyms.GCMS.Class.IntensityMatrix

    @author: Sean O'Callaghan
    @author: Vladimir Likic
    """
    
    n_scan, n_mz = im.get_size()
    
    im_smooth = copy.deepcopy(im)
    
    for ii in range(n_mz):
        ic = im_smooth.get_ic_at_index(ii)
        ic_smooth = window_smooth(ic, window, median)
        im_smooth.set_ic_at_index(ii, ic_smooth)
        
    return im_smooth

def __mean_window(ia, wing_length):

    """
    @summary: Applies mean-window averaging on the array of intensities.

    @param ia: Intensity array
    @type ia: nympy.core.ndarray
    @param wing_length: An integer value representing the number of
        points on either side of a point in the ion chromatogram
    @type wing_length: IntType

    @return: Smoothed intensity array
    @rtype: nympy.core.ndarray

    @author: Vladimir Likic
    """

#print " -> Window smoothing (mean): the wing is %d point(s)" % (wing_length)

    ia_denoise = numpy.repeat([0], ia.size)

    index = 0
    end = ia.size - 1

    while index <= end:
        left = index - wing_length
        right = index + wing_length + 1
        if left < 0: left = 0
        slice = ia[left:right]
        ia_denoise[index] = slice.mean()
        index = index + 1

    return ia_denoise

def __median_window(ia, wing_length):

    """
    @summary: Applies median-window averaging on the array of intensities.

    @param ia: Intensity array
    @type ia: nympy.core.ndarray
    @param wing_length: An integer value representing the number of
        points on either side of a point in the ion chromatogram
    @type wing_length: IntType

    @return: Smoothed intensity array
    @rtype: nympy.core.ndarray

    @author: Vladimir Likic
    """

#print " -> Window smoothing (median): the wing is %d point(s)" % (wing_length)

    ia_denoise = numpy.repeat([0], ia.size)

    index = 0
    end = ia.size - 1

    while index <= end:
        left = index - wing_length
        right = index + wing_length + 1
        if left < 0: left = 0
        slice = ia[left:right]
        ia_denoise[index] = median(slice)
        index = index + 1

    return ia_denoise

