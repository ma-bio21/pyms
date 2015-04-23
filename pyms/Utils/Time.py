"""
Time conversion and related functions
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

from pyms.Utils.Error import error
from pyms.Utils.Utils import is_int, is_str, is_str_num

def time_str_secs(time_str):

    """
    @summary: Resolves time string of the form "<NUMBER>s" or "<NUMBER>m",
        returns time in seconds

    @param time_str: A time string, which must be of the form
        "<NUMBER>s" or "<NUMBER>m" where "<NUMBER>" is a valid number
    @type time_str: StringType

    @return: Time in seconds
    @rtype: FloatType

    @author: Vladimir Likic
    """

    if not is_str(time_str):
        error("time string not a string")

    time_number = time_str[:-1]
    time_spec = time_str[-1].lower()

    if not is_str_num(time_number):
       print " --> received time string '%s'" % (time_number)
       error("improper time string")

    if not time_spec == "s" and not time_spec == "m":
        error("time string must end with either 's' or 'm'")

    time = float(time_number)

    if time_spec == "m":
        time = time*60.0

    return time

def window_sele_points(ic, window_sele, half_window=False):

    """
    @summary: Converts window selection parameter into points based
        on the time step in an ion chromatogram

    @param ic: ion chromatogram object relevant for the conversion
    @type ic: pyms.GCMS.Class.IonChromatogram

    @param window_sele: The window selection parameter. This can be
        an integer or time string. If integer, taken as the number
        of points. If a string, must of the form "<NUMBER>s" or
        "<NUMBER>m", specifying a time in seconds or minutes,
        respectively
    @type window_sele: IntType or StringType

    @param half_window: Specifies whether to return half-window
    @type half_window: BooleanType

    @return: The number of points in the window
    @rtype: IntType

    @author: Vladimir Likic
    """

    if not is_int(window_sele) and not is_str(window_sele):
        error("'window' must be an integer or a string")

    if is_int(window_sele):
        if half_window:
            if window_sele % 2 == 0:
                error("window must be an odd number of points")
            else:
                points = int(math.floor(window_sele*0.5))
        else:
            points = window_sele
    else:
        time = time_str_secs(window_sele)
        time_step = ic.get_time_step()

        if half_window:
            time = time*0.5

        points = int(math.floor(time/time_step))

    if half_window:
        if points < 1: error("window too small (half window=%d)" % (points))
    else:
        if points < 2: error("window too small (window=%d)" % (points))

    return points

