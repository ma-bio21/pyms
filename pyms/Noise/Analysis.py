"""
Noise analysis functions
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

import random, math

from pyms.Utils.Error import error
from pyms.GCMS.Function import is_ionchromatogram
from pyms.Utils.Time import window_sele_points
from pyms.Utils.Math import MAD

_DEFAULT_WINDOW = 256
_DEFAULT_N_WINDOWS = 1024

def window_analyzer(ic, window=_DEFAULT_WINDOW, n_windows=_DEFAULT_N_WINDOWS, rand_seed=None ):

    """
    @summary: A simple estimator of the signal noise based on randomly
        placed windows and median absolute deviation

        The noise value is estimated by repeatedly and picking random
        windows (of a specified width) and calculating median absolute
        deviation (MAD). The noise estimate is given by the minimum MAD.

    @param ic: An IonChromatogram object
    @type ic: pyms.IO.Class.IonCromatogram
    @param window: Window width selection
    @type window: IntType or StringType
    @param n_windows: The number of windows to calculate
    @type n_windows: IntType
    @param rand_seed: Random seed generator
    @type rand_seed: IntType

    @return: The noise estimate
    @rtype: FloatType

    @author: Vladimir Likic
    """

    if not is_ionchromatogram(ic):
        error("argument must be an IonChromatogram object")

    ia = ic.get_intensity_array() # fetch the intensities

    # create an instance of the Random class
    if rand_seed != None:
        generator = random.Random(rand_seed)
    else:
        generator = random.Random()

    window_pts = window_sele_points(ic, window)

    maxi = ia.size - window_pts
    noise_level = math.fabs(ia.max()-ia.min())
    best_window_pos = None
    seen_positions = []

    cntr = 0

    while cntr < n_windows:
        # generator.randrange(): last point not included in range
        try_pos = generator.randrange(0, maxi+1)
        # only process the window if not analyzed previously
        if try_pos not in seen_positions:
            end_slice = try_pos + window_pts
            slice = ia[try_pos:end_slice]
            crnt_mad = MAD(slice)
            if crnt_mad < noise_level:
                noise_level = crnt_mad
                best_window_pos = try_pos
        cntr = cntr + 1
        seen_positions.append(try_pos)

    return noise_level

