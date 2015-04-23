"""
Models a GC-MS experiment represented by a list of signal peaks
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

from pyms.Utils.Error import error
from pyms.Utils.Utils import is_str
from pyms.Peak.Class import Peak
from pyms.Peak.List.Utils import is_peak_list, sele_peaks_by_rt

class Experiment:

    """
    @summary: Models an experiment object

    @author: Vladimir Likic
    @author: Andrew Isaac
    """

    def __init__(self, expr_code, peak_list):

        """
        @summary: Models an experiment

        @param expr_code: Unique identifier for the experiment
        @type expr_code: StringType
        @param peak_list: A list of peak objects
        @type peak_list: ListType
        """

        if not is_str(expr_code):
            error("'expr_code' must be a string")
        if not is_peak_list(peak_list):
            error("'peak_list' must be a list of Peak objects")

        self.__expr_code = expr_code
        self.__peak_list = peak_list

    def get_expr_code(self):

        """
        @summary: Returns the expr_code of the experiment

        @return: The expr_code of the experiment
        @rtype:  StringType
        """

        return self.__expr_code

    def get_peak_list(self):

        """
        @summary: Returns the peak list

        @return: A list of peak objects
        @rtype: ListType
        """

        return self.__peak_list

    def sele_rt_range(self, rt_range):

        """
        @summary: Discards all peaks which have the retention time outside
        the specified range

        @param rt_range: Min, max retention time given as a list [rt_min,rt_max]
        @type rt_range: ListType

        @return: none
        @rtype: NoneType
        """

        peaks_sele = sele_peaks_by_rt(self.__peak_list, rt_range)
        self.__peak_list = peaks_sele
