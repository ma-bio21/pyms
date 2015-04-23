"""
Functions related to storing and loading a list of Peak objects
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

import string, cPickle

from pyms.Utils.Error import error
from pyms.Peak.Class import Peak
from pyms.Utils.Utils import is_str, is_list

def store_peaks(peak_list, file_name):

    """
    @summary:Store the list of peak objects

    @param peak_list: A list of peak objects
    @type peak_list: pyms.Peaks.Class.Peak
    @param file_name: File name to store peak list
    @type file_name: StringType

    @author: Andrew Isaac
    """

    if not is_str(file_name):
        error("'file_name' must be a string")

    fp = open(file_name,'w')
    cPickle.dump(peak_list, fp, 1)
    fp.close()

def load_peaks(file_name):

    """
    @summary: Loads the peak_list stored with 'store_peaks'

    @param file_name: File name of peak list
    @type file_name: StringType

    @return: The list of Peak objects
    @rtype: ListType

    @author: Andrew Isaac
    """

    if not is_str(file_name):
        error("'file_name' not a string")

    fp = open(file_name,'r')
    peak_list = cPickle.load(fp)
    fp.close()

    if not is_list(peak_list):
        error("'file_name' is not a List")
    if not len(peak_list) > 0 and not isinstance(peak_list[0], Peak):
        error("'peak_list' must be a list of Peak objects")

    return peak_list
