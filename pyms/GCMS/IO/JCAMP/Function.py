"""
Functions for I/O of data in JCAMP-DX format
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

from pyms.GCMS.Class import GCMS_data
from pyms.GCMS.Class import Scan
from pyms.Utils.IO import file_lines
from pyms.Utils.Utils import is_str
from pyms.Utils.Error import error

import numpy

def JCAMP_reader(file_name):

    """
    @summary: Generic reader for JCAMP DX files, produces GC-MS data
       object

    @author: Qiao Wang
    @author: Andrew Isaac
    @author: Vladimir Likic
    """

    if not is_str(file_name):
        error("'file_name' not a string")

    print " -> Reading JCAMP file '%s'" % (file_name)
    lines_list = open(file_name,'r')
    data = []
    page_idx = 0
    xydata_idx = 0
    time_list = []
    scan_list = []

    for line in lines_list:
        if not len(line.strip()) == 0:
            prefix = line.find('#')
            # key word or information
            if prefix == 0:
                fields = line.split('=')
                if fields[0].find("##PAGE") >= 0:
                    time = float(fields[2].strip()) #rt for the scan to be submitted
                    time_list.append(time)
                    page_idx = page_idx + 1
                elif fields[0].find("##DATA TABLE") >= 0:
                    xydata_idx = xydata_idx + 1
            # data
            elif prefix == -1:
                if page_idx > 1 or xydata_idx > 1:
                    if len(data) % 2 == 1:
                        error("data not in pair !")
                    mass = []
                    intensity = []
                    for i in range(len(data) / 2):
                        mass.append(data[i * 2])
                        intensity.append(data[i * 2 + 1])
                    if not len(mass) == len(intensity):
                        error("len(mass) is not equal to len(intensity)")
                    scan_list.append(Scan(mass, intensity))
                    data = []
                    data_sub = line.strip().split(',')
                    for item in data_sub:
                        if not len(item.strip()) == 0:
                            data.append(float(item.strip()))
                    if page_idx > 1:
                        page_idx = 1
                    if xydata_idx > 1:
                        xydata_idx = 1
                else:
                    data_sub = line.strip().split(',')
                    for item in data_sub:
                        if not len(item.strip()) == 0:
                            data.append(float(item.strip()))

    if len(data) % 2 == 1:
        error("data not in pair !")
    # get last scan
    mass = []
    intensity = []
    for i in range(len(data) / 2):
        mass.append(data[i * 2])
        intensity.append(data[i * 2 + 1])

    if not len(mass) == len(intensity):
        error("len(mass) is not equal to len(intensity)")
    scan_list.append(Scan(mass, intensity))

    # sanity check
    if not len(time_list) == len(scan_list):
        error("number of time points does not equal the number of scans")

    data = GCMS_data(time_list, scan_list)

    return data
