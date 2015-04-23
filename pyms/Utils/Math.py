"""
Provides mathematical functions
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

import copy, math

from pyms.Utils.Error import error
from pyms.Utils.Utils import is_list, is_number

def median(v):

    """
    @summary: Returns a median of a list or numpy array

    @param v: Input list or array
    @type v: ListType or numpy.core.ndarray
    @return: The median of the input list
    @rtype: FloatType

    @author: Vladimir Likic
    """

    if not is_list(v):
        error("argument neither list nor array")

    local_data = copy.deepcopy(v)
    local_data.sort()
    N = len(local_data)

    if (N % 2) == 0:
        # even number of points
        K = N/2 - 1 
        median = (local_data[K] + local_data[K+1])/2.0
    else:
	    # odd number of points
        K = (N - 1)/2 - 1
        median = local_data[K+1]

    return median

def vector_by_step(vstart,vstop,vstep):

    """
    @summary: generates a list by using start, stop, and step values

    @param vstart: Initial value 
    @type vstart: A number
    @param vstop: Max value
    @type vstop: A number
    @param vstep: Step
    @type vstep: A number
   
    @return: A list generated
    @rtype: ListType

    @author: Vladimir Likic
    """

    if not is_number(vstart) or not is_number(vstop) or not is_number(vstep):
        error("parameters start, stop, step must be numbers")

    v = []

    p = vstart 
    while p < vstop:
        v.append(p)
        p = p + vstep

    return v

def MAD(v):

    """
    @summary: median absolute deviation

    @param v: A list or array
    @type v: ListType, TupleType, or numpy.core.ndarray

    @return: median absolute deviation
    @rtype: FloatType

    @author: Vladimir Likic
    """

    if not is_list(v):
        error("argument neither list nor array")

    m = median(v)
    m_list = []

    for xi in v:
        d = math.fabs(xi - m)
        m_list.append(d)

    mad = median(m_list)/0.6745

    return mad

def amin(v):

    """
    @summary: Finds the minimum element in a list or array

    @param v: A list or array
    @type v: ListType, TupleType, or numpy.core.ndarray

    @return: Tuple (maxi, maxv), where maxv is the minimum 
        element in the list and maxi is its index
    @rtype: TupleType

    @author: Vladimir Likic
    """

    if not is_list(v):
        error("argument neither list nor array")

    minv = max(v) # built-in max() function
    mini = None

    for ii in range(len(v)):
        if v[ii] < minv:
            minv = v[ii]
            mini = ii

    if mini == None:
        error("finding maximum failed")

    return mini, minv

def mean(v):

    """
    @summary: Calculates the mean

    @param v: A list or array
    @type v: ListType, TupleType, or numpy.core.ndarray

    @return: Mean
    @rtype: FloatType

    @author: Vladimir Likic
    """

    if not is_list(v):
        error("argument neither list nor array")

    s = 0.0
    for e in v:
        s = s + e 
    s_mean = s/float(len(v))

    return s_mean

def std(v):

    """
    @summary: Calculates standard deviation

    @param v: A list or array
    @type v: ListType, TupleType, or numpy.core.ndarray

    @return: Mean
    @rtype: FloatType

    @author: Vladimir Likic
    """

    if not is_list(v):
        error("argument neither list nor array")

    v_mean = mean(v)

    s = 0.0 
    for e in v:
        d = e - v_mean
        s = s + d*d
    s_mean = s/float(len(v)-1)
    v_std = math.sqrt(s_mean)

    return v_std


def rmsd(list1, list2):

    """
    @summary: Calculates RMSD for the 2 lists

    @param list1: First data set
    @type list1: ListType, TupleType, or numpy.core.ndarray 
    @param list2: Second data set
    @type list2: ListType, TupleType, or numpy.core.ndarray 
    @return: RMSD value
    @rtype: FloatType

    @author: Qiao Wang
    @author: Andrew Isaac
    @author: Vladimir Likic
    """

    if not is_list(list1):
        error("argument neither list nor array")

    if not is_list(list2):
        error("argument neither list nor array")

    sum = 0.0
    for i in range(len(list1)):
        sum = sum + (list1[i] - list2[i]) ** 2
    rmsd = math.sqrt(sum / len(list1))
    return rmsd

