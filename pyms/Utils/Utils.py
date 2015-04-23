"""
General utility functions
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

import types, re

import numpy

from pyms.Utils.Error import error

def is_str(arg):

    """
    @summary: Returns True if the argument is a string, False otherwise

    @param arg: The argument to be evaluated as a string
    @type arg: arbitrary

    @return: A boolean indicator True or False
    @rtype: BooleanType

    @author: Vladimir Likic
    """

    if isinstance(arg,types.StringType):
        return True 
    else:
        return False

def is_int(arg):

    """
    @summary: Returns True if the argument is an integer, False
        otherwise

    @param arg: The argument to be evaluated as an integer
    @type arg: arbitrary

    @return: A boolean indicator True or False
    @rtype: BooleanType

    @author: Vladimir Likic
    """

    if isinstance(arg,types.IntType) or isinstance(arg,types.LongType):
        return True
    else:
        return False

def is_float(arg):

    """
    @summary: Returns True if the argument is a float, False otherwise

    @param arg: The argument to be evaluated as a float
    @type arg: arbitrary

    @return: A boolean indicator True or False
    @rtype: BooleanType

    @author: Vladimir Likic
    """

    if isinstance(arg,types.FloatType):
        return True
    else:
        return False

def is_number(arg):

    """
    @summary: Returns True if the argument is a number (integer or
        float), False otherwise
   
    @param arg: The argument to be evaluated as a number
    @type arg: arbitrary

    @return: A boolean indicator True or False
    @rtype: BooleanType

    @author: Vladimir Likic
    """

    if is_int(arg) or is_float(arg):
        return True 
    else:
        return False

def is_list(arg):

    """
    @summary: Returns True if the argument is a list, tuple, or numpy
        array, False otherwise

    @param arg: The argument to be evaluated as a list
    @type arg: arbitrary

    @return: A boolean indicator True or False
    @rtype: BooleanType

    @author: Vladimir Likic
    """

    if isinstance(arg,types.ListType) or isinstance(arg,types.TupleType) \
            or isinstance(arg, numpy.core.ndarray):
        return True 
    else:
        return False

def is_array(arg):

    """
    @summary: Returns True if the argument is a numpy array, False
        otherwise

    @param arg: The argument to be evaluated as a numpy array
    @type arg: arbitrary

    @return: A boolean indicator True or False
    @rtype: BooleanType 

    @author: Vladimir Likic
    """

    if isinstance(arg, numpy.core.ndarray):
        return True 
    else:
        return False


def is_boolean(arg):

    """
    @summary: Returns true of the argument is booleean, False otherwise

    @param arg: The argument to be evaluated as boolean
    @type arg: arbitrary

    @return: A boolean indicator True or False
    @rtype: BooleanType 

    @author: Vladimir Likic
    """

    if isinstance(arg,types.BooleanType):
        return True
    else:
        return False 

def is_str_num(arg):

    """
    @summary: Determines if the argument is a string in the format of a number

    The number can be an integer, or alternatively floating point in scientific
    or engineering format.

    @param arg: A string to be evaluate as a number
    @type arg: StringType

    @return: A boolean indicator True or False
    @rtype:  BooleanType

    @author: Gyro Funch (from Active State Python Cookbook)
    """

    NUM_RE = re.compile(r'^[-+]?([0-9]+\.?[0-9]*|\.[0-9]+)([eE][-+]?[0-9]+)?$')

    if NUM_RE.match(str(arg)):
        return True
    else:
        return False

def is_positive_int(arg):

    """
    @summary: Determines if the argument is an integer greater than zero

    @param arg: A string to be evaluate as a postive integer
    @type arg: types.StringType

    @return: A boolean indicator True or False
    @rtype:  BooleanType

    @author: Milica Ng
    """

    if not is_int(arg):
        return False
    elif not (arg > 0):
        return False
    else:
        return True

def is_list_of_dec_nums(arg):

    """
    @summary: Determines if the argument is a list of decimal numbers

    @param arg: A string to be evaluate as a list of decimal numbers
    @type arg: types.StringType

    @return: A boolean indicator True or False
    @rtype:  BooleanType

    @author: Milica Ng
    """

    if not(isinstance(arg, types.ListType)):
        return False
    elif (arg == []):
        return False
    else:
        for q in arg:
           if not(isinstance(q, types.FloatType)):
               return False
    return True

