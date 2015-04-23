"""
Functions related to experiment input/output
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
from pyms.Experiment.Class import Experiment
from pyms.Utils.Utils import is_str
from pyms.Utils.IO import file_lines
from pyms.GCMS.Class import IntensityMatrix

def load_expr(file_name):

    """
    @summary: Loads an experiment saved with 'store_expr'

    @param file_name: Experiment file name
    @type file_name: StringType

    @return: The experiment intensity matrix and peak list
    @rtype: pyms.Experiment.Class.Experiment

    @author: Vladimir Likic
    @author: Andrew Isaac
    """

    if not is_str(file_name):
        error("'file_name' not a string")

    fp = open(file_name,'r')
    expr = cPickle.load(fp)
    fp.close()

    if not isinstance(expr, Experiment):
        error("'file_name' is not an Experiment object")

    return expr

def store_expr(file_name, expr):

    """
    @summary: stores an expriment to a file

    @param file_name: The name of the file
    @type file_name: StringType
    @param expr: An experiment object
    @type expr: pyms.Experiment.Class.Experiment

    @return: none
    @rtype: NoneType

    @author: Vladimir Likic
    @author: Andrew Isaac
    """

    if not isinstance(expr, Experiment):
        error("argument not an instance of the class 'Experiment'")

    if not is_str(file_name):
        error("'file_name' not a string")

    fp = open(file_name,'w')
    cPickle.dump(expr, fp, 1)
    fp.close()

def read_expr_list(file_name):

    """
    @summary: Reads the set of experiment files and returns a list of
    Experiment objects

    @param file_name: The name of the file which lists experiment
        dump file names, one file per line
    @type file_name: StringType

    @return: A list of Experiment instances
    @rtype: ListType

    @author: Vladimir Likic
    """

    if not is_str(file_name):
        error("file_name argument must be a string")
    try:
        fp = open(file_name, 'r')
    except IOError:
        error("error opening file '%s' for reading" % file_name)

    exprfiles = fp.readlines()
    fp.close()

    exprl = []

    for exprfile in exprfiles:

        exprfile = string.strip(exprfile)
        expr = load_expr(exprfile)

        exprl.append(expr)

    return exprl
