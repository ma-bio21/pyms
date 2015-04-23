"""
General I/O functions
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

import types, os, string, cPickle

from pyms.Utils.Error import error
from pyms.Utils.Utils import is_number, is_str, is_list

def dump_object(object, file_name):

    """
    @summary: Dumps an object to a file through cPickle.dump()

    @param object: Object to be dumpted
    @type object: An instance of an arbitrary class

    @param file_name: Name of the file for the object dump
    @type file_name: StringType

    @author: Vladimir Likic
    """

    fp = open_for_writing(file_name)
    cPickle.dump(object, fp)
    close_for_writing(fp)

def load_object(file_name):

    """
    @summary: Loads an object previously dumped with dump_object()

    @param file_name: Name of the object dump file
    @type file_name: StringType

    @return: Object contained in the file 'file_name'
    @rtype: An instance of an arbitrary class

    @author: Vladimir Likic
    """

    fp = open_for_reading(file_name)
    object = cPickle.load(fp)
    close_for_reading(fp)

    return object

def open_for_reading(file_name):

    """
    @summary: Opens file for reading, returns file pointer

    @param file_name: Name of the file to be opened for reading
    @type file_name: StringType

    @return: Pointer to the opened file
    @rtype: FileType

    @author: Vladimir Likic
    """

    if not is_str(file_name):
        error("'file_name' is not a string")
    try:
        fp = open(file_name)
    except IOError:
        error("'%s' does not exist" % (file_name))

    return fp

def open_for_writing(file_name):

    """
    @summary: Opens file for writing, returns file pointer

    @param file_name: Name of the file to be opened for writing
    @type file_name: StringType

    @return: Pointer to the opened file
    @rtype: FileType

    @author: Vladimir Likic
    """

    if not is_str(file_name):
        error("'file_name' is not a string")
    try:
        fp = open(file_name, "w")
    except IOError:
        error("Cannot open '%s' for writing" % (file_name))

    return fp

def close_for_reading(fp):

    """
    @summary: Closes file pointer open for reading

    @param fp: A file pointer, previously opened for reading
    @type fp: FileType

    @return: none
    @rtype: NoneType

    @author: Vladimir Likic
    """

    fp.close()

def close_for_writing(fp):

    """
    @summary: Closes file pointer open for writing

    @param fp: A file pointer, previously opened for writing
    @type fp: FileType

    @return: none
    @rtype: NoneType

    @author: Vladimir Likic
    """

    fp.close()

def file_lines(file_name, filter=False):

    """
    @summary: Returns lines from a file, as a list

    @param file_name: Name of a file
    @type: StringType
    @param filter: If True, lines are pre-processes. Newline character
        if removed, leading and taling whitespaces are removed, and lines
        starting with '#' are discarded
    @type: BooleanType 

    @return: A list of lines
    @rtype: ListType

    @author: Vladimir Likic
    """

    if not is_str(file_name):
        error("'file_name' is not a string")

    fp = open_for_reading(file_name)
    lines = fp.readlines()
    close_for_reading(fp)

    if filter:
        # strip leading and talining whitespaces
        lines_filtered = []
        for line in lines:
            line = line.strip()
            lines_filtered.append(line)

        # discard comments
        lines_to_discard = []
        for line in lines_filtered:
            # remove empty lines and comments
            if len(line) == 0 or line[0] == "#":
                lines_to_discard.append(line)
        for line in lines_to_discard:
            lines_filtered.remove(line)
        lines = lines_filtered

    return lines

def save_data(file_name, data, format_str="%.6f", prepend="", sep=" ",
	compressed=False):

    """
    @summary: Saves a list of numbers or a list of lists of numbers
    to a file with specific formatting

    @param file_name: Name of a file
    @type: StringType
    @param data: A list of numbers, or a list of lists
    @type: ListType
    @param format_str: A format string for individual entries
    @type: StringType
    @param prepend: A string, printed before each row
    @type: StringType
    @param sep: A string, printed after each number
    @type: StringType
    @param compressed: A boolean. If True, the output will be gzipped
    @type: BooleanType

    @return: none
    @rtype: NoneType

    @author: Vladimir Likic
    """

    if not is_str(file_name):
        error("'file_name' is not a string")

    if not is_list(data):
        error("'data' is not a list")

    if not is_str(prepend):
        error("'prepend' is not a string")

    if not is_str(sep):
        error("'sep' is not a string")

    fp = open_for_writing(file_name)

    # decide whether data is a vector or matrix
    if is_number(data[0]):
        for item in data:
            if not is_number(item):
                error("not all elements of the list are numbers")
        data_is_matrix = 0
    else:
        for item in data:
            if not is_list(item):
                error("not all elements of the list are lists")
        data_is_matrix = 1

    if data_is_matrix:
        for ii in range(len(data)):
            fp.write(prepend)
            for jj in range(len(data[ii])):
                if is_number(data[ii][jj]):
                    fp.write(format_str % (data[ii][jj]))
                    if (jj<(len(data[ii])-1)): fp.write(sep)
                else:
                    error("datum not a number")
            fp.write("\n")
    else:
        for ii in range(len(data)):
            fp.write(prepend)
            fp.write(format_str % (data[ii]))
            fp.write("\n")

    close_for_writing(fp)

    if compressed:
        status = os.system('gzip %s' % (file_name))
        if status != 0:
            error("gzip compress failed")

