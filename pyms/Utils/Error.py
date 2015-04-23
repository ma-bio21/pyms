"""
Defines error()
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

import sys

def error(message=None):

    """
    @summary: PyMS wide error function

    Prints out the error message (if supplied) and terminates execution
    with sys.exit(1).

    @param message: The error message to be printed
    @type message: StringType

    @return: none
    @rtype: NoneType

    @author: Lewis Lee
    @author: Vladimir Likic
    """

    sys.stdout = sys.__stderr__

    if message == None:
        message = "(no message)"
    else:
        message.rstrip('\n')

    # Retrieve details of the caller function from the caller
    # function's stack frame (i.e. the second most recent stack
    # frame, denoted by "1").
    funame = sys._getframe(1).f_code.co_name
    fargcount = sys._getframe(1).f_code.co_argcount
    fmodule = sys._getframe(1).f_code.co_filename
    fmodule = fmodule.split('/')
    fmodule = fmodule[-1]
    fmodule = fmodule.split('.')[0]

    if funame == '?':
        fmessage = "ERROR: 'main' in module %s.\n" % (fmodule)
    else:
        fmessage = "ERROR: %s(%d) in module %s.\n" % \
                (funame, fargcount, fmodule)

    fmessage = fmessage + " %s" % ( message )

    message_list = fmessage.split("\n")

    n = 0
    for line in message_list:
        if len(line) > n:
            n = len(line)

    cstr = ""
    for ii in range(n):
        cstr = cstr + "=" 
    print "\n %s" % (cstr)
    print " %s" % (fmessage)
    print " %s\n" % (cstr)

    sys.exit(1)

def stop(message=None):

    """
    @summary: A simple termination of execution

    @param message: The message to be printed
    @type message: StringType

    @return: none
    @rtype: NoneType
    """

    if message != None:
        print message 

    raise RuntimeError

