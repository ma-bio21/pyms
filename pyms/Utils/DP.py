"""
Dynamic Programming routine
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

import numpy
from Error import error
#from mpi4py import MPI

def dp(S, gap_penalty):
    
    """ 
    @summary: Solves optimal path in score matrix based on global sequence
    alignment

    @param S: Score matrix
    @type S: numpy.
    @param gap_penalty: Gap penalty
    @type gap_penalty: FloatType

    @return: A dictionary of results
    @rtype: DictType

    @author: Tim Erwin
    """
    #comm = MPI.COMM_WORLD
    #rank = comm.Get_rank()
    #print " In DP.py, I am rank", rank

    try:
        row_length = len(S[:,0])
    except(IndexError):
        error('Zero length alignment found: Samples with no peaks \
               cannot be aligned')
    col_length = len(S[0,:])

    #D contains the score of the optimal alignment
    D = numpy.zeros((row_length+1,col_length+1), dtype='d')
    for i in range(1, row_length+1):
        D[i,0] = gap_penalty*i
    for j in range(1, col_length+1):
        D[0,j] = gap_penalty*j
    D[0,0] = 0.0
    D[1:(row_length+1), 1:(col_length+1)] = S.copy();

    # Directions for trace
    # 0 - match               (move diagonal)
    # 1 - peaks1 has no match (move up)
    # 2 - peaks2 has no match (move left)
    # 3 - stop
    trace_matrix = numpy.zeros((row_length+1,col_length+1))
    trace_matrix[:,0] = 1; 
    trace_matrix[0,:] = 2;
    trace_matrix[0,0] = 3;
   
    for i in range(1,row_length+1):
        for j in range(1,col_length+1):
 
            #
            # Needleman-Wunsch Algorithm assuming a score function S(x,x)=0
            #
            #              | D[i-1,j-1] + S(i,j)
            # D[i,j] = min | D(i-1,j] + gap
            #              | D[i,j-1] + gap
            #

            darray = [D[i-1,j-1]+S[i-1,j-1], D[i-1,j]+gap_penalty, D[i,j-1]+gap_penalty]
            D[i,j] = min(darray)
            #Store direction in trace matrix
            trace_matrix[i,j] = darray.index(D[i,j])

    # Trace back from bottom right
    trace = []
    matches = []
    i = row_length
    j = col_length
    direction = trace_matrix[i,j]
    p = [row_length-1]
    q = [col_length-1]
    
    while direction != 3:
        
        if direction == 0: #Match
            i = i-1
            j = j-1
            matches.append([i,j])
        elif direction == 1: #peaks1 has no match
            i = i-1
        elif direction == 2: #peaks2 has no match
            j = j-1
        p.append(i-1)
        q.append(j-1)
        trace.append(direction)
        direction=trace_matrix[i,j]

    #remove 'stop' entry
    p.pop()
    q.pop()
    # reverse the trace back
    p.reverse()
    q.reverse()
    trace.reverse()
    matches.reverse()

    return {'p':p, 'q':q, 'trace':trace, 'matches':matches, 'D':D, 'phi':trace_matrix}

