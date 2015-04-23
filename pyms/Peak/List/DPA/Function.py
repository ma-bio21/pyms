"""
Functions for peak alignment by dynamic programming
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

import copy
import numpy
#from math import sqrt, log
import math

from pyms.Utils.Error import error, stop
from pyms.Utils.Utils import is_list
from pyms.Utils.DP import dp
from pyms.Experiment.Class import Experiment

import Class
import Utils

# If psyco is installed, use it to speed up running time
try:
    import psyco
    psyco.full()
except:
    pass

def align_with_tree(T, min_peaks=1):

    """
    @summary: Aligns a list of alignments using the supplied guide tree

    @param T: The pairwise alignment object
    @type: pyms.Peak.List.DPA.Class.PairwiseAlignment
    @return: The final alignment consisting of aligned input alignments
    @rtype: pyms.Peak.List.DPA.Class.Alignment
    @author: Woon Wai Keen
    @author: Vladimir Likic
    """

    print " Aligning %d items with guide tree (D=%.2f, gap=%.2f)" % \
            (len(T.algts), T.D, T.gap)

    # For everything else, we align according to the guide tree provided by
    # Pycluster. From Pycluster documentation:
    #   Each item and subnode is represented by an integer. For hierarchical
    #   clustering of n items, we number the original items {0, ... , n-1},
    #   nodes are numbered {-1, ... , -(n-1)}. Note that the number of nodes
    #   is one less than the number of items.

    # extend As to length 2n to hold the n items, n-1 nodes, and 1 root
    As = copy.deepcopy(T.algts) + [ None for _ in range(len(T.algts)) ]

    # align the alignments into positions -1, ... ,-(n-1)
    total = len(T.tree)
    index = 0

    for node in T.tree:
        index = index - 1
        As[index] = align(As[node.left], As[node.right], T.D, T.gap)
        total = total - 1
        print " -> %d item(s) remaining" % total

    # the final alignment is in the root. Filter min peaks and return
    final_algt =  As[index]

    # useful for within state alignment only
    if min_peaks > 1:
        final_algt.filter_min_peaks(min_peaks)

    return final_algt

def exprl2alignment(exprl):

    """
    @summary: Converts experiments into alignments

    @param exprl: The list of experiments to be converted into an alignment
        objects
    @type exprl: ListType

    @author: Vladimir Likic
    """

    if not is_list(exprl):
        error("the argument is not a list")

    algts = []

    for item in exprl:
        if not isinstance(item, Experiment):
            error("list items must be 'Experiment' instances")
        else:
            algt = Class.Alignment(item)
        algts.append(algt)

    return algts

def align(a1, a2, D, gap):

    """
    @summary: Aligns two alignments

    @param a1: The first alignment
    @type a1: pyms.Peak.List.Class.Alignment
    @param a2: The second alignment
    @type a2: pyms.Peak.List.Class.Alignment
    @param D: Retention time tolerance
    @type D: FloatType
    @param gap: Gap penalty
    @type gap: FloatType

    @return: Aligned alignments
    @rtype: pyms.Peak.List.Class.Alignment

    @author: Woon Wai Keen
    @author: Vladimir Likic
    """

    # calculate score matrix for two alignments
    M = score_matrix(a1, a2, D)

    # run dynamic programming
    result = dp(M, gap)

    # make composite alignment from the results
    ma = merge_alignments(a1, a2, result['trace'])

    # calculate the similarity score
    ma.similarity = alignment_similarity(result['trace'], M, gap)

    return ma

def merge_alignments(A1, A2, traces):

    """
    @summary: Merges two alignments with gaps added in from DP traceback

    @param A1: First alignment
    @param A2: Second alignment
    @param traces: DP traceback

    @return: A single alignment from A1 and A2
    @author: Woon Wai Keen
    @author: Vladimir Likic
    @author: Qiao Wang
    """

    # Create object to hold new merged alignment and fill in its expr_codes
    ma = Class.Alignment(None)
    ma.expr_code = A1.expr_code + A2.expr_code

    # create empty lists of dimension |A1| + |A2|
    dimension = len(A1.peakpos) + len(A2.peakpos)
    merged = [ [] for _ in range(dimension) ]
    A1 = A1.peakpos
    A2 = A2.peakpos

    idx1 = idx2 = 0

    # trace can either be 0, 1, or 2
    # if it is 0, there are no gaps. otherwise, if it is 1 or 2,
    # there is a gap in A2 or A1 respectively.

    for trace in traces:

        if trace == 0:
            for i in range(len(A1)):
                merged[i].append(A1[i][idx1])

            for j in range(len(A2)):
                merged[1+i+j].append(A2[j][idx2])

            idx1 = idx1 + 1
            idx2 = idx2 + 1

        elif trace == 1:
            for i in range(len(A1)):
                merged[i].append(A1[i][idx1])

            for j in range(len(A2)):
                merged[1+i+j].append(None)

            idx1 = idx1 + 1

        elif trace == 2:
            for i in range(len(A1)):
                merged[i].append(None)

            for j in range(len(A2)):
                merged[1+i+j].append(A2[j][idx2])

            idx2 = idx2 + 1

    ma.peakalgt = numpy.transpose(merged)
    # sort according to average peak
    ma.peakalgt = list(ma.peakalgt)
    ma.peakalgt.sort(Utils.alignment_compare)
    ma.peakpos = numpy.transpose(ma.peakalgt)

    return ma

def alignment_similarity(traces, score_matrix, gap):

    """
    @summary: Calculates similarity score between two alignments (new method)

    @param traces: Traceback from DP algorithm
    @param score_matrix: Score matrix of the two alignments
    @param gap: Gap penalty

    @return: Similarity score (i.e. more similar => higher score)
    @author: Woon Wai Keen
    @author: Vladimir Likic
    """

    score_matrix = 1. - score_matrix
    similarity = 0.
    idx1 = idx2 = 0

    # Trace can either be 0, 1, or 2
    # If it is 0, there is a match and we add to the sum the score between
    # these two aligned peaks.
    #
    # Otherwise, if it is 1 or 2, and there is a gap in A2 or A1
    # respectively. We then subtract the gap penalty from the sum.
    for trace in traces:
        if trace == 0:
            similarity = similarity + score_matrix[idx1][idx2]
            idx1 = idx1 + 1
            idx2 = idx2 + 1
        elif trace == 1:
            similarity = similarity - gap
            idx1 = idx1 + 1
        elif trace == 2:
            similarity = similarity - gap
            idx2 = idx2 + 1

    return similarity

def score_matrix(a1, a2, D):

    """
    @summary: Calculates the score matrix between two alignments

    @param a1: The first alignment
    @type a1: pyms.Peak.List.Class.Alignment
    @param a2: The second alignment
    @type a2: pyms.Peak.List.Class.Alignment
    @param D: Retention time tolerance
    @type D: FloatType

    @return: Aligned alignments
    @rtype: pyms.Peak.List.Class.Alignment

    @author: Qiao Wang
    @author: Andrew Isaac
    """

    score_matrix = numpy.zeros((len(a1.peakalgt), len(a2.peakalgt)))
    row = 0
    col = 0
    sim_score = 0

    for algt1pos in a1.peakalgt:
        for algt2pos in a2.peakalgt:
            sim_score = position_similarity(algt1pos, algt2pos, D)
            score_matrix[row][col] = sim_score
            col = col+1
        row = row+1
        col = 0

    return score_matrix

def position_similarity(pos1, pos2, D):

    """
    @summary: Calculates the similarity between the two alignment
        positions.  A score of 0 is best and 1 is worst.

    @param pos1: The position of the first alignment
    @param pos2: The position of the second alignment
    @param D: Rentention time tolerance

    @return: The similarity value for the current position
    @rtype: FloatType

    @author: Qiao Wang
    @author: Vladimir Likic
    @author: Andrew Isaac
    """

    score = 0.0
    count = 0

    ## Attempt to speed up by only calculating 'in-range' values
    ## set tollerance to 1/1000
    _TOL = 0.001
    cutoff = D*math.sqrt(-2.0*math.log(_TOL))

    for a in pos1:
        if a is not None:
            aspec = a.mass_spec
            art = a.rt
            once = True
            for b in pos2:
                if b is not None:
                    brt = b.rt
                    # in range?
                    if abs(art-brt) > cutoff:
                        score += 1.0  # NB score of 1 is worst
                    else:
                        # Once per b-loop
                        if once:
                            mass_spect1 = numpy.array(aspec, dtype='d')
                            mass_spect1_sum = numpy.sum(mass_spect1**2, axis=0)
                            once = False
                        bspec = b.mass_spec
                        mass_spect2 = numpy.array(bspec, dtype='d')
                        mass_spect2_sum = numpy.sum(mass_spect2**2, axis=0)
                        try:
                            top = numpy.dot(mass_spect1, mass_spect2)
                        except(ValueError):
                            error("Mass Spectra are of different length\n\n" +  
                                 " Use IntensityMatrix.crop_mass() to set\n" 
                                  + " same length for all Mass Spectra""")
                        bot = numpy.sqrt(mass_spect1_sum*mass_spect2_sum)
                        if bot > 0:
                            cos = top/bot
                        else:
                            cos = 0
                        rtime = numpy.exp(-((art-brt)/float(D))**2 / 2.0)
                        score = score + (1.0 - (cos*rtime))
                    count = count + 1

    if count == 0:
        score = 1.0  # NB score of 1 is worst
    else:
        score = score/float(count)

    return score
