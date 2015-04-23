"""
Provides functions for simulation of GCMS data
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
import math
import sys

sys.path.append("/x/PyMS/")
from pyms.GCMS.Class import IonChromatogram, IntensityMatrix


def gcms_sim(time_list, mass_list, peak_list):
    
    """
    @summary: Simulator of GCMS data
    
    @param time_list: the list of scan times
    @type time_list: listType
    
    @param mass_list: the list of m/z channels
    @type mass_list: listType
    
    @param peak_list: A list of peaks
    @type peak_list: list of pyms.Peak.Class.Peak
    
    @return: A simulated Intensity Matrix object
    @rtype: pyms.GCMS.Class.IntenstityMatrix
    
    @author: Sean O'Callaghan
    """
    n_mz = len(mass_list)
    n_scan = len(time_list)
        
    t1 = time_list[0]
    period = time_list[1] - t1
    
    # initialise a 2D numpy array for intensity matrix
    i_array = numpy.zeros((n_scan, n_mz), 'd') 
    
    
    for peak in peak_list:
        print "-", 
        index = int((peak.get_rt() - t1)/period) 
        height = sum(peak.get_mass_spectrum().mass_spec)
        # standard deviation = area/(height * sqrt(2/pi)) 
        sigma = peak.get_area() / (height * (math.sqrt(2*math.pi)))
        print "width", sigma
        for i in range(len(peak.get_mass_spectrum().mass_list)):
            ion_height = peak.get_mass_spectrum().mass_spec[i]
            ic = chromatogram(n_scan, index, sigma, ion_height)
            i_array[:, i] += ic
            
    im = IntensityMatrix(time_list, mass_list, i_array)
    
    return im
        



def chromatogram(n_scan, x_zero, sigma, peak_scale):

    """
    @summary: Returns a simulated ion chromatogram of a pure component
              The ion chromatogram contains a single gaussian peak.
    
    @param n_scan: the number of scans
    @type n_scan: intType
    
    @param x_zero: The apex of the peak
    @type x_zero: intType
    
    @param sigma: The standard deviation of the distribution
    @type sigma: floatType
    
    @param: peak_scale: the intensity of the peak at the apex
    @type peak_scale: floatType
    
    @return: a list of intensities
    @rtype: listType
    
    @author: Sean O'Callaghan
    """

   
    ic = numpy.zeros((n_scan), 'd')                     
    
    for i in range(n_scan):
	x = float(i)
	
        ic[i] = gaussian(x,x_zero,sigma,peak_scale)

    return ic

def gaussian(point, mean, sigma, scale):
    """
    @summary: calculates a point on a gaussian density function
    
    f = s*exp(-((x-x0)^2)/(2*w^2));
    
    @param point: The point currently being computed
    @type point: floatType
    
    @param mean: The apex of the peak
    @type mean: intType
    
    @param sigma: The standard deviation of the gaussian
    @type sigma: floatType
     
    @param scale: The height of the apex
    @type scale: floatType
    
    @return: a single value from a normal distribution
    @rtype: floatType
    
    @author: Sean O'Callaghan
    """
     
    return scale*math.exp((-(point-mean)**2)/(2*(sigma**2)))

def add_gaussc_noise(im, scale):
    """
    @summary: adds noise to an IntensityMatrix object
     
    @param im: the intensity matrix object
    @param im: pyms.GCMS.Class.IntensityMatrix
    
    @param scale: the scale of the normal distribution from
                  which the noise is drawn
    @type scale: floatType
    
    @author: Sean O'Callaghan
    """
    n_scan, n_mz = im.get_size()
    
    for i in range(n_mz):
        ic = im.get_ic_at_index(i)
        add_gaussc_noise_ic(ic, scale)
        im.set_ic_at_index(i, ic)
    

def add_gaussv_noise(im, scale, cutoff, prop):
    """
    @summary: adds noise to an IntensityMatrix object
     
    @param im: the intensity matrix object
    @param im: pyms.GCMS.Class.IntensityMatrix
    
    @param scale: the scale of the normal distribution from
                  which the noise is drawn
    @type scale: floatType
    
    @param cutoff: The level below which the intensity of the ic at that point
                   has no effect on the scale of the noise distribution
    @type cutoff: intType
    
    @param scale: The scale of the normal distribution for ic values 
    @type scale: intType
    
    @param prop: For intensity values above the cutoff, the scale is 
                 multiplied by the ic value multiplied by prop
    @type prop: floatType
    
    @author: Sean O'Callaghan
    """
    n_scan, n_mz = im.get_size()
    
    for i in range(n_mz):
        ic = im.get_ic_at_index(i)
        add_gaussv_noise_ic(ic, scale, cutoff, prop)
        im.set_ic_at_index(i, ic)
    
    
    
def add_gaussc_noise_ic(ic, scale):
    """
    @summary: adds noise drawn from a normal distribution 
    with constant scale to an ion chromatogram
    
    @param ic: The ion Chromatogram
    @type ic: pyms.GCMS.IonChromatogram
    
    @param scale: The scale of the normal distribution
    @type scale: intType
        
    @author: Sean O'Callaghan
    """
    
    noise = numpy.random.normal(0.0, scale, (len(ic)))
     
    i_array_with_noise = ic.get_intensity_array() + noise
    ic.set_intensity_array(i_array_with_noise)
    
    
def add_gaussv_noise_ic(ic, scale, cutoff, prop):
    """
    @summary: adds noise to an ic. The noise value is drawn from a normal
              distribution, the scale of this distribution depends on the 
              value of the ic at the point where the noise is being added
     
    @param ic: The IonChromatogram
    @type ic: pyms.GCMS.IonChromatogram
    
    @param cutoff: The level below which the intensity of the ic at that point
                   has no effect on the scale of the noise distribution
    @type cutoff: intType
    
    @param scale: The scale of the normal distribution for ic values below the cutoff
                 is modified for values above the cutoff
    @type scale: intType
    
    @param prop: For ic values above the cutoff, the scale is multiplied by the ic
                 value multiplied by prop
    @type prop: floatType
    
    @author: Sean O'Callaghan
    """
    
    noise = numpy.zeros(len(ic))
    
    i_array = ic.get_intensity_array()
    time_list = ic.get_time_list()
    
    for i in range(len(ic)):
        if i_array[i] < cutoff:
            noise[i] = numpy.random.normal(0.0, scale, 1)
        else:
            noise[i] = numpy.random.normal(0.0, scale*i_array[i]*prop, 1)
        
    i_array_with_noise = noise + i_array
    ic.set_intensity_array(i_array_with_noise)
    
    
