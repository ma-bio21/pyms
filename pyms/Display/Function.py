"""Display.Function.py
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


import matplotlib.pyplot as plt

import sys
sys.path.append('/x/PyMS/')

from pyms.Utils.Error import error
from pyms.GCMS.Class import IonChromatogram 


def plot_ic(ic, line_label=" ", plot_title=" "):
    """
    @summary: Plots an Ion Chromatogram or List of same
    
    @param ic: The ion chromatogram
    @type ic: pyms.GCMS.Class.IonChromatogram
    
    @param line_label: plot legend
    @type line_label: stringType
    
    @param plot_title: A label for the plot
    @type plot_title: String Type
    
    @author: Sean O'Callaghan
    """
    			
    #Plotting Variables
    fig = plt.figure()
    ax = fig.add_subplot(111)
    
   
    if not isinstance(ic, IonChromatogram):
	error("ics argument must be an IonChromatogram\
            or a list of Ion Chromatograms")

    time_list = ic.get_time_list()
			
	
    
    intensity_list = ic.get_intensity_array()
    	
    ic_plot = plt.plot(time_list, intensity_list, label=line_label)
        
    t = ax.set_title(plot_title)
    l = ax.legend()
			
    fig.canvas.draw
    plt.show()
    
    
    
def plot_ms(mass_spec, plot_title=" "):
        
    """ 
    @summary: Plots the mass spec given a list of masses and intensities
        
    @param mass_spec: The mass spectrum at a given time/index
    @type mass_spec: GCMS.Class.MassSpectrum
        
    @param plot_title: A label for the plot
    @type plot_title: String Type
    
    @author: Sean O'Callaghan
    """
        
    fig = plt.figure()
    ax = fig.add_subplot(111)
        
    mass_list = mass_spec.mass_list
    intensity_list = mass_spec.mass_spec
        
    # to set x axis range find minimum and maximum m/z channels
    max_mz = mass_list[0]
    min_mz = mass_list[0]
        
    for i in range(len(mass_list)):
        if mass_list[i] > max_mz:
            max_mz = mass_list[i]
                
    for i in range(len(mass_list)):
        if mass_list[i] < min_mz:
            min_mz = mass_list[i]
        
    mass_spec_plot = plt.bar(mass_list, intensity_list,\
        width=0.01)
        
    x_axis_range = plt.xlim(min_mz, max_mz)
    t = ax.set_title(plot_title)
    plt.show()
