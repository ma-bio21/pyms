"""
Provides a class to model signal peak
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

import math
import copy

from pyms.GCMS.Class import IntensityMatrix, MassSpectrum
from pyms.Utils.Error import error
from pyms.Utils.Utils import is_int, is_number, is_list, is_boolean, is_str
from pyms.Utils.IO import open_for_writing, close_for_writing

class Peak:

    """
    @summary: Models a signal peak

    A signal peak object
    A peak object is initialised with retention time and
    Either an ion mass, a mass spectrum or None

    @author: Vladimir Likic
    @author: Andrew Isaac
    """

    def __init__(self, rt=0.0, ms=None, minutes=False):

        """
        @param rt: Retention time
        @type rt: FloatType
        @param ms: A ion mass, or spectra of maximising ions
        @type ms: FloatType, pyms.GCSM.Class.MassSpectrum
        @param minutes: Retention time units flag. If True, retention time
            is in minutes; if False retention time is in seconds
        @type minutes: BooleanType
        """

        if not is_number(rt):
            error("'rt' must be a number")

        if not ms == None and \
            not isinstance(ms, MassSpectrum) and \
            not is_number(ms):
            error("'ms' must be a Float or a MassSpectrum object")

        if minutes:
            rt = rt*60.0

        self.__minutes = minutes

        # basic peak attributes
        self.__rt = float(rt)
        # these two attributes are required for
        # setting the peak mass spectrum
        if not ms == None:
            if isinstance(ms, MassSpectrum):
                # mass spectrum
                self.__mass_spectrum = ms
                self.__ic_mass = None
                self.make_UID()

                # TEST: to test if this speeds things up
                self.mass_spec = ms.mass_spec
                self.ms = ms

            else:
                # single ion chromatogram properties
                self.__ic_mass = ms
                self.__mass_spectrum = None
                self.make_UID()

                # TEST: to test if this speeds things up
                self.mass_spec = None


        self.__pt_bounds = None
        self.__area = None
        self.__ion_areas = {}

        # TEST: to test if this speeds things up
        self.rt = self.__rt

    def make_UID(self):

        """
        @summary: Create a unique peak ID (UID) based on:
            Int masses of top two intensities, their ratio and RT (as
            Mass1-Mass2-Ratio*100-RT); or the single mass as int and RT.

        @author: Andrew Isaac
        """

        minutes = 1.0
        if self.__minutes == True:
            minutes = 60.0
        if self.__mass_spectrum != None:
            mass_list = self.__mass_spectrum.mass_list
            mass_spec = self.__mass_spectrum.mass_spec
            # find top two masses
            best = 0
            best_ii = 0
            best2_ii = 0
            for ii in range(len(mass_spec)):
                if mass_spec[ii] > best:
                    best = mass_spec[ii]
                    best2_ii = best_ii
                    best_ii = ii
            ratio = int(100*mass_spec[best2_ii]/best)
            UID = "%d-%d-%d-%.2f" % (int(mass_list[best_ii]),
                    int(mass_list[best2_ii]), ratio, self.__rt/minutes)
        elif self.__ic_mass != None:
            UID = "%d-%.2f" % (int(self.__ic_mass), self.__rt/minutes)
        else:
            UID =  "%.2f" % self.__rt/minutes

        self.__UID = UID

    def get_UID(self):

        """
        @summary: Return the unique peak ID (UID) based on:
            Int masses of top two intensities and their ratio (as
            Mass1-Mass2-Ratio*100); or the single mass as int.

        @return: UID string
        @rtype: StringType

        @author: Andrew Isaac
        """

        return self.__UID
    
    def get_third_highest_mz(self):

        """
        @summary: returns the mz value with the third highest
                  intensity

        @return: the ion with the third highest intensity
        @rtype: intType

        """

        if self.__mass_spectrum != None:
            mass_list = self.__mass_spectrum.mass_list
            mass_spec = self.__mass_spectrum.mass_spec
            # find top two masses
            best = 0
            best_ii = 0
            best2_ii = 0
            best3_ii = 0
            for ii in range(len(mass_spec)):
                if mass_spec[ii] > best:
                    best = mass_spec[ii]
                    best3_ii = best2_ii
                    best2_ii = best_ii
                    best_ii = ii

        return int(mass_list[best3_ii])
        

    def set_pt_bounds(self, pt_bounds):

        """
        @summary: Sets peak boundaries in points

        @param pt_bounds: A list containing left, apex, and right
            peak boundaries in points, left and right are offsets
        @type pt_bounds: ListType

        @return: none
        @rtype: NoneType
        """

        if not is_list(pt_bounds):
            error("'pt_bounds' must be a list")

        if not len(pt_bounds) == 3:
            error("'pt_bounds' must have exactly 3 elements")
        else:
            for item in pt_bounds:
                if not is_int(item):
                    error("'pt_bounds' element not an integer")

        self.__pt_bounds = pt_bounds

    def get_pt_bounds(self):

        """
        @summary: Gets peak boundaries in points

        @return: A list containing left, apex, and right
            peak boundaries in points, left and right are offsets
        @rtype: ListType

        @author: Andrew Isaac
        """

        return copy.deepcopy(self.__pt_bounds)

    def get_area(self):

        """
        @summary: Gets the area under the peak

        @return: The peak area
        @rtype: FloatType

        @author: Andrew Isaac
        """

        return self.__area

    def set_area(self, area):

        """
        @summary: Sets the area under the peak

        @param area: The peak area
        @type area: FloatType

        @author: Andrew Isaac
        """

        if not is_number(area) or area <= 0:
            error("'area' must be a positive number")

        self.__area = area
        
    def get_ion_area(self, ion):
        """
        @summary: gets the area of a single ion chromatogram
                  under the peak

        @param ion: the ion who's IC is to be returned
        @type ion: IntType

        @return: The area of the ion under this peak
        @rtype: FloatType
        """
        try:
            return self.__ion_areas[ion]
        except KeyError:
            return None
            
    def get_ion_areas(self):
        """
        @summary: returns a copy of the ion areas dict

        @return: The dictionary of ion:ion area pairs
        @rtype: dictType
        """
        if len(self.__ion_areas) == 0:
            error("no ion areas set")

        return copy.deepcopy(self.__ion_areas)

    
    def set_ion_area(self, ion, area):
        """
        @summary: sets the area for a single ion

        @param ion: the ion who's area is being entered
        @param area: the area under the IC of ion

        @author: Sean O'Callaghan
        """

        self.__ion_areas[ion] = area

    def set_ion_areas(self, ion_areas):
        """
        @summary: set the ion:ion area pair dictionary

        @param ion_areas: the dictionary of ion:ion_area pairs
        @type ion_areas: dictType
        """

        self.__ion_areas = ion_areas
    
    def get_int_of_ion(self, ion):
        """
        @summary: returns the intensity of a given ion

        @param ion: the m/z value of the ion of interest
        @type ion: intType

        @return: The intensity of the given ion in this peak
        @rtype: intType
        """

        for i,mass in enumerate(self.__mass_spectrum.mass_list):
            if mass == ion:
                index = i

        try:
            intensity = self.__mass_spectrum.mass_spec[index]
        except:
            intensity = None

        return intensity  

 

    def set_ic_mass(self, mz):

        """
        @summary: Sets the mass for a single ion chromatogram peak
            Clears the mass spectrum

        @param mz: The mass of the ion chromatogram that the peak is from
        @type mz: FloatType

        @return: none
        @rtype: NoneType
        """

        if not is_number(mz):
            error("'mz' must be a number")
        self.__ic_mass = mz
        # clear mass spectrum
        self.__mass_spectrum = None
        self.make_UID()

        # TEST: to test if this speeds things up
        self.mass_spec = None

    def set_mass_spectrum(self, ms):

        """
        @summary: Sets the mass spectrum
            Clears the mass for a single ion chromatogram peak

        @param ms: The mass spectrum at the apex of the peak
        @type ms: pyms.GCSM.Class.MassSpectrum

        @return: none
        @rtype: NoneType
        """

        if not isinstance(ms, MassSpectrum):
            error("'ms' must be a MassSpectrum object")

        self.__mass_spectrum = ms
        # clear ion mass
        self.__ic_mass = None
        self.make_UID()

        # TEST: to test if this speeds things up
        self.mass_spec = ms.mass_spec

    def get_rt(self):

        """
        @summary: Return the retention time

        @return: Retention time
        @rtype: FloatType
        """

        return self.__rt

    def get_ic_mass(self):

        """
        @summary: Gets the mass for a single ion chromatogram peak

        @return: The mass of the single ion chromatogram that the peak is from
        @rtype: FloatType or IntType
        """

        return self.__ic_mass

    def get_mass_spectrum(self):

        """
        @summary: Gets the mass spectrum at the apex of the peak

        @return: The mass spectrum at the apex of the peak
        @rtype: pyms.GCSM.Class.MassSpectrum
        """

        return copy.deepcopy(self.__mass_spectrum)

## TODO: What is this?
    def find_mass_spectrum(self, data, from_bounds=False):

        """
        @summary: Sets peak mass spectrum from the data
            Clears the single ion chromatogram mass

        @param data: An IntensityMatrix object
        @type data: pyms.GCMS.IntensityMatrix
        @param from_bounds: Indicator whether to use the attribute
            'pt_bounds' or to find the peak apex from the peak
            retention time
        @type from_bounds: BooleanType

        @return: none
        @rtype: NoneType
        """

        if not isinstance(data, IntensityMatrix):
            error("'data' must be an IntensityMatrix")

        if not is_boolean(from_bounds):
            error("'from_bounds' must be boolean")

        if from_bounds:
            if self.__pt_bounds == None:
                error("pt_bounds not set for this peak")
            else:
                pt_apex = self.__pt_bounds[1]
        else:
            # get the index of peak apex from peak retention time
            pt_apex = data.get_index_at_time(self.__rt)

        # set the mass spectrum
        self.__mass_spectrum = data.get_ms_at_index(pt_apex)

        # clear single ion chromatogram mass
        self.__ic_mass = None
        self.make_UID()

        # TEST: to test if this speeds things up
        self.mass_spec = self.__mass_spectrum.mass_spec

    def crop_mass(self, mass_min, mass_max):

        """
        @summary: Crops mass spectrum

        @param mass_min: Minimum mass value
        @type mass_min: IntType or FloatType
        @param mass_max: Maximum mass value
        @type mass_max: IntType or FloatType

        @return: none
        @rtype: NoneType

        @author: Andrew Isaac
        """

        if not is_number(mass_min) or not is_number(mass_max):
            error("'mass_min' and 'mass_max' must be numbers")
        if mass_min >= mass_max:
            error("'mass_min' must be less than 'mass_max'")

        mass_list = self.__mass_spectrum.mass_list

        if mass_min < min(mass_list):
            error("'mass_min' is less than the smallest mass: %d" \
              % min(mass_list))
        if mass_max > max(mass_list):
            error("'mass_max' is greater than the largest mass: %d" \
              % max(mass_list))

        # pre build mass_list and list of indecies
        new_mass_list = []
        new_mass_spec = []
        mass_spec = self.__mass_spectrum.mass_spec
        for ii in range(len(mass_list)):
             mass = mass_list[ii]
             if mass >= mass_min and mass <= mass_max:
                 new_mass_list.append(mass)
                 new_mass_spec.append(mass_spec[ii])

        self.__mass_spectrum.mass_list = new_mass_list
        self.__mass_spectrum.mass_spec = new_mass_spec

        if len(new_mass_list) == 0:
            error("mass spectrum is now empty")
        elif len(new_mass_list) < 10:
            print " WARNING: peak mass spectrum contains < 10 points"

        # update UID
        self.make_UID()

        # TEST: to test if this speeds things up
        self.mass_spec = self.__mass_spectrum.mass_spec

    def null_mass(self, mass):

        """
        @summary: Ignore given mass in spectra

        @param mass: Mass value to remove
        @type mass: IntType or FloatType

        @author: Andrew Isaac
        """

        if self.__mass_spectrum == None:
            error("mass spectrum not set for this peak")

        if not is_number(mass):
            error("'mass' must be numbers")

        mass_list = self.__mass_spectrum.mass_list

        if mass < min(mass_list) or mass > max(mass_list):
            error("'mass' not in mass range:", min(mass_list), "to", \
                max(mass_list))

        best = max(mass_list)
        ix = 0
        for ii in range(len(mass_list)):
            tmp = abs(mass_list[ii] - mass)
            if tmp < best:
                best = tmp
                ix = ii

        self.__mass_spectrum.mass_spec[ix] = 0

        # update UID
        self.make_UID()

        # TEST: to test if this speeds things up
        self.mass_spec = self.__mass_spectrum.mass_spec
