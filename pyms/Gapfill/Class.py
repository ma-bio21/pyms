'''
Provides a class for handling Missing Peaks in an output file (i.e. area.csv)
'''
import string

class MissingPeak(object):
    '''
    @summary:    Class to encapsulate a peak object identified as missing in 
        the output area matrix fom PyMS.
    
    @author: Jairus Bowne
    @author: Sean O'Callaghan
    '''
    
    def __init__(self, ci, qual_ion_1, qual_ion_2, rt=0.0):
        '''
        @summary:

        
        @param ci:    Common ion for the peak across samples in an experiment
        @type ci:    IntType
        
        @param UID:    Unique IDentifier for peak, listing top ions, ratio and
            retention time for the peak across an experiment (row in area.csv)
        @type UID:    StringType
        
        # TODO: Determine which of these will be used.
        @param SAMPLE_NAME:    The 
                                   sample name (filename) 
                                   experiment code (group-number-filename) 
                                                                           of
            the originating data file (column in area.csv)
        @type SAMPLE_NAME:    StringType
        
        @param rt:    Retention time of the peak. May or may not be set
        @type rt:    FloatType
        '''
        
        self.__ci = ci
        self.__qual_1 = qual_ion_1
        self.__qual_2 = qual_ion_2
        self.__rt = rt
        self.__exact_rt = 'na'
        self.__ci_area = 'na'
        
    
    def get_ci(self):
        '''
        @summary:    Returns the common ion for the peak object across an 
            experiment
        
        @return:    Common ion for the peak
        @rtype:    IntType
        
        @author:    Jairus Bowne
        '''
        
        return self.__ci
    
    def get_qual_ion1(self):
        '''
        @summary:    Returns the top (most abundant) ion for the peak object
        
        @return:    Most abundant ion
        @rtype:    IntType
        
        @author:    Jairus Bowne
        '''
        
        '''
        # TODO: Consider the abundance of ions when some (i.e. 73, 147) have 
            been im.null_mass()'d. Is there a way to determine whether that
            has been done to generate the original peak list?
        '''
        
        return self.__qual_1
        #return int(string.split(self.__UID, '-')[0])
    
    def get_qual_ion2(self):
        '''
        @summary:    Returns the second most abundant ion for the peak object
        
        @return:    Second most abundant ion
        @rtype:    IntType
        
        @author:    Jairus Bowne
        '''
        
        return self.__qual_2
        #return int(string.split(self.__UID, '-')[1])

    def set_ci_area(self, ci_area):
        """
        @summary: sets the common ion area calculated by the gap fill
                  algorithm
        @param ci_area: The area of the common ion
        @type ci_area: intType

        """
        self.__ci_area = ci_area

    def get_ci_area(self):
        """
        @summary: returns the common ion area

        @return ci_area: The area of the common ion
        @rtype: intType
        """
        return self.__ci_area
        

    def get_rt(self):
        """
        @summary: returns the retention time of the peak

        @return: the retention time of the peak
        @rtype: floatType
        """
        return self.__rt

    def set_exact_rt(self, rt):
        """
        @summary: sets the retention time of a peak

        @param rt: The retention time of the apex of the peak
        @type rt: floatType
        """
        self.__exact_rt = rt

    def get_exact_rt(self):
         """
        @summary: returns the retention time of the peak

        @return: the retention time of the peak
        @rtype: floatType
        """
         return self.__exact_rt
    




class Sample(object):
    """
    @summary: A collection of MissingPeak objects

    @author: Sean O'Callaghan
    """

    def __init__(self, sample_name, matrix_position):
        """
        @summary: A collection of MissingPeak objects

        @param sample_name: the experiment code/name
        @type sample_name: stringType

        @param matrix_position: position along x-axis
                                where sample is located
        @type matrix_position: intType

        """
        self.__sample_name = sample_name
        self.__matrix_position = matrix_position

        self.__missing_peak_list = []

    def get_name(self):
        """
        @summary: Returns the sample name

        @return: The name of the sample
        @rtype: stringType
        """
        return self.__sample_name
    
    def add_missing_peak(self, missing_peak):
        """
        @summary: Add a new MissingPeak object to the Sample

        @param missing_peak: The missing peak object to be added
        @type missing_peak: pyms.GapFilling.Class.MissingPeak
        """
        ###
        # Do some checking here!!!
        ###
        self.__missing_peak_list.append(missing_peak)

    def get_missing_peaks(self):
        """
        @summary: Returns a list of the MissingPeak objects
                  in the Sample object
        @return: list of pyms.GapFilling.Class.MissingPeak
        @rtype: listType
        """
        return self.__missing_peak_list

    def get_mp_rt_area_dict(self):
        """
        @summary: returns a dictionary containing rt:area pairs

        @return: a dict containing rt:area pairs
        @rtype: dictType
        """
        rt_area_dict = {}        
        for peak in self.__missing_peak_list:
            rt = peak.get_rt()
            area = peak.get_ci_area()

            rt_area_dict[rt] = area
        

        return rt_area_dict

    def get_mp_rt_exact_rt_dict(self):
        """
        @summary:returns a dictionary containing average_rt:exact_rt pairs

        @return: a dict of average_rt:exact_rt pairs
        @rtype: dictType
        """

        rt_exact_rt_dict = {}
        for peak in self.__missing_peak_list:
            rt = peak.get_rt()
            exact_rt = peak.get_exact_rt()

            rt_exact_rt_dict[rt] = exact_rt

        return rt_exact_rt_dict
        
        
