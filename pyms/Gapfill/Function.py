'''
@summary:    # Functions to fill missing peak objects

@author:    Jairus Bowne
@author:    Sean O'Callaghan
'''

import csv
import string
import sys, os, errno, string, numpy
sys.path.append("/x/PyMS")

from pyms.GCMS.IO.ANDI.Function import ANDI_reader
from pyms.GCMS.IO.MZML.Function import mzML_reader
from pyms.GCMS.Function import build_intensity_matrix_i
from pyms.Noise.SavitzkyGolay import savitzky_golay
from pyms.Baseline.TopHat import tophat
from pyms.Deconvolution.BillerBiemann.Function import get_maxima_list_reduced
from pyms.Peak.Function import ion_area

from Class import MissingPeak, Sample


# .csv reader (cloned from gcqc project)
def file2matrix(filename):
    '''
    @summary:    Convert a .csv file to a matrix (list of lists)
    
    @param filename:    Filename (.csv) to convert (area.csv, area_ci.csv)
    @type filename:    StringType
    
    @return:    Data matrix
    @rtype:    ListType (List of lists)
    '''
    
    # open(filename, 'rb')? Or unnecessary?
    with open(filename) as fp:
        reader = csv.reader(fp, delimiter=",",quotechar="\"")
        matrix = []
        for row in reader:
            newrow = []
            for each in row:
                try:
                    each = float(each)
                except:
                    pass
                newrow.append(each)
            matrix.append(newrow)
    
    return matrix


def mp_finder(inputmatrix):
    """
    @summary: setup sample objects with missing peak objects
              Finds the 'NA's in the transformed area_ci.csv file and
              makes Sample objects with them

    @param inputmatrix: Data matrix derived from the area_ci.csv file
    @type inputmatrix: listType

    @return: list of Sample objects
    @rtype: list of pyms.MissingPeak.Class.Sample

    """

    sample_list = []

    try:
        ci_pos = inputmatrix[0].index(' "Quant Ion"')
    except ValueError:
        ci_pos = inputmatrix[0].index('"Quant Ion"')
    
    uid_pos = inputmatrix[0].index('UID')
   

    # Set up the sample objects
    # All entries on line 1 beyond the Qual Ion position are sample names
    for i, sample_name in enumerate(inputmatrix[0][ci_pos:]):

        print sample_name
        sample = Sample(sample_name, i+3) #add 4 to allow for UID, RT,QualIon
        sample_list.append(sample)

    for line in inputmatrix[1:]:
        uid = line[uid_pos]
        common_ion = line[ci_pos]

        qual_ion_1 = uid.split("-")[0]
        qual_ion_2 = uid.split("-")[1]
        rt = uid.split("-")[-1]
        #print rt
        
        for i, area in enumerate(line[ci_pos:]):
            if area == 'NA':
                missing_peak = MissingPeak(common_ion, qual_ion_1, \
                                               qual_ion_2, rt)
                sample_list[i].add_missing_peak(missing_peak)

    return sample_list


def missing_peak_finder(sample, filename, points=13, null_ions=[73, 147],\
                            crop_ions=[50,540], threshold=1000, rt_window=1, filetype='mzml'):
    """
    @summary: Integrates raw data around missing peak locations
              to fill in NAs in the data matrix

    @param  sample: The sample object containing missing peaks
    @type sample: pyms.MissingPeak.Class.Sample

    @param  andi_file: Name of the raw data file
    @type andi_file: stringType

    @param  points: Peak finding - Peak if maxima over 'points' \
                    number of scans (Default 3) 
    @type points: intType

    @param  null_ions: Ions to be deleted in the matrix
    @type null_ions: listType

    @param crop_ions: Range of Ions to be considered
    @type crop_ions: listType 

    @param threshold: Minimum intensity of IonChromatogram allowable to fill\
                      missing peak
    @type threshold: intType

    @param  rt_window: Window in seconds around average RT to look for \
                       missing peak
    @type rt_window: floatType

    @param filetype: either mzml or netcdf
    @type filetype: stringType

    @author: Sean O'Callaghan
    """

    ### some error checks on null and crop ions

    ### a for root,files,dirs in os.path.walk(): loop
    print "Sample:", sample.get_name(), "File:", filename
    
    if filetype.lower() == 'cdf':
        data = ANDI_reader(filename)
    elif filetype.lower() == 'mzml':
        data = mzML_reader(filename)
    else:
        print "file type not valid"
    

    # build integer intensity matrix
    im = build_intensity_matrix_i(data)

    for null_ion in null_ions:
        im.null_mass(null_ion)

    im.crop_mass(crop_ions[0], crop_ions[1])

    # get the size of the intensity matrix
    n_scan, n_mz = im.get_size()

    # smooth data
    for ii in range(n_mz):
        ic = im.get_ic_at_index(ii)
        ic1 = savitzky_golay(ic, points)
        ic_smooth = savitzky_golay(ic1, points)
        ic_base = tophat(ic_smooth, struct="1.5m")
        im.set_ic_at_index(ii, ic_base)

    for mp in sample.get_missing_peaks():

        mp_rt = mp.get_rt()
        common_ion = mp.get_ci()
        qual_ion_1 = float(mp.get_qual_ion1())
        qual_ion_2 = float(mp.get_qual_ion2())
        

        ci_ion_chrom = im.get_ic_at_mass(common_ion)
        print "ci = ",common_ion
        qi1_ion_chrom = im.get_ic_at_mass(qual_ion_1)
        print "qi1 = ", qual_ion_1
        qi2_ion_chrom = im.get_ic_at_mass(qual_ion_2)
        print "qi2 = ", qual_ion_2
        ######
        # Integrate the CI around that particular RT
        #######

        #Convert time to points
        # How long between scans?
        
        points_1 = ci_ion_chrom.get_index_at_time(float(mp_rt))
        points_2 = ci_ion_chrom.get_index_at_time(float(mp_rt)-rt_window)
        print "rt_window = ", points_1 - points_2

        rt_window_points = points_1 - points_2

        maxima_list = get_maxima_list_reduced(ci_ion_chrom, mp_rt, \
                                                  rt_window_points)

        large_peaks = []

        for rt, intens in maxima_list:
            if intens > threshold:
                q1_index = qi1_ion_chrom.get_index_at_time(rt)
                q2_index = qi2_ion_chrom.get_index_at_time(rt)

                q1_intensity = qi1_ion_chrom.get_intensity_at_index(q1_index)
                q2_intensity = qi2_ion_chrom.get_intensity_at_index(q2_index)

                if q1_intensity > threshold/2 and q2_intensity > threshold/2:
                    large_peaks.append([rt, intens])
                
        print('found %d peaks above threshold'%len(large_peaks))

        areas = []
        for peak in large_peaks:
            apex = ci_ion_chrom.get_index_at_time(peak[0])
            ia = ci_ion_chrom.get_intensity_array().tolist()
            area, left, right, l_share, r_share = ion_area(ia, apex, 0)
            areas.append(area)
        ########################
        areas.sort()
        if len(areas)>0:
            biggest_area = areas[-1]
            mp.set_ci_area(biggest_area)
            mp.set_exact_rt("{:.3f}".format(float(mp_rt)/60.0))
            print "found area:", biggest_area, "at rt:", mp_rt
        else:
            print "Missing peak at rt = ", mp_rt
            mp.set_ci_area('na')

def transposed(lists):
   """
   @summary: transposes a list of lists

   @param lists: the list of lists to be transposed
   @type lists: listType
   """
   
   if not lists: return []
   return map(lambda *row: list(row), *lists)


def write_filled_csv(sample_list, area_file, filled_area_file):
    """
    @summary: creates a new area_ci.csv file, replacing NAs with
              values from the sample_list objects where possible
    @param sample_list: A list of sample objects
    @type sample_list: list of Class.Sample

    @param area_file: the file 'area_ci.csv' from PyMS output
    @type area_file: stringType

    @param filled_area_file: the new output file which has NA
                             values replaced
    @type filled_area_file: stringType
    """

    old_matrix = file2matrix(area_file)
    

    #Invert it to be a little more efficent
    invert_old_matrix = zip(*old_matrix)
    #print invert_old_matrix[0:5]

    uid_list = invert_old_matrix[0][1:]
    rt_list = []
    for uid in uid_list:
        rt = uid.split('-')[-1]
        rt_list.append(rt)
        
    
    #print rt_list

    #start setting up the output file
    invert_new_matrix = []
    for line in invert_old_matrix[0:2]:
        invert_new_matrix.append(line)
    
    for line in invert_old_matrix[3:]:
        sample_name = line[0]
        
        new_line = []
        new_line.append(sample_name)
        for sample in sample_list:
            if sample_name in sample.get_name():
                rt_area_dict = sample.get_mp_rt_area_dict()
                #print rt_area_dict

        for i, part in enumerate(line[1:]):
            #print part
            if part == 'NA':
                try:
                    area = rt_area_dict[str(rt_list[i])]
                    new_line.append(area)
                except(KeyError):
                    pass
                    #print 'missing peak not found for rt =', rt_list[i], \
                    #    "in sample:", sample_name
                
            else:
                new_line.append(part)
                
        invert_new_matrix.append(new_line)
    #print invert_new_matrix

    #print len(invert_new_matrix[0]), len(invert_new_matrix)
    

    fp_new = open(filled_area_file, 'w')
 
            
    #    new_matrix = numpy.empty(matrix_size)
    new_matrix = transposed(invert_new_matrix)

    for i, line in enumerate(new_matrix):
        for j, part in enumerate(line):
            fp_new.write(str(part) +',')
        fp_new.write("\n")
            
    fp_new.close()

def write_filled_rt_csv(sample_list, rt_file, filled_rt_file):
    """
    @summary: creates a new rt.csv file, replacing NAs with
              values from the sample_list objects where possible
    @param sample_list: A list of sample objects
    @type sample_list: list of Class.Sample

    @param area_file: the file 'rt.csv' from PyMS output
    @type area_file: stringType

    @param filled_area_file: the new output file which has NA
                             values replaced
    @type filled_area_file: stringType
    """

    old_matrix = file2matrix(rt_file)
    

    #Invert it to be a little more efficent
    invert_old_matrix = zip(*old_matrix)
    #print invert_old_matrix[0:5][0:5]

    uid_list = invert_old_matrix[0][1:]
    rt_list = []
    for uid in uid_list:
        rt = uid.split('-')[-1]
        rt_list.append(rt)
        
    
    #print rt_list

    #start setting up the output file
    invert_new_matrix = []
    for line in invert_old_matrix[0:1]:
        invert_new_matrix.append(line)
    
    for line in invert_old_matrix[2:]:
        sample_name = line[0]
        #print "Sample Name:", sample_name
        
        new_line = []
        new_line.append(sample_name)
        for sample in sample_list:
            #print "sample name:", sample.get_name()
            if sample_name in sample.get_name():
                
                rt_exact_rt_dict = sample.get_mp_rt_exact_rt_dict()
                #print "Got it"
                #print rt_area_dict

        for i, part in enumerate(line[1:]):
            #print part
            if part == 'NA':
                try:
                    rt_new = rt_exact_rt_dict[str(rt_list[i])]
                    new_line.append(rt_new)
                except(KeyError):
                    pass
                    #print 'missing peak not found for rt =', rt_list[i], \
                    #    "in sample:", sample_name
                
            else:
                new_line.append(part)
                
        invert_new_matrix.append(new_line)
    #print invert_new_matrix

    #print len(invert_new_matrix[0]), len(invert_new_matrix)
    

    fp_new = open(filled_rt_file, 'w')
 
            
    #    new_matrix = numpy.empty(matrix_size)
    new_matrix = transposed(invert_new_matrix)

    for i, line in enumerate(new_matrix):
        for j, part in enumerate(line):
            fp_new.write(str(part) +',')
        fp_new.write("\n")
            
    fp_new.close()
        
                
                


    
    
             
