import matplotlib.pyplot as plt 
from Bio import SeqIO
from collections import defaultdict
import collections
import numpy as np
import math
from collections import OrderedDict 
from Bio.SeqIO.QualityIO import FastqGeneralIterator
from collections import Counter


ENCODINGS = {
            "Sanger Phred+33": (33, 0, 40),  # raw reads typically (0, 40)
            "Solexa Solexa+64": (64, -5, 40),  # raw reads typically (-5, 40)
            # raw reads typically (0, 40)
            "Illumina 1.3+ Phred+64": (64, 0, 40),
            # raw reads typically (3, 41) - 0,1 unused and 2 - 'B' - quality control indicator
            "Illumina 1.5+ Phred+64": (64, 2, 41),
            # raw reads typically (0, 41)
            "Illumina 1.8+ Phred+33": (33, 0, 41)
        }


def findQuality(filename):
 fname=filename
 max_value=-9999
 min_value=9999
 with open(filename) as handle:
    for (title, sequence, quality) in FastqGeneralIterator(handle):
        ascii_score=[ord(number) for number in quality]
        if min(ascii_score)<min_value:
            min_value=min(ascii_score)
        if max(ascii_score)>max_value:
            max_value=max(ascii_score)
    return(min_value,max_value)


def identify_encoding(filename):
        ENCODINGS = {
            "Sanger Phred+33": (33, 0, 40),
            "Solexa Solexa+64": (64, -5, 40),  
            "Illumina 1.3+ Phred+64": (64, 0, 40),
            "Illumina 1.5+ Phred+64": (64, 2, 41),
            "Illumina 1.8+ Phred+33": (33, 0, 41)
        }
        fname=filename

        min_and_max = findQuality(fname)
        min_val = min_and_max[0]
        max_val = min_and_max[1]
        diff_sums_dict = {}

        for encoding in ENCODINGS:
            diff_from_expected_min = abs(
                min_val - ENCODINGS[encoding][0] - ENCODINGS[encoding][1])
            diff_from_expected_max = abs(
                max_val - ENCODINGS[encoding][0] - ENCODINGS[encoding][2])
            diff_sum = diff_from_expected_min + diff_from_expected_max
            diff_sums_dict[encoding] = diff_sum

        smallest_diff_value_encoding = min(
            diff_sums_dict, key=diff_sums_dict.get)

        print("Found encoding: %s" % smallest_diff_value_encoding)

def find_percentage(filename):
    f=open("results.txt","w")    
    fname=filename
    resultlist=[]
    seqcounterlist=[]
    seqlist=[]
    seqcounter=1
    with open(fname) as handle:
     for (title, sequence, quality) in FastqGeneralIterator(handle):
       # print(sequence)
        counter=0
        for char in sequence:
            if(char=="C" or char=="G"):
             counter=counter+1
        result=round(counter/len(sequence),3)
       



        f.write(str(seqcounter)+" ")
        seqcounterlist.append(seqcounter)
        seqcounter=seqcounter+1
        resultlist.append(result)
        seqlist.append(sequence)
        f.write(str(result))
        f.write("\n")
        f.write(str(sequence))
        f.write("\n")

     # x axis values 
# corresponding y axis values
    
   
    plt.hist(resultlist,bins = 75)
    plt.xticks(np.arange(0,1,0.05))
# plotting the points
  
# naming the x axis 
    plt.xlabel('G/C dažnis') 
# naming the y axis 
    plt.ylabel('seku skaičius') 
  
# giving a title to my graph 
    plt.title('Grafikas') 
  
# function to show the plot 
    plt.show() 
       
            


identify_encoding("reads_for_analysis.fastq")
find_percentage("reads_for_analysis.fastq")
