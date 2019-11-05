#!/usr/bin/env python
import gzip
import argparse
import numpy as np
import matplotlib.pyplot as plt

def get_args():
    """ get file names for read 1, read 2, index 1,index 2, barcode file, and quality score"""
    parser = argparse.ArgumentParser(description=' File name to be used')
    parser.add_argument('-R', '--Read', type=str,help='File Name for Read 1 as a gzip file', required=True)
    parser.add_argument('-r', '--read', type=str, help='File Name for Read 2 as a gzip file', required=True)
    parser.add_argument('-I', '--Index', type=str, help='File Name for index 1 as a gzip file', required=True)
    parser.add_argument('-i', '--index', type=str, help='File Name for index 2 as a gzip file', required=True)
    parser.add_argument('-b', '--barcodes', type=str, help='Known barcodes file as tsv with code and index on each line', required=True)
    parser.add_argument('-q', '--quality', type=str, help='quality score threshold', required=True)
    return parser.parse_args()
args= get_args()
r1 = args.Read
r2 = args.read
i1 = args.Index
i2 = args.index
b = args.barcodes
q = args.quality

def revcomp(seq):
    """ revcomp creates a reverse compliment of a given string of nucleotides"""
    rev = seq[::-1]
    new = ''
    for n in rev:
        if n == 'A':
            new += 'T'
        elif n == 'T':
            new += 'A'
        elif n == 'C':
            new += 'G'
        elif n == 'G':
            new += 'C'
        elif n == 'N':
            new += 'N'

    return new

def barcodes(b):
    """takes barcodes from a tab separated file and creates an array of known barcodes"""
    with open(b,'r') as bar:
        barcodes = []
        for line in bar:
            x = line.split('\t')[1]
            y = x.strip()
            barcodes.append(y)
        return barcodes
bar = barcodes(b)
def qual(q,seq):
    """takes a quality score and checks if the mean phred score of the score line is greater

     q -- quality score cutoff
     seq -- encoded quality score line
     returns value -- true if mean score is greater and false if not"""
    sum=0
    for score in seq:
        # sums all quality scores
        n = (((ord(score) - 33)))
        sum += n
    if (int(sum)/len(seq))>=int(q):
        return True
    else:
        return False

def header(header,i1,i2):
    """appends i1_i2 to the end of the header line, returns the adjusted header"""
    sub = '_'+i1+'_'+i2
    x = header + sub
    return x
def demult(r1,r2,i1,i2,bar,q):
    """Takes read and index files record by record, checks if indexes contain an N, and if it does places it in low
    quality, checks if average read quality is above or equal to user input threshold, if not places record in low
    quality files,checks if index 1 is the reverse compliment of index two, if not places in unmatched index files,
    if true places in that barcode file. All headers are appended to include the barcodes at the start. Also includes
    counters for all reads.All headers are appended to include the barcodes at the start. Returns counts for each index,
    unknown indexes,and index hopped indices as well as creating files """

    with gzip.open(r1, 'rt') as r1, gzip.open(r2, 'rt') as r2, gzip.open(i1, 'rt') as i1, gzip.open(i2, 'rt') as i2, \
    open('unknown_fw.fq', 'a') as uf, open('unknown_rv.fq', 'a') as ur, open('index_hopped_fw.fq', 'a') as ihf, \
    open('index_hopped_rv.fq', 'a') as ihr:
        # opens all read files and unknown and index hopped files to write to
        codes = []
        for barcode in bar:
            # creates a list containing file names for all fw and rv reads for all barcodes
            codes.append(barcode + '_fw.fq')
            codes.append(barcode + '_rv.fq')
        files = {code: open(code,'a') for code in codes}
        # opens all files contained in codes
        counts = {'unknown': 0, 'index_hopped': 0}
        # Initializes counts dictionary with unknown and index_hopped being set to zer
        swap = {}
        # Dictionary containing every barcode set to zero
        swapped = {}
        # Dictionary of dictionaries where each key is one barcode, and each value is a dictionary containing a second
        # Barcode and a count of how many times we've the first barcode paired with the second
        for index in bar:
            swap[index]=0
        for index in bar:
            swapped[index]={}
            for item in swap:
                swapped[index][item]=0
        # Generates a dictionary of dictionaries containg all possible pairs
        ln = 0
        r1temp = []
        r2temp = []
        i1temp = []
        i2temp = []
        i2_rev = []
        # initializes temp variables to contain all lines of a record
        for r1_line , r2_line , i1_line , i2_line in zip(r1,r2,i1,i2):
            ln += 1
            r1temp.append(r1_line.strip())
            r2temp.append(r2_line.strip())
            i1temp.append(i1_line.strip())
            i2temp.append(i2_line.strip())
            # reads one line from each of the four files and adds it to the temp variable
            if ln % 4 == 2:
                i2_rev = revcomp(i2_line.strip())
                # creates a reverse compliment of the i2 barcode if its reading the sequence line
            if ln % 4 == 0:
                # checks to see if there should a full record in each of the temp variables
                if qual(q,r1temp[3]) and qual(q,r2temp[3]) and qual(q,i1temp[3]) and qual(q,i2temp[3]) is True:
                    # checks if the all quality score lines pass user defined quality score cutoffs
                    if str(i1temp[1]) == str(i2_rev) and str(i1temp[1]) in bar:
                        # Checks if the barcodes are equal and in known barcodes
                        headerf = header(r1temp[0], i1temp[1], i2temp[1])
                        headerr = header(r1temp[0], i1temp[1], i2temp[1])
                        # creates new headers containing barcode information
                        files[(i1temp[1] + "_fw.fq")].write(headerf + '\n' + r1temp[1] + '\n' + r1temp[2] + '\n' + r1temp[3] + '\n')
                        files[(i1temp[1] + "_rv.fq")].write(headerr + '\n' + r2temp[1] + '\n' + r2temp[2] + '\n' + r2temp[3] + '\n')
                        # writes the records to the appropriate read files
                        if i1temp[1] in counts:
                            counts[i1temp[1]]+=1
                            # checks if the barcode is in counts, and increments by one if it is
                            r1temp = []
                            r2temp = []
                            i1temp = []
                            i2temp = []
                            # reset temp variables
                        else:
                            counts[i1temp[1]] = 1
                            # adds barcode to counts if not present and sets it equal to 1
                            r1temp = []
                            r2temp = []
                            i1temp = []
                            i2temp = []
                            #resets temp variables
                    elif str(i1temp[1]) != str(i2_rev) and str(i1temp[1]) in bar and str(i2_rev) in bar:
                        # checks if the barcodes aren't the same but are both known barcodes
                        headerf = header(r1temp[0], i1temp[1], i2temp[1])
                        headerr = header(r1temp[0], i1temp[1], i2temp[1])
                        # generates new headers
                        ihf.write(headerf + '\n' + r1temp[1] + '\n' + r1temp[2] + '\n' + r1temp[3] + '\n')
                        ihr.write(headerr + '\n' + r2temp[1] + '\n' + r2temp[2] + '\n' + r2temp[3] + '\n')
                        # writes the records to the index hopped files
                        counts['index_hopped'] += 1
                        swapped[i1temp[1]][i2_rev] += 1
                        # adds to counts of index hopped and the the count of the specific index pair in swapped dict
                        r1temp = []
                        r2temp = []
                        i1temp = []
                        i2temp = []
                        #resets temp variables
                    else:
                        headerf = header(r1temp[0], i1temp[1], i2temp[1])
                        headerr = header(r1temp[0], i1temp[1], i2temp[1])
                        # generates new headers
                        uf.write(headerf + '\n' + r1temp[1] + '\n' + r1temp[2] + '\n' + r1temp[3] + '\n')
                        ur.write(headerr + '\n' + r2temp[1] + '\n' + r2temp[2] + '\n' + r2temp[3] + '\n')
                        #writes to unknown file
                        counts['unknown'] += 1
                        # increments unknown counter by 1
                        r1temp = []
                        r2temp = []
                        i1temp = []
                        i2temp = []
                        # resets temp variables
                else:
                    headerf = header(r1temp[0],i1temp[1],i2temp[1])
                    headerr = header(r1temp[0],i1temp[1],i2temp[1])
                    # generates new headers
                    uf.write(headerf + '\n' + r1temp[1] + '\n' + r1temp[2] + '\n' + r1temp[3] + '\n')
                    ur.write(headerr + '\n' + r2temp[1] + '\n' + r2temp[2] + '\n' + r2temp[3] + '\n')
                    # writes to unknown files
                    counts['unknown'] += 1
                    # increments unknown counter
                    r1temp = []
                    r2temp = []
                    i1temp = []
                    i2temp = []
                    # resets temp variables
    for file in files.values():
        file.close()
    # close all files not opened using with open
    for item in counts:
        if item != 'unknown' and item != 'index_hopped':
            swapped[item][item] = counts[item]
    # adds the correctly matched barcodes to swapped dictionary for heatmap generation
    print("Barcode\tCount")
    for item in counts:
        print(str(item)+'\t'+str(counts[item]))
    # prints counts of each barcode, unknown reads, and index hopped reads
    return counts, swapped


counts, swapped = demult(r1, r2, i1, i2, bar, q)

def countplot(counts):
    """Uses matplotlib to generate a basic bar graph of counts for each barcode as a png

    Counts -- a dictionary object where barcodes are keys and counts are values
    Return value: counts.png """
    y = []
    x = []
    for item in counts:
        x.append(item)
        y.append(counts[item])
    plt.bar(x, y, align='center', alpha=0.5)
    plt.ylabel('Count')
    plt.xlabel('Barcode')
    plt.xticks(np.arange(26), x, rotation=90)
    plt.title('Counts for each barcode')
    plt.savefig('counts.png', bbox_inches = "tight")
    plt.close()
    return
def swapplot(swapped):
    """Generates a heatmap of barcode combinations using matplotlib

    swapped -- a dictionary of dictionaries containing counts of index hopping
    return value -- a heatmap as a .png """
    z = []
    twod = []
    for d in swapped:
        for k in swapped[d]:
            if len(z) < 24:
                z.append(swapped[d][k])
            else:
                twod.append(z)
                z = [swapped[d][k]]
    # generates an array of barcode values
    bars = barcodes(b)
    x = range(24)
    y = range(24)
    intensity = np.array(twod)
    x,y = np.meshgrid(x,y)
    plt.pcolormesh(x, y, intensity)
    plt.colorbar()
    plt.xticks(np.arange(24),bars,rotation=90)
    plt.yticks(np.arange(24),bars)
    plt.title('Counts for each barcode combo')
    plt.gcf().subplots_adjust(left=0.15)
    plt.savefig('heatmap.png', bbox_inches="tight")
    plt.close()
    return



countplot(counts)
swapplot(swapped)







