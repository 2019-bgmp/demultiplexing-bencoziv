#!/usr/bin/env python
import gzip
import argparse
import numpy as np
import matplotlib.pyplot as plt

def get_args():
    # get file names for read 1, read 2, index 1, and index 2
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
    # revcomp creates a reverse compliment of a given string of nucleotides
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
    #takes barcodes from a tab separated file and creates an array of known barcodes
    with open(b,'r') as bar:
        barcodes = []
        for line in bar:
            x = line.split('\t')[1]
            y = x.strip()
            barcodes.append(y)
        return barcodes
bar = barcodes(b)
def qual(q,seq):
    #takes a quality score and checks if the mean phred score of a given score line is greater and returns true if it is
    #and false if not
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
    #appends i1_i2 to the end of the header line, returns the adjusted header
    sub = '_'+i1+'_'+i2
    x = header + sub
    return x
def demult(r1,r2,i1,i2,bar,q):
    #takes read and index files record by record, checks if indexes contain an N, and if it does places it in low quality,
    #checks if average read quality is above or equal to user input threshold, if not places record in low quality files,
    #checks if index 1 is the reverse compliment of index two, if not places in unmatched index files, if true places
    #in that barcode file. All headers are appended to include the barcodes at the start. Also includes counters  All
    # headers are appended to include the barcodes at the start. Returns counts for each index, unkown indexes,
    # and index hopped indices as well as creating files
    with gzip.open(r1, 'rt') as r1, gzip.open(r2, 'rt') as r2, gzip.open(i1, 'rt') as i1, gzip.open(i2 , 'rt') as i2:
        counts = {'unknown' : 0 ,'index_hopped' : 0 }
        ln = 0
        r1temp = []
        r2temp = []
        i1temp = []
        i2temp = []
        i2_rev = []
        for r1_line , r2_line , i1_line , i2_line in zip(r1,r2,i1,i2):
            ln += 1
            r1temp.append(r1_line.strip())
            r2temp.append(r2_line.strip())
            i1temp.append(i1_line.strip())
            i2temp.append(i2_line.strip())
            if ln % 4 == 2:
                i2_rev = revcomp(i2_line.strip())
            if ln % 4 == 0:
                if qual(q,r1temp[3]) and qual(q,r2temp[3]) and qual(q,i1temp[3]) and qual(q,i2temp[3]) is True:
                    if str(i1temp[1]) == str(i2_rev) and str(i1temp[1]) in bar:
                        with open(i1temp[1]+'_fw.fq', 'a') as fr, open(i1temp[1]+'_rv.fq', 'a') as rr:
                            headerf = header(r1temp[0], i1temp[1], i2temp[1])
                            headerr = header(r1temp[0], i1temp[1], i2temp[1])
                            fr.write(headerf + '\n' + r1temp[1] + '\n' + r1temp[2] + '\n' + r1temp[3] + '\n')
                            rr.write(headerr + '\n' + r2temp[1] + '\n' + r2temp[2] + '\n' + r2temp[3] + '\n')
                            if i1temp[1] in counts:
                                counts[i1temp[1]]+=1
                                r1temp = []
                                r2temp = []
                                i1temp = []
                                i2temp = []
                            else:
                                counts[i1temp[1]] = 1
                                r1temp = []
                                r2temp = []
                                i1temp = []
                                i2temp = []
                    elif str(i1temp[1]) != str(i2_rev) and str(i1temp[1]) in bar and str(i2_rev) in bar:
                        with open('Index_hopped_fw.fq', 'a') as fr, open('Index_hopped_rv.fq', 'a') as rr:
                            headerf = header(r1temp[0], i1temp[1], i2temp[1])
                            headerr = header(r1temp[0], i1temp[1], i2temp[1])
                            fr.write(headerf + '\n' + r1temp[1] + '\n' + r1temp[2] + '\n' + r1temp[3] + '\n')
                            rr.write(headerr + '\n' + r2temp[1] + '\n' + r2temp[2] + '\n' + r2temp[3] + '\n')
                            counts['index_hopped'] += 1
                            r1temp = []
                            r2temp = []
                            i1temp = []
                            i2temp = []
                    else:
                        with open('unknown_fw.fq', 'a') as fr, open('unknown_rv.fq', 'a') as rr:
                            headerf = header(r1temp[0], i1temp[1], i2temp[1])
                            headerr = header(r1temp[0], i1temp[1], i2temp[1])
                            fr.write(headerf + '\n' + r1temp[1] + '\n' + r1temp[2] + '\n' + r1temp[3] + '\n')
                            rr.write(headerr + '\n' + r2temp[1] + '\n' + r2temp[2] + '\n' + r2temp[3] + '\n')
                            counts['unknown'] += 1
                            r1temp = []
                            r2temp = []
                            i1temp = []
                            i2temp = []
                else:
                    with open ('unknown_fw.fq','a') as fr, open ('unknown_rv.fq','a') as rr:
                        headerf = header(r1temp[0],i1temp[1],i2temp[1])
                        headerr = header(r1temp[0],i1temp[1],i2temp[1])
                        fr.write(headerf + '\n' + r1temp[1] + '\n' + r1temp[2] + '\n' + r1temp[3] + '\n')
                        rr.write(headerr + '\n' + r2temp[1] + '\n' + r2temp[2] + '\n' + r2temp[3] + '\n')
                        counts['unknown'] += 1
                        r1temp = []
                        r2temp = []
                        i1temp = []
                        i2temp = []
    print ("Barcode\tCount")
    for item in counts:
        print (str(item)+'\t'+str(counts[item]))

    return counts
counts = demult(r1,r2,i1,i2,bar,q)
y=[]
x=[]
for item in counts:
    x.append(item)
    y.append(counts[item])
y =  y
x =  x
plt.bar(x, y, align='center', alpha=0.5)
plt.ylabel('Count')
plt.xlabel('Barcode')
plt.title('Counts for each barcode')
plt.savefig('counts.png')







