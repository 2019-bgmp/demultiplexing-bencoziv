#!/usr/bin/env python
import gzip
import argparse
import numpy as np
import matplotlib.pyplot as plt

def get_args():
    #get the input file name
    parser = argparse.ArgumentParser(description=' File name to be used')
    parser.add_argument('-f', '--file', type=str,help='File Name to be evaluated',required=True)
    return parser.parse_args()
args= get_args()
f = args.file
def length(f):
    # f is a file, determine the length of the read and output it
    with gzip.open(f, 'rt') as fh:
        ln=0
        for line in fh:
            ln+=1
            if ln%4==0:
                #checks if the line is a score line
                s=line.strip('\n')
                x=len(s)
                return x

x = length(f)
def scores(f,x):
    #f is a file, x is the length determined by the previous function, ouputs the mean quality score per nucleotide
    with gzip.open(f, 'rt') as fh:
        ln=0
        mean= np.zeros((x),dtype=float)
        #creates an array for containing mean scores based on x
        for line in fh:
            ln+=1
            if ln%4==0:
                #checks if the line is a quality score line
                x=0
                phred_score = line.strip('\n')
                for score in phred_score:
                    #looks at each quality score individually
                    n=(((ord(score)-33)))
                    mean[x]+=n
                    x+=1
        l=(ln/4)
        means=[]
        for i in range(x):
            #creates a new array containing the mean at each nucleotide position
            means.append(mean[i]/l)
    return means


mean = scores(f,x)


def fileFromPath(f):
    # removes file path, will replace with os when were allowed to use it to remove path from any file
    return f.strip('/projects/bgmp/shared/2017_sequencing/')
# graphs the mean scores and saves the table as file name .png to directory the file was run from


fp = fileFromPath(f)
y = mean
x = range(x)
plt.bar(x, y, align='center', alpha=0.5)
plt.ylabel('Mean Score')
plt.xlabel('Nucleotide Index')
plt.title('Mean Scores at Nucleotide Indexes');
plt.savefig(fp+'.png')
