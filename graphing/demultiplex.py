#!/usr/bin/env python3

import argparse
import gzip
import numpy as np

def get_arguments():
	parser = argparse.ArgumentParser (description = "graphing fastq scores")
	parser.add_argument("-f", "--file", help="input file", required=True, type=str)
	parser.add_argument("-t", "--fType", help="type of input file (index or read)", type=str, default = "read")
	parser.add_argument("-n", "--name", help="name of graph",required=True, type=str)
	return parser.parse_args()
args = get_arguments()


def convert_phred(letter):
    """Converts a single character into a phred score"""
    return(ord(letter)-33)

if args.fType == "index":
    bp = 8
else:
    bp = 101

assert convert_phred("A") == 32, "Phred score incorrect"
assert convert_phred("@") == 31, "Phred score incorrect"
assert convert_phred("#") == 2, "Phred score incorrect"

stats = np.zeros((bp, 2), dtype = float)

with gzip.open (args.file, "rt") as fh:
    i = 1
    LN = 0
    for line in fh:
        scoreCount = 0
        if i%4 == 0:
            line = line.strip('\n')
            while scoreCount < len(line):
                stats[scoreCount][0] = scoreCount
                stats[scoreCount][1] += ((convert_phred(line[scoreCount])))
                scoreCount+=1 
            LN += 1
        i+=1
	
fileName = args.name + ".txt"
filePath = "/projects/bgmp/maddyg/demultiplex/files/%s" % fileName
count = 0
	
with open (filePath, "w") as out:
	out.write(args.name)
	out.write("\n")
	while count < bp:
		out.write(str(int(stats[count][0])))
		out.write("\t")
		out.write(str(stats[count][1]/LN))
		out.write("\n")
		count += 1




