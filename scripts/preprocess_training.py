#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  1 10:25:06 2020

This script takes as input a illumina sam file.
The file already have te be filtered like this

$samtools view -F 1294 -h -q 60 file.bam > $filtered.sam

get rid of discordant reads, and make two samfiles, one with  COVID reads 
another with non_covid readsx

Run the script like:

$python3 preprocess_shotReads.py samfile <output_name> <type_file>

file_type is single or paired

@author: labuser
"""

import sys
import os

samfile = sys.argv[1]
typefile = sys.argv[2]

name = samfile[:10]

outfile_covid = "".join(samfile.split('.')[:-1])+'_covid.sam'

outfile_non_covid = "".join(samfile.split('.')[:-1])+'_non_covid.sam'


# If output file exists delete them
if os.path.exists(outfile_covid):
    os.remove(outfile_covid)
    
if os.path.exists(outfile_non_covid):
    os.remove(outfile_non_covid)

with open(samfile, 'r') as file_in:
    for line in file_in:
        line = line.rstrip().split('\t')
        # check the read map to the same contig
        if typefile == 'paired':
            if len(line) > 6:
                if line[6] != '=':
                    continue
        if line[0][0] == '@':
            with open(outfile_covid, 'a') as file_covid:
                file_covid.write("\t".join(line)+'\n')
            with open(outfile_non_covid, 'a') as file_non_covid:
                file_non_covid.write("\t".join(line)+'\n')
                continue
        if line[2] == 'MN908947.3':
            with open(outfile_covid, 'a') as file_covid:
                file_covid.write("\t".join(line)+'\n')
        else:
            with open(outfile_non_covid, 'a') as file_non_covid:
                file_non_covid.write("\t".join(line)+'\n')

