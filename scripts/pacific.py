#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  1 15:45:58 2020

RNN classifier for COVID-19 reads

@author: labuser
"""

from Bio import SeqIO
first_record = next(SeqIO.parse("example.fasta", "fasta"))



