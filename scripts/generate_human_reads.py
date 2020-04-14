#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Apr 12 11:40:36 2020

@author: labuser
"""

import Bio
import random
from Bio.Seq import Seq


def parse_fastq(trancriptome):
    '''
    function will take tranciprtome and make reads
    '''
    fasta_sequences = SeqIO.parse(open(trancriptome),'fasta')
    sequences = []
    for fasta in fasta_sequences:
        name, sequence = fasta.id, str(fasta.seq)
        sequences.append(sequence)
    return sequences

def rev_compl(sequence):
    seq = Seq(sequence)
    return str(seq.reverse_complement())


def random_reads(sequences, number_reads, length, file):
    '''
    '''
    counter = 0
    for read in range(number_reads):
        transcript = random.choice(sequences)
        try:
            start = random.randint(0, len(transcript) - length)
            read = transcript[start:start+length]
        except:
            continue
        if random.randint(0,1) == 1:
            read = rev_compl(read)
        with open(file, 'a') as file_out:
            print('>'+str(counter), file=file_out)
            print(read, file=file_out)
        counter +=1
    return True



def main(path, depth, out):
        
    path_in = path
    genome = parse_fastq(path_in)
    
    total_lenght = int(len(''.join(genome))/150)
    coverage = 5
    total_reads  = coverage * total_lenght
    
    path_out = out
    random_reads(human_cds, total_reads, 150, path_out)
    
    return True


    
##### Corniviridae

path_human = '/media/labuser/Data/COVID-19_classifier/pacific/data/custom_references/Rhinovirus/custom_reference_Rhinovirus.fasta'
human_cds = parse_fastq(path_human)

total_lenght = int(len(''.join(human_cds))/150)
coverage = 150
total_reads  = coverage * total_lenght
print('Total number of reads', total_reads)

path_out = '/media/labuser/Data/COVID-19_classifier/pacific/data/synthetic_reads/Human/human_synthetic.fasta'

random_reads(human_cds, total_reads, 150, path_out)



main('/media/labuser/Data/COVID-19_classifier/pacific/data/custom_references/Cornidovirineae/custom_reference_Cornidovirineae.fasta',
     50,
      '/media/labuser/Data/COVID-19_classifier/pacific/data/new_synthetic_reads/Cornidovirineae/Cornidovirineae_synthetic.fasta')


main('/media/labuser/Data/COVID-19_classifier/pacific/data/custom_references/Influenza/custom_reference_Influenza.fasta',
     80,
      '/media/labuser/Data/COVID-19_classifier/pacific/data/new_synthetic_reads/Influenza/Influenza_synthetic.fasta')


main('/media/labuser/Data/COVID-19_classifier/pacific/data/custom_references/Rhinovirus/custom_reference_Rhinovirus.fasta',
     150,
      '/media/labuser/Data/COVID-19_classifier/pacific/data/new_synthetic_reads/Rhinovirus/Rhinovirus_synthetic.fasta')






























