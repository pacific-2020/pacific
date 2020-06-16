#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 15 09:32:02 2020

PACIFIC takes a FASTA/FASTQ input file and predicts the presence of the following viruses and their relative sample proportions:
        SARS-CoV-2;
        128 taxonomic units from Influenza,
        5 species from Metapneumovirus,
        130 species from Rhinovirus, and
        11 species from Coronaviridae (non-SARS-CoV-2).                                

@author: Pablo Acera

"""

import argparse

parser = argparse.ArgumentParser(prog='PACIFIC v0.1', description=
                                 """ 
                                 PACIFIC takes a FASTA/FASTQ input file and predicts the presence of the following viruses and their relative sample proportions:
                                 SARS-CoV-2;
                                 128 taxonomic units from Influenza,
                                 5 species from Metapneumovirus,
                                 130 species from Rhinovirus, and
                                 11 species from Coronaviridae (non-SARS-CoV-2).
                                 
                                 We recommend that users use default parameters to ensure high accuracy.
                                 """, usage='python PACIFIC.py [options] -i <in.fa>|<in.fq> -m <model> -t <tokenizer> -l <label-maker>\nversion: %(prog)s')

OPTIONAL = parser._action_groups.pop()
REQUIRED = parser.add_argument_group('required arguments')

#Inputs
## CHANGE -m  -t -l -f to OPTIONAL and CREATE RELATIVE PATHS FOR THESE FILES

REQUIRED.add_argument("-i", "--input_file",
                      help="FASTA/FASTQ input file path",
                      metavar='\b',
                      required=True)

REQUIRED.add_argument("-m", "--model",
                      help="PACIFIC model file path",
                      metavar='\b',
                      required=True)

REQUIRED.add_argument("-t", "--tokenizer",
                      help="Tokenizer file path",
                      metavar='\b',
                      required=True)

REQUIRED.add_argument("-l", "--label_maker",
                      help="Label maker object file path",
                      metavar='\b',
                      required=True)

#arguments
OPTIONAL.add_argument("-f", "--file_type",
                      help='FASTA or FASTQ training file format [fasta]',
                      metavar='<fasta/fastq>',
                      default='fasta',
                      )

OPTIONAL.add_argument("-o", "--outputdir",
                      help='Path to output directory [.]',
                      metavar='<dir>',
                      default="")

#OPTIONAL.add_argument("-k", "--k_mers",
#                      help='K-mer number use to train the model [9]',
#                      default=9,
#                      type=int)

OPTIONAL.add_argument("-T", "--prediction_threshold",
                      help='Threshold/cutoff for predictions [0.95]',
                      metavar='<float>',
                      default=0.95,
                      type=int
                      )

OPTIONAL.add_argument("-c", "--chunk_size",
                      help='Number of reads per chunk [10000]',
                      metavar='<int>',
                      default=50000,
                      type=int
                      )                      

OPTIONAL.add_argument("-O", "--output_fasta",
                      help='If this option is "True", a FASTA file containing predictions for each read will be provided [False]',
                      default=False,
                      action='store_true'
                      )

OPTIONAL.add_argument('-v', '--version', 
                        action='version', 
                        version='%(prog)s')                        


parser._action_groups.append(OPTIONAL)

ARGS = parser.parse_args()

# Inputs
FILE_IN = ARGS.input_file
MODEL = ARGS.model
TOKENIZER = ARGS.tokenizer
LABEL_MAKER = ARGS.label_maker


# Arguments
K_MERS = ARGS.k_mers
MODEL = ARGS.model
FILE_TYPE = ARGS.file_type
OUTPUTDIR = ARGS.outputdir
THRESHOLD_PREDICTION = ARGS.prediction_threshold
OUTPUT_FASTA = ARGS.output_fasta
CHUNK_SIZE = ARGS.chunk_size

# import other packages
from Bio import SeqIO

import pickle
from keras.models import load_model
import random
import numpy as np
import pandas as pd
import tensorflow as tf
import sys
import os

def prepare_read(trancriptome, file_type):
    '''
    function will take tranciprtome and make reads
    '''
    fasta_sequences = SeqIO.parse(open(trancriptome),file_type)
    sequences = []
    names = []
    for fasta in fasta_sequences:
        name, sequence = fasta.id, str(fasta.seq)
        sequences.append(sequence)
        names.append(name)
    return sequences, names


def process_reads(sequences, length, kmer, names):
    '''
    '''
    r_reads = []
    new_names = []
    for i in enumerate(sequences):
        # check the reads does not contain weird characters
        if all(c in 'AGCT' for c in i[1].upper()) and len(i[1]) >= 150:
            read = i[1][:150]
            r_reads.append(' '.join(read[x:x+kmer].upper() for x in range(len(read) - kmer + 1)))
            new_names.append(names[i[0]])
    return r_reads, new_names


def main(file, size_lenght, k_mer_size, file_type):
    '''
    '''
   
    all_transcripts, names = prepare_read(file,
                                          file_type)
    reads, names = process_reads(all_transcripts, 
                                 size_lenght,
                                 k_mer_size,
                                 names)

    return all_transcripts, reads, names

def accuracy(labels, predictions):
    '''
    calculate accuracy
    '''
    try:
        if labels.shape != predictions.shape:
            print('labels and predictions does not have same dimentions')
            return False
        
        correct = 0
        for i in range(len(labels)):
            if labels[i] == predictions[i]:
                correct +=1
    except:
        return 0
    
    return correct/len(labels)


if __name__ == '__main__':

    seed_value = 42
    random.seed(seed_value)# 3. Set `numpy` pseudo-random generator at a fixed value
    np.random.seed(seed_value)# 4. Set `tensorflow` pseudo-random generator at a fixed value
    tf.random.set_seed(seed_value)# 5. For layers that introduce randomness like dropout, make sure to set seed values 
    
    
    config = tf.compat.v1.ConfigProto()
    config.gpu_options.allow_growth = True
    sess = tf.compat.v1.Session(config=config)
    
    
    model = load_model(MODEL)
    
    # Keras loading sequences tokenizer 
    with open(TOKENIZER, 'rb') as handle:
        tokenizer = pickle.load(handle)
        
    # loading label maker
    with open(LABEL_MAKER, 'rb') as handle:
        label_maker = pickle.load(handle)
    
    print()    
    print('Converting reads into k-mers...')
    print()
    
    # Convert reads into k-mers
    reads, kmer_sequences, names = main(FILE_IN,
                                            150,
                                            K_MERS,
                                            FILE_TYPE)
    
    sequences = tokenizer.texts_to_sequences(kmer_sequences)
    
    if not sequences:
        sys.exit('All sequences are smaler than 150bp. No predictions made')
          
        
    print()
    print('Making predictions...')
       
    predictions = model.predict(np.array(sequences))
    labels = label_maker.inverse_transform(np.array(predictions), threshold=THRESHOLD_PREDICTION)

    if OUTPUT_FASTA is True:
        print()
        fasta_name_out = 'output_'+os.path.split(FILE_IN)[1]
        print('Output FASTA file '+fasta_name_out)
        with open(fasta_name_out,'w') as output:
            for i in enumerate(names):
                print(i[1]+':'+str(max(predictions[i[0]]))+':'+labels[i[0]], file=output)
                print(reads[i[0]], file=output)
                
    
    print()
    print('Using '+str(THRESHOLD_PREDICTION)+' thresholds to filter predictions...')
    
    predictions_high_acc = []
    names_high_acc = []
    for i in enumerate(predictions):
        if max(i[1]) > THRESHOLD_PREDICTION:
            predictions_high_acc.append(i[1])
            names_high_acc.append(names[i[0]])

    labels = label_maker.inverse_transform(np.array(predictions_high_acc), threshold=THRESHOLD_PREDICTION)
    
    df = pd.DataFrame(predictions_high_acc, columns = ['Coronaviridae',
                                                       'Human',
                                                       'Influenza',
                                                       'Metapneumovirus',
                                                       'Rhinovirus',
                                                       'Sars_cov_2'])
    
    df['Read_id'] = names_high_acc
    
    df['Labels'] = labels
    
    cols = ['Read_id',
            'Coronaviridae',
            'Human',
            'Influenza',
            'Metapneumovirus',
            'Rhinovirus',
            'Sars_cov_2',
            'Labels'
            ]
    
    df = df[cols]
    
    Coronaviridae = len(labels[labels == 'Coronaviridae']) / len(labels) * 100
    Human = len(labels[labels == 'Human']) / len(labels) * 100
    Influenza =  len(labels[labels == 'Influenza']) / len(labels) * 100
    Metapneumovirus =  len(labels[labels == 'Metapneumovirus']) / len(labels) * 100
    Rhinovirus =  len(labels[labels == 'Rhinovirus']) / len(labels) * 100
    Sars_cov_2 =  len(labels[labels == 'Sars_cov_2']) / len(labels) * 100
    
    results = {'Influenza':   Influenza,
               'Coronaviridae': Coronaviridae,
               'Metapneumovirus': Metapneumovirus,
               'Rhinovirus':  Rhinovirus,
               'Sars_cov_2':  Sars_cov_2
               }
    
    
    print()
    print('Saving output file to ', FILE_OUT)
    
    df.to_csv(FILE_OUT, sep='\t')
    
    print('From a total of '+str(len(kmer_sequences))+' 150bp reads, '+
          str(len(kmer_sequences) - len(predictions_high_acc))+' predictions are below the'+\
          'threshold and were discarted from the results')
    
    
    print()
    print('Relative proportion of virus in the sample')
    print()
    
    # specify the number of discarted reads
    
    
    print('Coronaviridae: ', Coronaviridae)
    
    print('Human: ', Human)
    
    print('Influenza', Influenza)
    
    print('Metapneumovirus', Metapneumovirus)
    
    print('Rhinovirus', Rhinovirus)
    
    print('Sars_cov_2', Sars_cov_2)
    print()
    print('Virus group proportions that overpass the empirical threshold are: ')
    
    limit_detection = {'Influenza': 0.001,
                       'Coronaviridae': 0.009,
                       'Metapneumovirus': 0.001,
                       'Rhinovirus': 0.0242,
                       'Sars_cov_2': 0.017 
                      }
    
    virus_positive = []
    for virus in results:
        if results[virus] > limit_detection[str(virus)]:
            virus_positive.append((str(virus)))
    
    if not virus_positive:
        print('None')
    else:
        print(' '.join(virus_positive))
    
    print()
    print('Thank you for using PACIFIC =^)')
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    

















