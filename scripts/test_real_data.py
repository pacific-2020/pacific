#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 15 09:32:02 2020

@author: labuser
"""

import pandas as pd
import os
from Bio import SeqIO

import pickle
from keras.models import load_model

import random
import os

import keras.backend as K

from sklearn.preprocessing import LabelBinarizer
import numpy as np
from keras.preprocessing.text import Tokenizer
from keras.preprocessing.sequence import pad_sequences
import tensorflow as tf
import pickle
import time
from keras.models import load_model

import seaborn as sns
import matplotlib.pyplot as plt


def prepare_read(trancriptome):
    '''
    function will take tranciprtome and make reads
    '''
    fasta_sequences = SeqIO.parse(open(trancriptome),'fastq')
    sequences = []
    for fasta in fasta_sequences:
        name, sequence = fasta.id, str(fasta.seq)
        sequences.append(sequence)
    return sequences


def process_reads(sequences, length, kmer):
    '''
    '''
    r_reads = []
    for i in enumerate(sequences):
        # check the reads does not contain weird characters
        if all(c in 'AGCT' for c in i[1].upper()) and len(i[1]) >= 150:
            read = i[1][:150]
            r_reads.append(' '.join(read[x:x+kmer].upper() for x in range(len(read) - kmer + 1)))
    return r_reads


def main(directory, size_lenght, k_mer_size):
    '''
    '''
    files = os.listdir(directory)
    reads = []
    for file in files:
        all_transcripts = prepare_read(directory+'/'+file)
        reads += process_reads(all_transcripts, 
                               size_lenght,
                               k_mer_size)
    return reads 

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

def select_high_acc(predictions, labels, threshold):
    '''
    '''
    labels_high_acc = []
    predictions_high_acc = []
    for i in enumerate(predictions):
        if max(i[1]) > 0.9:
            predictions_high_acc.append(i[1].tolist())
            labels_high_acc.append(labels[i[0]].tolist())
    try:
        predictions_high_acc = np.argmax(predictions_high_acc, axis=1)
        labels_high_acc = np.argmax(labels_high_acc, axis=1)
    except:
        predictions_high_acc = [0]
        labels_high_acc = [0]
    return predictions_high_acc, labels_high_acc


if __name__ == '__main__':

    seed_value = 42
    random.seed(seed_value)# 3. Set `numpy` pseudo-random generator at a fixed value
    np.random.seed(seed_value)# 4. Set `tensorflow` pseudo-random generator at a fixed value
    tf.set_random_seed(seed_value)# 5. For layers that introduce randomness like dropout, make sure to set seed values 
       
    config = tf.ConfigProto()
    config.gpu_options.allow_growth = True
    sess = tf.Session(config=config)
    
    model = load_model("/media/labuser/Data/COVID-19_classifier/pacific/model/pacific.01.pacific_9mers.h5")
    
    # Keras loading sequences tokenizer 
    with open('/media/labuser/Data/COVID-19_classifier/pacific/model/tokenizer.01.pacific_9mers.pickle', 'rb') as handle:
        tokenizer = pickle.load(handle)
        
    # loading label maker
    with open('/media/labuser/Data/COVID-19_classifier/pacific/model/label_maker.01.pacific_9mers.pickle', 'rb') as handle:
        label_maker = pickle.load(handle)
    

    # Whole FASTQ experiment
    for folder in ['SRR11412228', 'SRR11412229', 'SRR11412230']:
    
        GISAID = main('/media/labuser/Data/COVID-19_classifier/pacific/data/non_synthetic/illumina/SARS_human/'+folder,
                      150, 
                      9)
        
        labels = list(np.repeat('Human',len(GISAID)))
        
        sequences_GISAID = tokenizer.texts_to_sequences(GISAID)
        labels_GISAID = label_maker.transform(labels)
        
        start = time.time()
        predictions_GISAID = model.predict(np.array(sequences_GISAID))
        print(time.time() - start)
        
        predictions_high_acc, labels_high_acc = select_high_acc(predictions_GISAID,
                                                                labels_GISAID,
                                                                0.95)
        
        discarted = (len(sequences_GISAID) - len(predictions_high_acc)) / len(predictions_high_acc) * 100
        Cornidovirineae = len(predictions_high_acc[predictions_high_acc == 0]) / len(predictions_high_acc) * 100
        Human =len(predictions_high_acc[predictions_high_acc == 1]) / len(predictions_high_acc) * 100
        Influenza = len(predictions_high_acc[predictions_high_acc == 2]) / len(predictions_high_acc) * 100
        Metapneumovirus = len(predictions_high_acc[predictions_high_acc == 3]) / len(predictions_high_acc) * 100
        Rhinovirus = len(predictions_high_acc[predictions_high_acc == 4]) / len(predictions_high_acc) * 100
        Sars_cov_2 =  len(predictions_high_acc[predictions_high_acc == 5]) / len(predictions_high_acc) * 100
        
        print('discarted: ', discarted,
              'Cornidovirineae: ', Cornidovirineae,
              'Human: ', Human,
              'Influenza', Influenza,
              'Metapneumovirus', Metapneumovirus,
              'Rhinovirus', Rhinovirus,
              'Sars_cov_2', Sars_cov_2)


    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    

















