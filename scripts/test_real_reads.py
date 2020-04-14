#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 13 14:42:04 2020

Load trainned model, grab illumina reads and test the model

@author: labuser
"""


from Bio import SeqIO
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



##### Functions to test with real illumina reads

def prepare_read_illumina(trancriptome):
    '''
    function will take tranciprtome and make reads
    '''
    fasta_sequences = SeqIO.parse(open(trancriptome),'fastq')
    sequences = []
    for fasta in fasta_sequences:
        name, sequence = fasta.id, str(fasta.seq)
        sequences.append(sequence)
    return sequences


def process_reads_illumina(sequences, number_reads, length, kmer):
    '''
    '''
    r_reads = []
    for i in enumerate(sequences):
        if len(i[1]) >= length and 'N' not in i[1]:
            read = i[1][:length]
            r_reads.append(' '.join(read[x:x+kmer].upper() for x in range(len(read) - kmer + 1)))
            if len(r_reads) == number_reads:
                break
    return r_reads


def main_illumina(file, number_reads, size_lenght, k_mer_size):
    '''
    '''
    all_transcripts = prepare_read_illumina(file)
    reads = process_reads_illumina(all_transcripts, 
                                   number_reads, 
                                   size_lenght,
                                   k_mer_size)
    return reads 


def accuracy(labels, predictions):
    '''
    calculate accuracy
    '''
    if labels.shape != predictions.shape:
        print('labels and predictions does not have same dimentions')
        return False
    
    correct = 0
    for i in range(len(labels)):
        if labels[i] == predictions[i]:
            correct +=1
    
    return correct/len(labels)



if __name__ == '__main__':

    # keras load model
    model = load_model("/media/labuser/Data/COVID-19_classifier/pacific/model/pacific.h5")
    
    # Keras loading sequences tokenizer 
    with open('/media/labuser/Data/COVID-19_classifier/pacific/model/tokenizer.pickle', 'rb') as handle:
        tokenizer = pickle.load(handle)
        
    # loading label maker
    with open('/media/labuser/Data/COVID-19_classifier/pacific/model/label_maker.pickle', 'rb') as handle:
        label_maker = pickle.load(handle)
    
    #### Test with real data
    
    # short illumina covid reads
    path = '/media/labuser/Data/COVID-19_classifier/pacific/data/non_synthetic/illumina/processed/'
    files = os.listdir(path)
    files_covid = [i for i in files if i.endswith('filtered_covid.sam.fastq')]
    
    # take covid reads
    covid = main_illumina(path + files_covid[0], 10000, 150, 4)
    covid2 = main_illumina(path + files_covid[1], 10000, 150, 4)
    covid3 = main_illumina(path + files_covid[2], 10000, 150, 4)
    covid4 = main_illumina(path + files_covid[3], 10000, 150, 4)
    
    covid_reads = covid +covid2 + covid3 + covid4
    
    # Test covid reads    
    labels_predict = list(np.repeat('Sars_cov_2',len(covid_reads)))
    labels_predict = label_maker.transform(labels_predict)
    
    covid_reads = tokenizer.texts_to_sequences(covid_reads)
    max_length =150
    covid_reads = pad_sequences(covid_reads, maxlen = max_length, padding = 'post')

    predictions_covid = model.predict(covid_reads)
    
    # get rid of predictions where the best prediction is lower than 0.5
    labels_high_acc = []
    predictions_high_acc = []
    for i in enumerate(predictions_covid):
        if max(i[1]) > 0.9:
            predictions_high_acc.append(i[1].tolist())
            labels_high_acc.append(labels_predict[i[0]].tolist())
       
    predictions_high_acc = np.argmax(predictions_high_acc, axis=1)
    labels_high_acc = np.argmax(labels_high_acc, axis=1)
    
    accuracy_validation = accuracy(labels_high_acc, predictions_high_acc)
    print('accuracy covid reads', accuracy_validation)
    print('Discarted reads ', predictions_covid.shape[0] - len(predictions_high_acc))
    
    ####################################==============================================
    
    ### test with Human data
    path = '/media/labuser/Data/COVID-19_classifier/pacific/data/non_synthetic/illumina/processed/'
    files = os.listdir(path)
    files_non_covid = [i for i in files if i.endswith('non_covid.sam.fastq')]
    
    # Lets try 20k reads for training each
    non_covid = main_illumina(path + files_non_covid[0], 5000, 150, 4)
    non_covid2 = main_illumina(path + files_non_covid[1], 5000, 150, 4)
    
    host_reads = non_covid + non_covid2
    
    # Test host reads
     # Test covid reads    
    labels_predict =list(np.repeat('Human',len(host_reads)))
    labels_predict = label_maker.transform(labels_predict)
    
    host_reads = tokenizer.texts_to_sequences(host_reads)
    
    host_reads = pad_sequences(host_reads, maxlen = max_length, padding = 'post')

    predictions_host = model.predict(host_reads)
    
    # get rid of predictions where the best prediction is lower than 0.5
    labels_high_acc = []
    predictions_high_acc = []
    for i in enumerate(predictions_host):
        if max(i[1]) > 0.9:
            predictions_high_acc.append(i[1].tolist())
            labels_high_acc.append(labels_predict[i[0]].tolist())
       
    predictions_high_acc = np.argmax(predictions_high_acc, axis=1)
    labels_high_acc = np.argmax(labels_high_acc, axis=1)
    
    accuracy_validation = accuracy(labels_high_acc, predictions_high_acc)
    print('accuracy human reads', accuracy_validation)
    print('Discarted reads ', predictions_host.shape[0] - len(predictions_high_acc))
    
    
    
    
    
    
    