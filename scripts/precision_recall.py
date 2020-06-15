#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Apr  5 21:42:12 2020

Calculate: accuracy, recall, precision, FPR, TPR

Accuracy = correctly classified / total
recall = true positives / (true positives + false negative) says about the sensitivity
precision = true positives / (true positives + false positives) says something about specificity

@author: labuser
"""

from Bio import SeqIO
import random
import os

import numpy as np
import pandas as pd
from keras.preprocessing.sequence import pad_sequences
import tensorflow as tf

from numpy.random import seed
from tensorflow import set_random_seed
import pickle

from keras.models import load_model
import matplotlib.pyplot as plt
import seaborn as sns
sns.set_style('dark')   
        
##### Functions to test with real illumina reads

def prepare_read_illumina(trancriptome, file_type):
    '''
    function will take tranciprtome and make reads
    '''
    fasta_sequences = SeqIO.parse(open(trancriptome),file_type)
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


def main_illumina(file, number_reads, size_lenght, k_mer_size, file_type):
    '''
    '''
    all_transcripts = prepare_read_illumina(file, file_type)
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


def recall(labels, predictions):
    '''
    calculate accuracy or recall for 1 class
    recall = true positives / (true positives + false negative) says about the sensitivity
    '''
    if labels.shape != predictions.shape:
        print('labels and predictions does not have same dimentions')
        return False
    
    correct = 0
    for i in range(len(labels)):
        if labels[i] == predictions[i]:
            correct +=1
    
    return correct/len(labels)



def precision(labels, predictions, other_virus):
    '''
    calculate precision
    precision = true positives / (true positives + false positives) says something about specificity
    '''
    if labels.shape != predictions.shape:
        print('labels and predictions does not have same dimentions')
        return False
    
    correct = 0
    for i in enumerate(predictions):
        if labels[0] == predictions[i[0]]:
            correct +=1
    
    rest_virus = 0
    for i in enumerate(other_virus):
        if labels[0] == other_virus[i[0]]:
            rest_virus +=1

    return correct/(correct + rest_virus)



if __name__ == '__main__':
    
    seed_value = 42
    random.seed(seed_value)# 3. Set `numpy` pseudo-random generator at a fixed value
    np.random.seed(seed_value)# 4. Set `tensorflow` pseudo-random generator at a fixed value
    tf.set_random_seed(seed_value)# 5. For layers that introduce randomness like dropout, make sure to set seed values 
       
    config = tf.ConfigProto()
    config.gpu_options.allow_growth = True
    sess = tf.Session(config=config)
 
    # keras load model
    model = load_model("/media/labuser/Data/COVID-19_classifier/pacific/model/pacific.pacific_9mers.01.h5")
    
    # Keras loading sequences tokenizer 
    with open('/media/labuser/Data/COVID-19_classifier/pacific/model/tokenizer.01.pacific_9mers.pickle', 'rb') as handle:
        tokenizer = pickle.load(handle)
        
    # loading label maker
    with open('/media/labuser/Data/COVID-19_classifier/pacific/model/label_maker.01.pacific_9mers.pickle', 'rb') as handle:
        label_maker = pickle.load(handle)
    
    ## illumina reads 
    Cornidovirineae_path = '/media/labuser/Data/COVID-19_classifier/pacific/data/InSilicoSeq_reads/Cornidovirineae/novaseq_reads_Cornidoviridae_1M.fastq'
    Influenza_path = '/media/labuser/Data/COVID-19_classifier/pacific/data/InSilicoSeq_reads/Influenza/novaseq_reads_Influenza_1M.fastq'
    Metapneumovirus_path  = '/media/labuser/Data/COVID-19_classifier/pacific/data/InSilicoSeq_reads/Metapneumovirus/novaseq_reads_Metapneumovirus_1M.fastq'
    Rhinovirus_path = '/media/labuser/Data/COVID-19_classifier/pacific/data/InSilicoSeq_reads/Rhinovirus/novaseq_reads_Rhinovirus_1M.fastq'
    SARS_CoV_2_path = '/media/labuser/Data/COVID-19_classifier/pacific/data/InSilicoSeq_reads/Sars-CoV-2/novaseq_reads_sars-cov-2_1M.fastq'
    Human_path = '/media/labuser/Data/COVID-19_classifier/pacific/data/InSiicoSeq_reads/Human/novaseq_reads_Human_1M.fastq'
    
    kmer = 9
    
    influenza = main_illumina(Influenza_path, 500000, 150, kmer, 'fastq')
    Cornidovirineae = main_illumina(Cornidovirineae_path, 500000, 150, kmer, 'fastq')
    Metapneumovirus = main_illumina(Metapneumovirus_path, 500000, 150, kmer, 'fastq')
    Rhinovirus = main_illumina(Rhinovirus_path, 500000, 150, kmer, 'fastq')
    SARS_CoV_2 = main_illumina(SARS_CoV_2_path, 500000, 150, kmer, 'fastq')
    Human = main_illumina(Human_path, 500000, 150, kmer, 'fastq')
    
    max_length = 142
    influenza_reads = pad_sequences(tokenizer.texts_to_sequences(influenza), maxlen = max_length, padding = 'post')
    predictinos_influenza = model.predict(influenza_reads)
    
    Cornidovirineae_reads = pad_sequences(tokenizer.texts_to_sequences(Cornidovirineae), maxlen = max_length, padding = 'post')
    predictinos_Cornidovirineae = model.predict(Cornidovirineae_reads)
    
    Metapneumovirus_reads = pad_sequences(tokenizer.texts_to_sequences(Metapneumovirus), maxlen = max_length, padding = 'post')
    predictinos_Metapneumovirus = model.predict(Metapneumovirus_reads)
    
    Rhinovirus_reads = pad_sequences(tokenizer.texts_to_sequences(Rhinovirus), maxlen = max_length, padding = 'post')
    predictinos_Rhinovirus = model.predict(Rhinovirus_reads)
    
    SARS_CoV_2_reads = pad_sequences(tokenizer.texts_to_sequences(SARS_CoV_2), maxlen = max_length, padding = 'post')
    predictinos_SARS_CoV_2 = model.predict(SARS_CoV_2_reads)
    
    Human_reads = pad_sequences(tokenizer.texts_to_sequences(Human), maxlen = max_length, padding = 'post')
    predictinos_Human = model.predict(Human_reads)

    recall_results = {'Human': [], 'Influenza': [], 'Sars_cov_2':[],  'Rhinovirus': [], 'Metapneumovirus': [], 'Cornidovirineae': []}
    precision_results = {'Human': [], 'Influenza': [], 'Sars_cov_2':[],  'Rhinovirus': [], 'Metapneumovirus': [], 'Cornidovirineae': []}
    
    classes = ['Human', 'Influenza', 'Sars_cov_2',  'Rhinovirus', 'Metapneumovirus', 'Cornidovirineae']
    
    # make 100 experiments measuring accuracy. precision and recall per class
    for i in range(100):
        # first select a random number of reads per class
        tmp_human = predictinos_Human[np.random.randint(len(Human), size=random.randrange(500000)),:]
        tmp_influenza = predictinos_influenza[np.random.randint(len(influenza), size=random.randrange(500000)),:]
        tmp_SARS_CoV_2 = predictinos_SARS_CoV_2[np.random.randint(len(SARS_CoV_2), size=random.randrange(500000)),:]
        tmp_Rhinovirus = predictinos_Rhinovirus[np.random.randint(len(Rhinovirus), size=random.randrange(500000)),:]
        tmp_Metapneumovirus = predictinos_Metapneumovirus[np.random.randint(len(Metapneumovirus), size=random.randrange(500000)),:]
        tmp_Cornidovirineae = predictinos_Cornidovirineae[np.random.randint(len(Cornidovirineae), size=random.randrange(500000)),:]
        
        tmp_list = [tmp_human, tmp_influenza, tmp_SARS_CoV_2, tmp_Rhinovirus, tmp_Metapneumovirus, tmp_Cornidovirineae]
        
        for j in enumerate(classes):
             labels_predict = list(np.repeat(j[1],len(tmp_list[j[0]])))
             labels_predict = label_maker.transform(labels_predict)
             recall_results[j[1]] += [accuracy(np.argmax(labels_predict, axis=1), np.argmax(tmp_list[j[0]], axis=1))] 
             other_viruses = np.array([])
             for y in enumerate(classes):
                 if y[0] == j[0]:
                     continue
                 else:
                     if len(other_viruses) == 0:
                         other_viruses =  tmp_list[y[0]]
                     else:
                         other_viruses = np.concatenate((other_viruses, tmp_list[y[0]]))
             precision_results[j[1]] += [precision(np.argmax(labels_predict, axis=1),
                                                   np.argmax(tmp_list[j[0]], axis=1), 
                                                   np.argmax(other_viruses, axis=1)
                                                   )
                                        ]

    columns_precision = precision_results['Human'] +\
                        precision_results['Influenza']+\
                        precision_results['Sars_cov_2']+\
                        precision_results['Rhinovirus']+\
                        precision_results['Metapneumovirus']+\
                        precision_results['Cornidovirineae']
        
    columns_virus = ['Human']*100+\
                    ['Influenza']*100+\
                    ['Sars_cov_2']*100+\
                    ['Rhinovirus']*100+\
                    ['Metapneumovirus']*100+\
                    ['Cornidovirineae']*100
    
    df_precision = pd.DataFrame({'precisiion': columns_precision, 'virus':columns_virus})
    df_precision.to_csv('/media/labuser/Data/COVID-19_classifier/pacific/results/9-mers/precision_9mers.csv')
    
    df_precision = pd.read_csv('/media/labuser/Data/COVID-19_classifier/pacific/results/9-mers/precision_9mers.csv')
    
    sns.set()
    f, ax = plt.subplots(figsize=(13,9))
    b = sns.boxplot(x='virus', y='precisiion', data=df_precision)
    b.axes.set_title('Precision per class', fontsize = 25)
    plt.ylim(0.978, 1)
    b.tick_params(axis='y', labelsize=25)
    b.tick_params(axis='x', labelsize=25, rotation=45)
    b.set_ylabel("Precision",fontsize=25)
    plt.ylabel('Precision', fontsize=25)
    plt.savefig('/media/labuser/Data/COVID-19_classifier/pacific/results/9-mers/precision_9mers_zoom.pdf',
                format='pdf',
                dpi=1200,
                bbox_inches='tight', pad_inches=0)
    
    columns_recall = recall_results['Human'] +\
                     recall_results['Influenza']+\
                     recall_results['Sars_cov_2']+\
                     recall_results['Rhinovirus']+\
                     recall_results['Metapneumovirus']+\
                     recall_results['Cornidovirineae']
        
    df_recall = pd.DataFrame({'recall': columns_recall, 'virus':columns_virus})
    df_recall.to_csv('/media/labuser/Data/COVID-19_classifier/pacific/results/9-mers/recall_9mers.csv')
    df_recall = pd.read_csv('/media/labuser/Data/COVID-19_classifier/pacific/results/9-mers/recall_9mers.csv')
    
    
    sns.set()
    f, ax = plt.subplots(figsize=(13,9))
    b = sns.boxplot(x='virus', y='recall', data=df_recall)
    b.axes.set_title('Recall per class', fontsize = 25)
    b.tick_params(axis='y', labelsize=25)
    b.tick_params(axis='x', labelsize=25, rotation=45)
    b.set_ylabel("Recall",fontsize=25)
    plt.ylabel('Recall', fontsize=25)
    plt.savefig('/media/labuser/Data/COVID-19_classifier/pacific/results/9-mers/recall_9mers_zoom.pdf',
                format='pdf',
                dpi=1200,
                bbox_inches='tight', pad_inches=0)
    
    
    
    






















