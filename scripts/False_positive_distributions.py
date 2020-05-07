#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Apr  5 21:42:12 2020

This script will take different FASTQ files and plot the porcentage of predicted reads per class
@author: labuser
"""

from Bio import SeqIO
import random
import os

import numpy as np
import pandas as pd

from keras.preprocessing.sequence import pad_sequences
import tensorflow as tf

import random
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


def test_false_positives(virus_class, predictions_class, predictions_outside_class, plot=False):
    '''
    '''
    label = np.argmax(label_maker.transform([virus_class]))
    true_positives = []
    for i in enumerate(predictions_class):
        if  np.argmax(i[1]) == label: # if True positive
            true_positives.append(max(i[1]))
    
    false_positives = []
    for i in enumerate(predictions_outside_class):
        if  np.argmax(i[1]) == label: # if True Negative
            false_positives.append(max(i[1]))

    if plot :
        f, ax = plt.subplots(figsize=(13,9))
        plt.title('True positives vs false positives '+str(virus_class))
        sns.distplot(true_positives, kde=False, bins=50, label='True positives')
        sns.distplot(false_positives, kde=False, bins=50, label='false positives')
        plt.xlabel('Predicted probabilities')
        plt.legend()
        plt.ylim(0, 1000)
        plt.xlim(0.7, 1)
        plt.savefig('/media/labuser/Data/COVID-19_classifier/pacific/results/9-mers/FPR_'+virus_class+'_0.5_illumina_synthetic_distributions_large.pdf',
                    format='pdf',
                    dpi=1200,
                    bbox_inches='tight', pad_inches=0)
    return true_positives, false_positives


def percentile_proportion(virus_class, predictions_outside_class, threshold):
    '''
    '''
    label = np.argmax(label_maker.transform([virus_class]))
    false_positives = 0
    Total = 0
    for i in enumerate(predictions_outside_class):
        if np.max(i[1]) >= threshold:
            if  np.argmax(i[1]) == label: # if True Negative
                false_positives += 1
            else:
                Total +=1

    return (100/Total)*false_positives

 
def proportion_distribution(virus_label, virus_group, iterations, threshold):
    '''
    '''
    proportions = []
    for i in range(iterations):
        #first total number of reads for predictions between 100K and 900K
        total_reads_n = random.randint(1000, 900000)
        idx = np.random.randint(900000, size=total_reads_n)
        total_reads =  virus_group[idx,:]
        proportions.append(percentile_proportion(virus_label,
                                                 total_reads,
                                                 threshold))
    return proportions


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
    Human_path = '/media/labuser/Data/COVID-19_classifier/pacific/data/InSilicoSeq_reads/Human/novaseq_reads_Human_1M.fastq'
    
    influenza = main_illumina(Influenza_path, 100000, 150, 4, 'fastq')
    Cornidovirineae = main_illumina(Cornidovirineae_path, 100000, 150, 4, 'fastq')
    Metapneumovirus = main_illumina(Metapneumovirus_path, 100000, 150, 4, 'fastq')
    Rhinovirus = main_illumina(Rhinovirus_path, 100000, 150, 4, 'fastq')
    SARS_CoV_2 = main_illumina(SARS_CoV_2_path, 100000, 150, 4, 'fastq')
    Human = main_illumina(Human_path, 100000, 150, 4, 'fastq')
    
    max_length = 142

    influenza_reads = pad_sequences(tokenizer.texts_to_sequences(influenza), maxlen = max_length, padding = 'post')
    Cornidovirineae_reads = pad_sequences(tokenizer.texts_to_sequences(Cornidovirineae), maxlen = max_length, padding = 'post')
    Metapneumovirus_reads = pad_sequences(tokenizer.texts_to_sequences(Metapneumovirus), maxlen = max_length, padding = 'post' )
    Rhinovirus_reads = pad_sequences(tokenizer.texts_to_sequences(Rhinovirus), maxlen = max_length, padding = 'post')
    SARS_CoV_2_reads = pad_sequences(tokenizer.texts_to_sequences(SARS_CoV_2), maxlen = max_length, padding = 'post')
    Human_reads = pad_sequences(tokenizer.texts_to_sequences(Human), maxlen = max_length, padding = 'post')
    
    predictinos_influenza = model.predict(influenza_reads)
    predictinos_Cornidovirineae = model.predict(Cornidovirineae_reads)
    predictinos_Metapneumovirus = model.predict(Metapneumovirus_reads)
    predictinos_Rhinovirus = model.predict(Rhinovirus_reads)
    predictinos_SARS_CoV_2 = model.predict(SARS_CoV_2_reads)
    predictinos_Human = model.predict(Human_reads)
    
    
    # Make plots True positive vs False positives
    true_influenza, false_influenza = test_false_positives('Influenza', 
                                                           predictinos_influenza, 
                                                           np.concatenate((
                                                                           predictinos_Cornidovirineae,
                                                                           predictinos_Metapneumovirus,
                                                                           predictinos_Rhinovirus,
                                                                           predictinos_SARS_CoV_2,
                                                                           predictinos_Human),axis=0),
                                                           'plot')
    
    true_Cornidovirineae, false_Cornidovirineae = test_false_positives('Cornidovirineae', 
                                                                       predictinos_Cornidovirineae, 
                                                                       np.concatenate((
                                                                           predictinos_influenza,
                                                                           predictinos_Metapneumovirus,
                                                                           predictinos_Rhinovirus,
                                                                           predictinos_SARS_CoV_2,
                                                                           predictinos_Human),axis=0),
                                                                        'plot')
    
    true_Metapneumovirus, false_Metapneumovirus = test_false_positives('Metapneumovirus', 
                                                                       predictinos_Metapneumovirus, 
                                                                       np.concatenate((
                                                                           predictinos_influenza,
                                                                           predictinos_Cornidovirineae,
                                                                           predictinos_Rhinovirus,
                                                                           predictinos_SARS_CoV_2,
                                                                           predictinos_Human),axis=0),
                                                                               'plot')
    
    true_Rhinovirus, false_Rhinovirus = test_false_positives('Rhinovirus', 
                                                              predictinos_Rhinovirus, 
                                                              np.concatenate((
                                                                           predictinos_influenza,
                                                                           predictinos_Cornidovirineae,
                                                                           predictinos_Metapneumovirus,
                                                                           predictinos_SARS_CoV_2,
                                                                           predictinos_Human),axis=0),
                                                                      'plot')
    
    true_SARS_CoV_2, false_SARS_CoV_2 = test_false_positives('Sars_cov_2',
                                                             predictinos_SARS_CoV_2,
                                                             np.concatenate((
                                                                           predictinos_influenza,
                                                                           predictinos_Cornidovirineae,
                                                                           predictinos_Metapneumovirus,
                                                                           predictinos_Rhinovirus,
                                                                           predictinos_Human),axis=0),
                                                                     'plot')
    
    true_Human, false_Human = test_false_positives('Human',
                                                   predictinos_Human,
                                                   np.concatenate((
                                                                   predictinos_influenza,
                                                                   predictinos_Cornidovirineae,
                                                                   predictinos_Metapneumovirus,
                                                                   predictinos_Rhinovirus,
                                                                   predictinos_SARS_CoV_2),axis=0),
                                                           'plot')
                                                   
                                                   
    ## 
    
    proportions_Influenza = proportion_distribution('Influenza',
                                                     predictinos_Human,
                                                     100,
                                                     0.95)
    
    proportions_Sars_cov_2 = proportion_distribution('Sars_cov_2',
                                                     predictinos_Cornidovirineae,
                                                     100,
                                                     0.95)
    
    
    proportions_Cornidovirineae = proportion_distribution('Cornidovirineae',
                                                           predictinos_Human,
                                                            100,
                                                            0.95)

    proportions_Rhinovirus = proportion_distribution('Rhinovirus',
                                                     predictinos_Human,
                                                     100,
                                                     0.95)
    
    proportions_Metapneumovirus = proportion_distribution('Metapneumovirus',
                                                           predictinos_Human,
                                                            100,
                                                            0.95)
    
    influenza_name = ['Influenza']*100
    Cornidovirineae_name = ['Cornidovirineae']*100
    Metapneumovirus_name = ['Metapneumovirus']*100
    Rhinovirus_name = ['Rhinovirus']*100
    Sars_cov_2_name = ['Sars_cov_2']*100
    
    virus = influenza_name + Cornidovirineae_name + Metapneumovirus_name + Rhinovirus_name + Sars_cov_2_name
    
    proportions = proportions_Influenza  + \
                  proportions_Cornidovirineae + \
                  proportions_Metapneumovirus +\
                  proportions_Rhinovirus +\
                  proportions_Sars_cov_2
    
    df_proportions = pd.DataFrame({'virus' : virus, 'FPR in the top 95 percentil': proportions})
    
    df_proportions.to_csv('/media/labuser/Data/COVID-19_classifier/pacific/results/FPR_0.95_illumina_in_silico_100_experiments_distributions.csv')
    
    f, ax = plt.subplots(figsize=(13,9))
    plt.title('FPR in-silico 100 experiments per class')
    sns.boxplot(x="virus", y='FPR in the top 95 percentil', data=df_proportions)
    sns.swarmplot(x="virus", y='FPR in the top 95 percentil', data=df_proportions, color=".25")
    
    plt.savefig('/media/labuser/Data/COVID-19_classifier/pacific/results/FPR_0.95_illumina_in_silico_100_experiments_distributions_boxplots.pdf',
                    format='pdf',
                    dpi=1200,
                    bbox_inches='tight', pad_inches=0)
    
    print(max(df_proportions[df_proportions['virus'] =='Influenza']['FPR in the top 95 percentil']))
    print(max(df_proportions[df_proportions['virus'] =='Cornidovirineae']['FPR in the top 95 percentil']))
    print(max(df_proportions[df_proportions['virus'] =='Metapneumovirus']['FPR in the top 95 percentil']))
    print(max(df_proportions[df_proportions['virus'] =='Rhinovirus']['FPR in the top 95 percentil']))
    print(max(df_proportions[df_proportions['virus'] =='Sars_cov_2']['FPR in the top 95 percentil']))

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
