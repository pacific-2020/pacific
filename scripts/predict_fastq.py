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

def test_predictions(virus_sequences, virus, threshold):
    '''
    '''
     # Test INfluenza reads
    labels_predict = list(np.repeat(virus,len(virus_sequences)))
    labels_predict = label_maker.transform(labels_predict)
    
    reads = tokenizer.texts_to_sequences(virus_sequences)
    max_length =150
    reads = pad_sequences(reads, maxlen = max_length, padding = 'post')

    predictions = model.predict(reads)
    
    # get rid of predictions where the best prediction is lower than 0.5
    labels_high_acc = []
    predictions_high_acc = []
    for i in enumerate(predictions):
        if max(i[1]) > threshold:
            predictions_high_acc.append(i[1].tolist())
            labels_high_acc.append(labels_predict[i[0]].tolist())
       
    predictions_high_acc = np.argmax(predictions_high_acc, axis=1)
    labels_high_acc = np.argmax(labels_high_acc, axis=1)
    
    accuracy_validation = accuracy(labels_high_acc, predictions_high_acc)
    print('Data from '+virus +' total reads ' + str(len(virus_sequences)))
    print('accuracy reads', accuracy_validation)
    print('Discarted reads ', predictions.shape[0] - len(predictions_high_acc))
    
    make_plot(predictions_high_acc, virus)
    
    return accuracy_validation, predictions.shape[0] - len(predictions_high_acc)
    

def make_plot(predictions_high_acc, virus_group):

    Cornidovirineae = len(predictions_high_acc[ predictions_high_acc == 0]) / len(predictions_high_acc) * 100
    Human =len(predictions_high_acc[ predictions_high_acc == 1]) / len(predictions_high_acc) * 100
    Influenza = len(predictions_high_acc[ predictions_high_acc == 2]) / len(predictions_high_acc) * 100
    Metapneumovirus = len(predictions_high_acc[ predictions_high_acc == 3]) / len(predictions_high_acc) * 100
    Rhinovirus = len(predictions_high_acc[ predictions_high_acc == 4]) / len(predictions_high_acc) * 100
    Sars_cov_2 =  len(predictions_high_acc[ predictions_high_acc == 5]) / len(predictions_high_acc) * 100
    
    
    x = np.array(['Human','Cornidovirineae', 'Influenza', 'Metapneumovirus',
       'Rhinovirus', 'Sars_cov_2']) 

    y = np.array([Human,  Cornidovirineae, Influenza, Metapneumovirus, Rhinovirus, Sars_cov_2])
  
    f, (ax, ax2) = plt.subplots(2, 1, sharex=True, figsize=(13,9))

    # plot the same data on both axes
    #ax.bar(x, y)
    #ax2.bar(x, y)
    ax = sns.barplot(x, y, ax=ax)
    ax2 = sns.barplot(x, y, ax=ax2)
    # zoom-in / limit the view to different portions of the data
    ax.set_ylim(90, 100)  # outliers only
    ax2.set_ylim(0, 8)  # most of the data
    
    # hide the spines between ax and ax2
    #ax.spines['bottom'].set_visible(False)
    #ax2.spines['top'].set_visible(False)
    ax.xaxis.tick_top()
    ax.tick_params(labeltop=False)  # don't put tick labels at the top
    ax2.xaxis.tick_bottom()
    
    # This looks pretty good, and was fairly painless, but you can get that
    # cut-out diagonal lines look with just a bit more work. The important
    # thing to know here is that in axes coordinates, which are always
    # between 0-1, spine endpoints are at these locations (0,0), (0,1),
    # (1,0), and (1,1).  Thus, we just need to put the diagonals in the
    # appropriate corners of each of our axes, and so long as we use the
    # right transform and disable clipping.
    
    d = .015  # how big to make the diagonal lines in axes coordinates
    # arguments to pass to plot, just so we don't keep repeating them
    kwargs = dict(transform=ax.transAxes, color='k', clip_on=False)
    ax.plot((-d, +d), (-d, +d), **kwargs)        # top-left diagonal
    ax.plot((1 - d, 1 + d), (-d, +d), **kwargs)  # top-right diagonal
    
    kwargs.update(transform=ax2.transAxes)  # switch to the bottom axes
    ax2.plot((-d, +d), (1 - d, 1 + d), **kwargs)  # bottom-left diagonal
    ax2.plot((1 - d, 1 + d), (1 - d, 1 + d), **kwargs)  # bottom-right diagonal
    
    # What's cool about this is that now if we vary the distance between
    # ax and ax2 via f.subplots_adjust(hspace=...) or plt.subplot_tool(),
    # the diagonal lines will move accordingly, and stay right at the tips
    # of the spines they are 'breaking'
   
    plt.title('FPR in '+virus_group+' synthetic reads')
    plt.ylabel('Percentage of predicted reads')
    plt.savefig('/media/labuser/Data/COVID-19_classifier/pacific/results/FPR_'+virus_group+'_0.5_illumina_synthetic.pdf',
                format='pdf',
                dpi=1200,
                bbox_inches='tight', pad_inches=0)
    print(y)
    return True

def make_heatmap(virus_sequences, virus, threshold):
    '''
    '''
    reads = tokenizer.texts_to_sequences(virus_sequences)
    max_length =150
    reads = pad_sequences(reads, maxlen = max_length, padding = 'post')
    predictions = model.predict(reads)
    
    # get rid of predictions where the best prediction is lower than 0.5
    predictions_high_acc = []
    for i in enumerate(predictions):
        if max(i[1]) > threshold:
            predictions_high_acc.append(i[1].tolist())
       
    predictions_high_acc = np.argmax(predictions_high_acc, axis=1)
    
    Cornidovirineae = len(predictions_high_acc[predictions_high_acc == 0]) / len(predictions_high_acc) * 100
    Human =len(predictions_high_acc[ predictions_high_acc == 1]) / len(predictions_high_acc) * 100
    Influenza = len(predictions_high_acc[ predictions_high_acc == 2]) / len(predictions_high_acc) * 100
    Metapneumovirus = len(predictions_high_acc[ predictions_high_acc == 3]) / len(predictions_high_acc) * 100
    Rhinovirus = len(predictions_high_acc[ predictions_high_acc == 4]) / len(predictions_high_acc) * 100
    Sars_cov_2 =  len(predictions_high_acc[ predictions_high_acc == 5]) / len(predictions_high_acc) * 100
    
    return [Human, Cornidovirineae, Influenza, Metapneumovirus, Rhinovirus, Sars_cov_2]



if __name__ == '__main__':
    
    seed_value = 42
    random.seed(seed_value)# 3. Set `numpy` pseudo-random generator at a fixed value
    np.random.seed(seed_value)# 4. Set `tensorflow` pseudo-random generator at a fixed value
    tf.set_random_seed(seed_value)# 5. For layers that introduce randomness like dropout, make sure to set seed values 
       
    config = tf.ConfigProto()
    config.gpu_options.allow_growth = True
    sess = tf.Session(config=config)
 
    # keras load model
    model = load_model("/media/labuser/Data/COVID-19_classifier/pacific/model/pacific.01.h5")
    
    # Keras loading sequences tokenizer 
    with open('/media/labuser/Data/COVID-19_classifier/pacific/model/tokenizer.01.pickle', 'rb') as handle:
        tokenizer = pickle.load(handle)
        
    # loading label maker
    with open('/media/labuser/Data/COVID-19_classifier/pacific/model/label_maker.01.pickle', 'rb') as handle:
        label_maker = pickle.load(handle)
    
    influenza_path = '/media/labuser/Data/COVID-19_classifier/pacific/data/custom_references/'+ \
                     'Influenza/all_genome/SRR7841658_filtered.fasta'
    
    '''
    ## Real reads
    Cornidovirineae_path = '/media/labuser/Data/COVID-19_classifier/pacific/data/non_synthetic/illumina/Cornidovirineae/alignment/SRR3742834_filter.fastq'
    Influenza_path = '/media/labuser/Data/COVID-19_classifier/pacific/data/non_synthetic/illumina/Influenza/alingment/SRR1577743_filtered.fastq'
    Metapneumovirus_path  = '/media/labuser/Data/COVID-19_classifier/pacific/data/non_synthetic/illumina/Metapneumovirus/alignment/SRR8787081.fastq'
    Rhinovirus_path = '/media/labuser/Data/COVID-19_classifier/pacific/data/non_synthetic/illumina/Rhinovirus/alignment/SRR8356904_filtered.fastq'
    SARS_CoV_2_path = '/media/labuser/Data/COVID-19_classifier/pacific/data/non_synthetic/illumina/SARS_human/processed/all_filtered_covid.sam.fastq'
    Human_path = '/media/labuser/Data/COVID-19_classifier/pacific/data/non_synthetic/illumina/SARS_human/processed/all_filtered_non_covid.sam.fastq'
    exp_SRR11412227 = '/media/labuser/Data/COVID-19_classifier/pacific/data/non_synthetic/illumina/SARS_human/SRR11412227.fastq'
    '''
    
    ## illumina reads 
    Cornidovirineae_path = '/media/labuser/Data/COVID-19_classifier/pacific/data/InSilicoSeq_reads/Cornidovirineae/novaseq_reads_Cornidoviridae_1M.fastq'
    Influenza_path = '/media/labuser/Data/COVID-19_classifier/pacific/data/InSilicoSeq_reads/Influenza/novaseq_reads_Influenza_1M.fastq'
    Metapneumovirus_path  = '/media/labuser/Data/COVID-19_classifier/pacific/data/InSilicoSeq_reads/Metapneumovirus/novaseq_reads_Metapneumovirus_1M.fastq'
    Rhinovirus_path = '/media/labuser/Data/COVID-19_classifier/pacific/data/InSilicoSeq_reads/Rhinovirus/novaseq_reads_Rhinovirus_1M.fastq'
    SARS_CoV_2_path = '/media/labuser/Data/COVID-19_classifier/pacific/data/InSilicoSeq_reads/Sars-CoV-2/novaseq_reads_sars-cov-2_1M.fastq'
    Human_path = '/media/labuser/Data/COVID-19_classifier/pacific/data/InSilicoSeq_reads/Human/novaseq_reads_Human_1M.fastq'
    
    influenza = main_illumina(Influenza_path, 300000, 150, 4, 'fastq')
    Cornidovirineae = main_illumina(Cornidovirineae_path, 300000, 150, 4, 'fastq')
    Metapneumovirus = main_illumina(Metapneumovirus_path, 300000, 150, 4, 'fastq')
    Rhinovirus = main_illumina(Rhinovirus_path, 300000, 150, 4, 'fastq')
    SARS_CoV_2 = main_illumina(SARS_CoV_2_path, 300000, 150, 4, 'fastq')
    Human = main_illumina(Human_path, 300000, 150, 4, 'fastq')
    
    
    column_Influenza = make_heatmap(influenza, 'Influenza', 0.5)
    column_Cornidovirineae = make_heatmap(Cornidovirineae, 'Cornidovirineae', 0.5)
    column_Metapneumovirus = make_heatmap(Metapneumovirus, 'Metapneumovirus', 0.5)
    column_Rhinovirus = make_heatmap(Rhinovirus, 'Rhinovirus', 0.5)
    column_Sars_cov_2 = make_heatmap(SARS_CoV_2, 'Sars_cov_2', 0.5)
    column_Human = make_heatmap(Human, 'Human', 0.5)
    
    df_PFR = pd.DataFrame({'Human':column_Human,
                           'Cornidovirineae': column_Cornidovirineae,
                           'Influenza': column_Influenza,
                           'Metapneumovirus':column_Metapneumovirus,
                           'Rhinovirus': column_Rhinovirus,
                           'Sars_cov_2':column_Sars_cov_2})
    
    
    

    test_predictions(influenza, 'Influenza', 0.5)
    test_predictions(Cornidovirineae, 'Cornidovirineae', 0.5)
    test_predictions(Metapneumovirus, 'Metapneumovirus', 0.5)
    test_predictions(Rhinovirus, 'Rhinovirus', 0.5)
    test_predictions(SARS_CoV_2, 'Sars_cov_2', 0.5)
    test_predictions(Human, 'Human', 0.5)

























