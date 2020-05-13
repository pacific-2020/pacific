#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 29 10:14:18 2020

@author: labuser
"""


from Bio import SeqIO
import random
import os

import numpy as np
import pandas as pd
from scipy import stats

from keras.preprocessing.sequence import pad_sequences
import tensorflow as tf

from numpy.random import seed
from tensorflow import set_random_seed
import pickle

from keras.models import load_model
import matplotlib.pyplot as plt
import seaborn as sns
sns.set_style('dark')   
sns.set()
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


def percentile_proportion(virus_class, predictions, threshold):
    '''
    '''
    label = np.argmax(label_maker.transform([virus_class]))
    false_positives = 0
    Total = 0
    for i in enumerate(predictions):
        if np.max(i[1]) >= threshold:
            if  np.argmax(i[1]) == label: # if True Negative
                false_positives += 1
            else:
                Total +=1
    
   
    return (100/Total)*false_positives




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
        
    '''
    ## illumina reads miseq
    Cornidovirineae_path = '/media/labuser/Data/COVID-19_classifier/pacific/data/InSilicoSeq_reads/Cornidovirineae/miseq/miseq_reads_Cornidovirineae.fastq'
    Influenza_path = '/media/labuser/Data/COVID-19_classifier/pacific/data/InSilicoSeq_reads/Influenza/miseq/miseq_reads_Influenza.fastq'
    Metapneumovirus_path  = '/media/labuser/Data/COVID-19_classifier/pacific/data/InSilicoSeq_reads/Metapneumovirus/miseq/miseq_reads_Metapneumovirus.fastq'
    Rhinovirus_path = '/media/labuser/Data/COVID-19_classifier/pacific/data/InSilicoSeq_reads/Rhinovirus/miseq/miseq_reads_rhinovirus.fastq'
    SARS_CoV_2_path = '/media/labuser/Data/COVID-19_classifier/pacific/data/InSilicoSeq_reads/Sars-CoV-2/miseq/miseq_reads_sars-cov-2.fastq'
    Human_path = '/media/labuser/Data/COVID-19_classifier/pacific/data/InSilicoSeq_reads/Human/miseq/miseq_reads_human.fastq'
    '''
    
    ## illumina reads novaseq
    Cornidovirineae_path = '/media/labuser/Data/COVID-19_classifier/pacific/data/InSilicoSeq_reads/Cornidovirineae/novaseq_reads_Cornidoviridae_1M.fastq'
    Influenza_path = '/media/labuser/Data/COVID-19_classifier/pacific/data/InSilicoSeq_reads/Influenza/novaseq_reads_Influenza_1M.fastq'
    Metapneumovirus_path  = '/media/labuser/Data/COVID-19_classifier/pacific/data/InSilicoSeq_reads/Metapneumovirus/novaseq_reads_Metapneumovirus_1M.fastq'
    Rhinovirus_path = '/media/labuser/Data/COVID-19_classifier/pacific/data/InSilicoSeq_reads/Rhinovirus/novaseq_reads_Rhinovirus_1M.fastq'
    SARS_CoV_2_path = '/media/labuser/Data/COVID-19_classifier/pacific/data/InSilicoSeq_reads/Sars-CoV-2/novaseq_reads_sars-cov-2_1M.fastq'
    Human_path = '/media/labuser/Data/COVID-19_classifier/pacific/data/InSilicoSeq_reads/Human/novaseq_reads_Human_1M.fastq'
    
    
    ## Real reads 
    Cornidovirineae_path = '/media/labuser/Data/COVID-19_classifier/pacific/data/non_synthetic/illumina/Cornidovirineae/alignment/SRR3742834_filter.fastq'
    Influenza_path = '/media/labuser/Data/COVID-19_classifier/pacific/data/non_synthetic/illumina/Influenza/alingment/SRR1577743_filtered.fastq'
    Metapneumovirus_path  = '/media/labuser/Data/COVID-19_classifier/pacific/data/non_synthetic/illumina/Metapneumovirus/alignment/SRR8787081.fastq'
    Rhinovirus_path = '/media/labuser/Data/COVID-19_classifier/pacific/data/non_synthetic/illumina/Rhinovirus/alignment/SRR8356904_filtered.fastq'
    SARS_CoV_2_path = '/media/labuser/Data/COVID-19_classifier/pacific/data/non_synthetic/illumina/SARS_human/processed/all_filtered_covid.sam.fastq'
    Human_path = '/media/labuser/Data/COVID-19_classifier/pacific/data/non_synthetic/illumina/SARS_human/processed/all_filtered_non_covid.sam.fastq'
    
    kmer = 9
    
    influenza = main_illumina(Influenza_path, 50000, 150, kmer, 'fastq')
    Cornidovirineae = main_illumina(Cornidovirineae_path, 50000, 150, kmer, 'fastq')
    SARS_CoV_2 = main_illumina(SARS_CoV_2_path, 50000, 150, kmer, 'fastq')
    Metapneumovirus = main_illumina(Metapneumovirus_path, 50000, 150, kmer, 'fastq')
    Rhinovirus = main_illumina(Rhinovirus_path, 50000, 150, kmer, 'fastq')
    Human = main_illumina(Human_path, 500000, 150, kmer, 'fastq')
    
    classes_dic = {'influenza': influenza,
                   'Cornidovirineae': Cornidovirineae, 
                   'Metapneumovirus': Metapneumovirus, 
                   'Rhinovirus': Rhinovirus,
                   'Sars_cov_2': SARS_CoV_2
                   }
         
    classes= ['influenza',
              'Cornidovirineae', 
              'Metapneumovirus', 
              'Rhinovirus',
              'Sars_cov_2'
               ]

    cutoffs = {'proportion_Cornidovirineae': 0.250501,
               'proportion_Metapneumovirus': 0.0335,
               'proportion_Rhinovirus': 0.38288,
               'proportion_SARS_CoV_2': 0.32697,
               'proportion_Influenza': 0.0751592
                }
    
    max_length = 142
    Human_reads = pad_sequences(tokenizer.texts_to_sequences(Human), maxlen = max_length, padding = 'post')
    predictinos_Human = model.predict(Human_reads)
    percentage_results = {}
    
    percentages = [ 2.5, 1, 0.5]

    for percentage in percentages:
        for virus in classes:
            print(virus)
            virus_amount = int(len(predictinos_Human)/100*percentage)
            virus_reads = classes_dic[virus][:virus_amount]
            virus_reads = pad_sequences(tokenizer.texts_to_sequences(virus_reads), maxlen = max_length, padding = 'post')
            predictions_virus = model.predict(virus_reads)
            idx_human = np.random.randint(len(predictinos_Human), size=len(predictinos_Human))
            total_predictions = np.concatenate((predictinos_Human[idx_human,:], predictions_virus), axis=0)
            
            ## look at all the proportions of reads per virus 
            proportion_Cornidovirineae = percentile_proportion('Cornidovirineae', total_predictions, 0.95)
            proportion_Influenza = percentile_proportion('Influenza', total_predictions, 0.95)
            proportion_SARS_CoV_2 = percentile_proportion('Sars_cov_2', total_predictions, 0.95)
            proportion_Metapneumovirus = percentile_proportion('Metapneumovirus', total_predictions, 0.95)
            proportion_Rhinovirus = percentile_proportion('Rhinovirus', total_predictions, 0.95)
            
            percentage_results[str(percentage)+'_'+virus] = [proportion_Cornidovirineae, 
                                                             proportion_Influenza, 
                                                             proportion_SARS_CoV_2,
                                                             proportion_Metapneumovirus,
                                                             proportion_Rhinovirus
                                                            ]
            print([proportion_Cornidovirineae, 
                 proportion_Influenza, 
                 proportion_SARS_CoV_2,
                 proportion_Metapneumovirus,
                 proportion_Rhinovirus
                ])
    
    
    # make a heatmap per percentage
    
    # Influenza
    influenza_percentages = [
                percentage_results['0.5_influenza'],
                percentage_results['1_influenza'],
                percentage_results['2.5_influenza'],
                percentage_results['5_influenza'],
                percentage_results['10_influenza']
                ]
    Corno = []
    Rhi = []
    Sars = []
    Influ = []
    Metap = []
    for percentage in enumerate(influenza_percentages):
        Corno.append(percentage[1][0])
        Rhi.append(percentage[1][4])
        Sars.append(percentage[1][2])
        Influ.append(percentage[1][1])
        Metap.append(percentage[1][3])
    
    f, ax = plt.subplots(figsize=(13,9))
    plt.title('miseq experiments Influenza + Human different proportions')
    sns.lineplot([0.5, 1, 2.5, 5, 10], Corno, label='Cornidovirineae')
    sns.lineplot([0.5, 1, 2.5, 5, 10], Rhi, label= 'Rhinovirus')
    sns.lineplot([0.5, 1, 2.5, 5, 10], Sars, label = 'Sars_cov_2')
    sns.lineplot([0.5, 1, 2.5, 5, 10], Metap, label = 'Metapneumovirus')
    sns.lineplot([0.5, 1, 2.5, 5, 10],Influ, label = 'Influenza')
    plt.xticks([0.5, 1, 2.5, 5, 10])
    plt.yticks([0.5, 1, 2.5, 5, 10])
    plt.savefig('/media/labuser/Data/COVID-19_classifier/pacific/results/9-mers/FPR_0.95_Novaseq_experiments_proportions_influenza.pdf',
                format='pdf',
                dpi=1200,
                bbox_inches='tight', pad_inches=0)
    
    # Cornidovirineae
    Cornidovirineae_percentages = [
                percentage_results['0.5_Cornidovirineae'],
                percentage_results['1_Cornidovirineae'],
                percentage_results['2.5_Cornidovirineae'],
                percentage_results['5_Cornidovirineae'],
                percentage_results['10_Cornidovirineae']
                ]
    Corno = []
    Rhi = []
    Sars = []
    Influ = []
    Metap = []
    for percentage in enumerate(Cornidovirineae_percentages):
        Corno.append(percentage[1][0])
        Rhi.append(percentage[1][4])
        Sars.append(percentage[1][2])
        Influ.append(percentage[1][1])
        Metap.append(percentage[1][3])
    
    f, ax = plt.subplots(figsize=(13,9))
    plt.title('miseq experiments Cornidovirineae + Human different proportions')
    sns.lineplot([0.5, 1, 2.5, 5, 10], Corno, label='Cornidovirineae')
    sns.lineplot([0.5, 1, 2.5, 5, 10], Rhi, label= 'Rhinovirus')
    sns.lineplot([0.5, 1, 2.5, 5, 10], Sars, label = 'Sars_cov_2')
    sns.lineplot([0.5, 1, 2.5, 5, 10], Metap, label = 'Metapneumovirus')
    sns.lineplot([0.5, 1, 2.5, 5, 10],Influ, label = 'Influenza')
    plt.xticks([0.5, 1, 2.5, 5, 10])
    plt.yticks([0.5, 1, 2.5, 5, 10])
    plt.savefig('/media/labuser/Data/COVID-19_classifier/pacific/results/9-mers/FPR_0.95_Novaseq_experiments_proportions_Cornidovirineae.pdf',
                format='pdf',
                dpi=1200,
                bbox_inches='tight', pad_inches=0)
    
    
    # Sars_cov_2
    Sars_cov_2_percentages = [
                percentage_results['0.5_Sars_cov_2'],
                percentage_results['1_Sars_cov_2'],
                percentage_results['2.5_Sars_cov_2'],
                percentage_results['5_Sars_cov_2'],
                percentage_results['10_Sars_cov_2']
                ]
    Corno = []
    Rhi = []
    Sars = []
    Influ = []
    Metap = []
    for percentage in enumerate(Sars_cov_2_percentages):
        Corno.append(percentage[1][0])
        Rhi.append(percentage[1][4])
        Sars.append(percentage[1][2])
        Influ.append(percentage[1][1])
        Metap.append(percentage[1][3])
    
    f, ax = plt.subplots(figsize=(13,9))
    plt.title('miseq experiments Sars_cov_2 + Human different proportions')
    sns.lineplot([0.5, 1, 2.5, 5, 10], Corno, label='Cornidovirineae')
    sns.lineplot([0.5, 1, 2.5, 5, 10], Rhi, label= 'Rhinovirus')
    sns.lineplot([0.5, 1, 2.5, 5, 10], Sars, label = 'Sars_cov_2')
    sns.lineplot([0.5, 1, 2.5, 5, 10], Metap, label = 'Metapneumovirus')
    sns.lineplot([0.5, 1, 2.5, 5, 10],Influ, label = 'Influenza')
    plt.xticks([0.5, 1, 2.5, 5, 10])
    plt.yticks([0.5, 1, 2.5, 5, 10])
    plt.savefig('/media/labuser/Data/COVID-19_classifier/pacific/results/9-mers/FPR_0.95_Novaseq_experiments_proportions_Sars_cov_2.pdf',
                format='pdf',
                dpi=1200,
                bbox_inches='tight', pad_inches=0)
    
    # Metapneumovirus
    Rhinovirus_percentages = [
                percentage_results['0.5_Rhinovirus'],
                percentage_results['1_Rhinovirus'],
                percentage_results['2.5_Rhinovirus'],
                percentage_results['5_Rhinovirus'],
                percentage_results['10_Rhinovirus']
                ]
    Corno = []
    Rhi = []
    Sars = []
    Influ = []
    Metap = []
    for percentage in enumerate(Rhinovirus_percentages):
        Corno.append(percentage[1][0])
        Rhi.append(percentage[1][4])
        Sars.append(percentage[1][2])
        Influ.append(percentage[1][1])
        Metap.append(percentage[1][3])
    
    f, ax = plt.subplots(figsize=(13,9))
    plt.title('miseq experiments Rhinovirus + Human different proportions')
    sns.lineplot([0.5, 1, 2.5, 5, 10], Corno, label='Cornidovirineae')
    sns.lineplot([0.5, 1, 2.5, 5, 10], Rhi, label= 'Rhinovirus')
    sns.lineplot([0.5, 1, 2.5, 5, 10], Sars, label = 'Sars_cov_2')
    sns.lineplot([0.5, 1, 2.5, 5, 10], Metap, label = 'Metapneumovirus')
    sns.lineplot([0.5, 1, 2.5, 5, 10],Influ, label = 'Influenza')
    plt.xticks([0.5, 1, 2.5, 5, 10])
    plt.yticks([0.5, 1, 2.5, 5, 10])
    plt.savefig('/media/labuser/Data/COVID-19_classifier/pacific/results/9-mers/FPR_0.95_Novaseq_experiments_proportions_Rhinovirus.pdf',
                format='pdf',
                dpi=1200,
                bbox_inches='tight', pad_inches=0)
    
    
    # Rhinovirus
    Metapneumovirus_percentages = [
                percentage_results['0.5_Metapneumovirus'],
                percentage_results['1_Metapneumovirus'],
                percentage_results['2.5_Metapneumovirus'],
                percentage_results['5_Metapneumovirus'],
                percentage_results['10_Metapneumovirus']
                ]
    Corno = []
    Rhi = []
    Sars = []
    Influ = []
    Metap = []
    for percentage in enumerate(Metapneumovirus_percentages):
        Corno.append(percentage[1][0])
        Rhi.append(percentage[1][4])
        Sars.append(percentage[1][2])
        Influ.append(percentage[1][1])
        Metap.append(percentage[1][3])
    
    f, ax = plt.subplots(figsize=(13,9))
    plt.title('miseq experiments Metapneumovirus + Human different proportions')
    sns.lineplot([0.5, 1, 2.5, 5, 10], Corno, label='Cornidovirineae')
    sns.lineplot([0.5, 1, 2.5, 5, 10], Rhi, label= 'Rhinovirus')
    sns.lineplot([0.5, 1, 2.5, 5, 10], Sars, label = 'Sars_cov_2')
    sns.lineplot([0.5, 1, 2.5, 5, 10], Metap, label = 'Metapneumovirus')
    sns.lineplot([0.5, 1, 2.5, 5, 10],Influ, label = 'Influenza')
    plt.xticks([0.5, 1, 2.5, 5, 10])
    plt.yticks([0.5, 1, 2.5, 5, 10])
    plt.savefig('/media/labuser/Data/COVID-19_classifier/pacific/results/9-mers/FPR_0.95_Novaseq_experiments_proportions_Metapneumovirus.pdf',
                format='pdf',
                dpi=1200,
                bbox_inches='tight', pad_inches=0)
    
    
    cutoffs = {'Cornidovirineae': 0.250501,
               'Metapneumovirus': 0.0335,
               'Rhinovirus': 0.38288,
               'Sars-cov-2': 0.32697,
               'Influenza': 0.0751592
               }
    
    df_percentages = pd.DataFrame(percentage_results)
    for name in classes:
        name_cols = [col for col in df_percentages.columns if name in col]
        df_temp = df_percentages[name_cols]
        df_temp.columns = ['10%', '5%', '2.5%', '1%', '0,5%']
        df_temp.index = ['Cornidovirineae','Influenza','Sars-cov-2','Metapneumovirus','Rhinovirus']
        plt.figure(figsize=(13,9))
        plt.title('Percentage of '+name+' in the sample')
        ax = sns.heatmap(df_temp,
                         annot=True,
                         cmap = sns.color_palette("Blues"))

        ax.tick_params(labelsize=15)
        for text in ax.texts:
            text.set_size(14)
            if text.get_position()[1] == 0.5:
                if float(text.get_text()) > 0.25050:
                    text.set_weight('bold')
                    text.set_size(18)
            if text.get_position()[1] == 1.5:
                if float(text.get_text()) > 0.07515:
                    text.set_weight('bold')
                    text.set_size(18)
            if text.get_position()[1] == 2.5:
                if float(text.get_text()) > 0.32697:
                    text.set_weight('bold')
                    text.set_size(18)
            if text.get_position()[1] == 3.5:
                if float(text.get_text()) > 0.0335:
                    text.set_weight('bold')
                    text.set_size(18)
            if text.get_position()[1] == 4.5:
                if float(text.get_text()) > 0.38288:
                    text.set_weight('bold')
                    text.set_size(18)
        ax.xaxis.tick_top() # x axis on top
        plt.savefig('/media/labuser/Data/COVID-19_classifier/pacific/results/9-mers/FPR_0.95_Novaseq_experiments_heatmap_'+name+'.pdf',
                format='pdf',
                dpi=1200,
                bbox_inches='tight', pad_inches=0)
    
    # correlations between real and predicted
    
    
    real = [10, 5, 2.5, 1, 0.5, 0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0,
            0,0,0,0,0, 10, 5, 2.5, 1, 0.5 ,0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0,
            0,0,0,0,0,  0,0,0,0,0, 10, 5, 2.5, 1, 0.5, 0,0,0,0,0, 0,0,0,0,0,
            0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0, 10, 5, 2.5, 1, 0.5, 0,0,0,0,0,
            0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0, 10, 5, 2.5, 1, 0.5 ]
    
    predicted = [10, 5, 2.5, 1, 0.5, 0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0, 0.0096,0.0086,0.0078,0.0075,0.009,
                 0.0013,0.0023,0.0012,0.0012,0.0008, 9.7,4.7,2.3,0.91,0.56, 0,0,0,0,0, 0,0,0,0,0,  0.0089, 0.0088,0.0092, 0.0087,0.0068,
                 0.0013,0.0017,0.0002,0.0012,0.0006, 0,0,0,0,0, 10, 5, 2.5, 1, 0.5, 0,0,0,0,0,  0.0067, 0.0087,0.01, 0.0085,0.0084,
                 0.00091,0.0015,0.00014,0.00099,0.00014, 0,0,0,0,0, 0,0,0,0,0,10, 5, 2.5, 1, 0.5, 0.0075, 0.013,0.0094, 0.011, 0.01,
                 0.0013,0.00075,0.00059,0.0014,0.0006, 0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0, 9.9, 5, 2.5, 1, 0.51]

    print(stats.pearsonr(real, predicted))
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    