#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 28 08:52:20 2020

In this script I am going to convert genomes into tables with k-mers and abundances and perform UMAP on them

@author: labuser
"""

from Bio import SeqIO
import itertools
import os
import pandas as pd
import umap
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np


def prepare_genome(trancriptome):
    '''
    function will take tranciprtome and make reads
    '''
    fasta_sequences = SeqIO.parse(open(trancriptome),'fasta')
    sequences = []
    for fasta in fasta_sequences:
        name, sequence = fasta.id, str(fasta.seq)
        sequences.append(sequence)
        name = name.split('\t')[0]
    return name, sequences

def make_k_mers(sequence, k_mer_lenght, dictionary):
    '''
    '''
    for i in range(0, len(sequence)-k_mer_lenght):
        kmer = sequence[i:i+k_mer_lenght]
        if all(c in 'AGCT' for c in kmer.upper()):
            dictionary[kmer] +=1 
    return dictionary

def main(genome, dictionary, k_mer):
    '''
    '''
    name, sequence = prepare_genome(genome)
    dictionary_done = make_k_mers(sequence[0], k_mer, dictionary)
    return dictionary_done, name


def file_to_df(path, max_files, dictionary, virus_group, kmer_lenght):
    '''
    '''
    genomes_k_mer = pd.DataFrame()
    virus_genomes  = os.listdir(path)
    virus_genomes.pop(virus_genomes.index(''.join([i for i in virus_genomes if i[-5:] =='fasta'])))
    
    for genome in virus_genomes:
        dictionary_temp, name = main(path+'/'+genome, dictionary.copy(), kmer_lenght)
        dictionary_temp['virus_id'] = name
        dictionary_temp['Virus_group'] = virus_group
        dictionary_temp = pd.DataFrame(dictionary_temp, index=[0])
        if len(genomes_k_mer) == 0:
            genomes_k_mer = dictionary_temp
        else:
            genomes_k_mer = pd.concat([genomes_k_mer, dictionary_temp], ignore_index=True)
        if genomes_k_mer.shape[0] >= max_files:
            return genomes_k_mer
    
    return genomes_k_mer



if __name__ == '__main__':

    for KMER in range(4,9):
    
        ## make an empty dictionary with all possible 4-mers
        kmers_keys = itertools.product(['A','G','T','C'], repeat=KMER) 
        dictionary_kmers_empty = {}
        for i in kmers_keys:
            i = ''.join(i)
            dictionary_kmers_empty[i] = 0
        
        counter = {}
        path_Metapneumovirus = '/media/labuser/Data/COVID-19_classifier/pacific/data/training_genomes/Metapneumovirus/'
        path_Cornidovirineae = '/media/labuser/Data/COVID-19_classifier/pacific/data/training_genomes/Cornidovirineae/'
        path_rhinovirus = '/media/labuser/Data/COVID-19_classifier/pacific/data/training_genomes/rhinovirus/'
        path_Influenza = '/media/labuser/Data/COVID-19_classifier/pacific/data/training_genomes/Influenza/'
        path_Sars_cov_2 = '/media/labuser/Data/COVID-19_classifier/pacific/data/training_genomes/Sars_cov_2/'
        
        df_Metapneumovirus = file_to_df('/media/labuser/Data/COVID-19_classifier/pacific/data/training_genomes/Metapneumovirus/',
                                        200, 
                                        dictionary_kmers_empty.copy(),
                                        0,
                                        KMER)
        
        df_Cornidovirineae = file_to_df('/media/labuser/Data/COVID-19_classifier/pacific/data/training_genomes/Cornidovirineae/',
                                        200,
                                        dictionary_kmers_empty.copy(),
                                        1,
                                        KMER)

        df_rhinovirus = file_to_df('/media/labuser/Data/COVID-19_classifier/pacific/data/training_genomes/rhinovirus/',
                                   200,
                                   dictionary_kmers_empty.copy(),
                                   2,
                                   KMER)

        df_Influenza = file_to_df('/media/labuser/Data/COVID-19_classifier/pacific/data/training_genomes/Influenza/',
                                  200,
                                  dictionary_kmers_empty.copy(),
                                  3,
                                  KMER)

        df_Sars_cov_2 = file_to_df('/media/labuser/Data/COVID-19_classifier/pacific/data/training_genomes/Sars_cov_2/',
                                   200,
                                   dictionary_kmers_empty.copy(),
                                   4,
                                   KMER)
        
        
        #genomes_numpy = genomes_k_mer.iloc[:, :-1].to_numpy()
        
        df_genomes = pd.concat((df_Metapneumovirus,df_Cornidovirineae, df_rhinovirus, df_Influenza, df_Sars_cov_2))
        
        numpy_genomes = df_genomes.iloc[:,:-2].to_numpy() 
        
        reducer = umap.UMAP()
        embedding = reducer.fit_transform(numpy_genomes)

        f, ax = plt.subplots(figsize=(13,9))
        plt.title('Genome embeddings with UMAP k-mer: '+str(KMER))
        sns.scatterplot(x=embedding[:5,0], y=embedding[:5,1], label='Metapneumovirus', s=100)
        sns.scatterplot(x=embedding[5:17,0], y=embedding[5:17,1], label='Cornidovirineae',s=100)
        sns.scatterplot(x=embedding[17:147,0], y=embedding[17:147,1], label='rhinovirus',s=100)
        sns.scatterplot(x=embedding[147:275,0], y=embedding[147:275,1], label='Influenza',s=100)
        sns.scatterplot(x=embedding[275:362,0], y=embedding[275:362,1], label='Sars_cov_2',s=100)
        plt.savefig('/media/labuser/Data/COVID-19_classifier/pacific/results/genome_embeddings_kmer'+str(KMER)+'.pdf',
                    format='pdf',
                    dpi=1200,
                    bbox_inches='tight',
                    pad_inches=0)

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    

