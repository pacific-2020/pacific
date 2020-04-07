#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Apr  5 21:42:12 2020

Prepare the training with synthetic data that will be generated from transcriptomes

@author: labuser
"""

from Bio import SeqIO
import random
import os

from sklearn.preprocessing import LabelBinarizer
import numpy as np
from keras.preprocessing.text import Tokenizer, text_to_word_sequence

from sklearn.model_selection import train_test_split
from keras.preprocessing.sequence import pad_sequences

from keras.models import Sequential
from keras.layers import Embedding, LSTM, Dense, Bidirectional, Conv1D
from keras.layers import Dropout, Activation, MaxPooling1D
import tensorflow as tf

from numpy.random import seed
from tensorflow import set_random_seed


def prepare_read(trancriptome):
    '''
    function will take tranciprtome and make reads
    '''
    fasta_sequences = SeqIO.parse(open(trancriptome),'fasta')
    sequences = []
    for fasta in fasta_sequences:
        name, sequence = fasta.id, str(fasta.seq)
        sequences.append(sequence)
    return sequences


def random_reads(sequences, number_reads, length, kmer):
    '''
    '''
    r_reads = []
    for i in range(number_reads):
        transcript_number = random.randrange(len(sequences))
        try:
            begining = random.randrange(len(sequences[transcript_number])-length)
            sequence_complete = sequences[transcript_number][begining:begining+length]
            r_reads.append(' '.join(sequence_complete[x:x+kmer].upper() for x in range(len(sequence_complete) - kmer + 1)))
        except:
            continue
    return r_reads


def main(directory, number_reads, size_lenght, k_mer_size):
    '''
    '''
    files = os.listdir(directory)
    read_per_file =int(number_reads/len(files))
    reads = []
    for file in files:
        all_transcripts = prepare_read(directory+'/'+file)
        reads += random_reads(all_transcripts, 
                              read_per_file, 
                              size_lenght,
                              k_mer_size)
    return reads

### Illumina reads functions
    
def get_reads_illumina(samfile, kmer):
    '''
    This function take a filtered samfile and get reads
    '''
    reads = []
    with open(samfile, 'r') as file_in:
        counter = 0
        for line in file_in:
            counter +=1
            line = line.rstrip()
            if counter == 2:
                # set limit read lenght
                if len(line) >= 100 and 'N' not in line:
                    #make kmer
                    reads.append(' '.join(line[x:x+kmer].upper() for x in range(len(line) - kmer + 1)))
            if counter == 4:
                counter = 0
    return reads


def read_preprocess_illumina(files, number, kmer_size):
    '''
    This function take a list of reads and preprocess a number of reads
    with tokenizer
    '''
    reads = []
    for file in enumerate(files):
        reads += get_reads_illumina(path+file[1], kmer_size)
        if len(reads) > number:
            reads = reads[:number]
            break
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

##### functions for nanopore process reads
    

def get_reads_nanopore(samfile, kmer):
    '''
    This function take a filtered samfile and get reads
    '''
    reads = []
    with open(samfile, 'r') as file_in:
        counter = 0
        for line in file_in:
            counter +=1
            line = line.rstrip()
            if counter == 2:
                # set limit read lenght
                if len(line) >= 100 and 'N' not in line:
                    #make kmer
                    line = line[:150]
                    reads.append(' '.join(line[x:x+kmer].upper() for x in range(len(line) - kmer + 1)))
            if counter == 4:
                counter = 0
    return reads


def read_preprocess_nanopore(files, number, kmer_size):
    '''
    This function take a list of reads and preprocess a number of reads
    with tokenizer
    '''
    reads = []
    for file in enumerate(files):
        reads += get_reads_nanopore(path+file[1], kmer_size)
        if len(reads) > number:
            reads = reads[:number]
            break
    return reads


if __name__ == '__main__':
    

    seed_value = 42
    random.seed(seed_value)# 3. Set `numpy` pseudo-random generator at a fixed value
    np.random.seed(seed_value)# 4. Set `tensorflow` pseudo-random generator at a fixed value
    tf.set_random_seed(seed_value)# 5. For layers that introduce randomness like dropout, make sure to set seed values 
       
    config = tf.ConfigProto()
    config.gpu_options.allow_growth = True
    sess = tf.Session(config=config)
    
    # Read lenght
    read_lenght = 100
    
    # make synthetic reads
    Cornidovirineae_reads = main('/media/labuser/Data/COVID-19_classifier/pacific/data/virus_ncbi/Cornidovirineae',
                                         50000, read_lenght, 4)
    Influenza_reads = main('/media/labuser/Data/COVID-19_classifier/pacific/data/virus_ncbi/Influenza',
                                         50000, read_lenght, 4)
    Metapneumovirus_reads = main('/media/labuser/Data/COVID-19_classifier/pacific/data/virus_ncbi/Metapneumovirus',
                                         50000, read_lenght, 4)
    human_reads = main('/media/labuser/Data/COVID-19_classifier/pacific/data/human/',
                                         80000, read_lenght, 4)
    Sars_cov_2 = main('/media/labuser/Data/COVID-19_classifier/pacific/data/virus_ncbi/Sars_cov_2',
                                         30000, read_lenght, 4)
    
    # get real illumina reads from experiment 1 and 4
    path = '/media/labuser/Data/COVID-19_classifier/pacific/data/'+\
           'non_synthetic/illumina/processed/experiment_1_4/'
    
    files = os.listdir(path)
    files_covid = [i for i in files if i.endswith('filtered_covid.sam.fastq')]
    files_non_covid = [i for i in files if i.endswith('non_covid.sam.fastq')]
    
    covid = read_preprocess_illumina(files_covid, 15000, 4)
    non_covid = read_preprocess_illumina(files_non_covid, 15000, 4)
    
    human_reads += non_covid
    Sars_cov_2 += covid
    
    total_sequences =  Cornidovirineae_reads + \
                       Influenza_reads +\
                       Metapneumovirus_reads +\
                       human_reads +\
                       Sars_cov_2
    
    labels_to_fit = ['Cornidovirineae','Influenza',"Metapneumovirus","human_reads","Sars_cov_2"]
    label_maker = LabelBinarizer()
    transfomed_label = label_maker.fit(labels_to_fit)
    
    labels = list(np.repeat('Cornidovirineae',len(Cornidovirineae_reads))) + \
             list(np.repeat('Influenza',len(Influenza_reads))) + \
             list(np.repeat('Metapneumovirus',len(Metapneumovirus_reads))) + \
             list(np.repeat('human_reads',len(human_reads))) + \
             list(np.repeat('Sars_cov_2',len(Sars_cov_2)))
             
    labels_proces = label_maker.fit_transform(labels)
    
    # Tokenize the vocabulary
    tokenizer = Tokenizer()
    tokenizer.fit_on_texts(total_sequences)
    sequences_preproces = tokenizer.texts_to_sequences(total_sequences)
    
    max_features = len(tokenizer.word_index)+1
    
    max_length = max([len(s.split()) for s in total_sequences])
    sequences_preproces = pad_sequences(sequences_preproces, maxlen = max_length, padding = 'post')
    
    #split dataset into training and testing
    X_train,X_test,y_train,y_test = train_test_split(sequences_preproces, 
                                                     labels_proces,
                                                     test_size=0.10, 
                                                     random_state=42)
    
    # Convolution
    kernel_size = 3
    filters = 64
    pool_size = 2
    
    # LSTM
    lstm_output_size = 70
    
    # Training
    batch_size = 30
    epochs = 2
    
    
    # Define the model the model
    model = Sequential()
    #Input = Input(shape=(None,)))
    model.add(Embedding(max_features, 50))
    model.add(Dropout(0.20))
    model.add(Conv1D(filters,
                     kernel_size,
                     padding='valid',
                     activation='relu',
                     strides=1))
    model.add(MaxPooling1D(pool_size=pool_size))
    model.add(Bidirectional(LSTM(lstm_output_size)))
    model.add(Dense(5))
    model.add(Activation('softmax'))
    
    model.compile(loss='categorical_crossentropy',
                  optimizer='adam',
                  metrics=['binary_accuracy', 
                           'categorical_accuracy',
                           ])
    
    model.summary()
    
    print('Train...')
    model.fit(X_train, y_train,
              batch_size=batch_size,
              epochs=epochs,
              validation_data=(X_test, y_test))
    score, binary_acc, categorical_acc = model.evaluate(X_test, y_test, batch_size=batch_size)
    print('Test accuracy:', binary_acc, categorical_acc, score)
    
    ### seegin the acuracy for sars
    
    '''
    predictions = model.predict(X_test)
    
    q = np.argmax(predictions, axis=1)

    accuracy_validation = accuracy(np.array(labels_predict), predictions_argmax)
    '''
    
    #### Test with real data
    
    # short illumina reads
    path = '/media/labuser/Data/COVID-19_classifier/pacific/data/non_synthetic/illumina/processed/experiment_2_3/'
    files = os.listdir(path)
    files_covid = [i for i in files if i.endswith('filtered_covid.sam.fastq')]
    files_non_covid = [i for i in files if i.endswith('non_covid.sam.fastq')]
    
    # Lets try 20k reads for training each
    covid = read_preprocess_illumina(files_covid, 5000, 4)
    non_covid = read_preprocess_illumina(files_non_covid, 5000, 4)
    
    total_sequences_predict =  covid + non_covid
    
    
    # labes order are ['Cornidovirineae', 'Influenza', 'Metapneumovirus', 'Sars_cov_2', 'human_reads'])
    
    labels_predict = list(np.repeat(3, len(covid))) + \
                     list(np.repeat(4, len(non_covid)))

    sequences_predict = tokenizer.texts_to_sequences(total_sequences_predict)
    
    sequences_preproces_predict = pad_sequences(sequences_predict, maxlen = max_length, padding = 'post')

    predictions = model.predict(sequences_preproces_predict)
    
    # get rid of predictions where the best prediction is lower than 0.5
    labels_high_acc = []
    predictions_high_acc = []
    for i in enumerate(predictions):
        if max(i[1]) > 0.9:
            predictions_high_acc.append(i[1])
            labels_high_acc.append(labels_predict[i[0]])
            
    predictions_high_acc = np.array(predictions_high_acc)
    predictions_argmax = np.argmax(predictions_high_acc, axis=1)

    accuracy_validation = accuracy(np.array(labels_high_acc), predictions_argmax)
    
    results = []
    for i in enumerate(labels_high_acc):
        results.append((i[1], predictions_argmax[i[0]]))
        
    
    #### Test with real data Nanopore
    
    # short nanopore reads
    path = '/media/labuser/Data/COVID-19_classifier/pacific/data/non_synthetic/nanopore/processed/'
    files = os.listdir(path)
    files_covid = [i for i in files if i.endswith('filtered_covid.fastq')]
    files_non_covid = [i for i in files if i.endswith('non_covid.fastq')]
    
    # Lets try 20k reads for training each
    covid = read_preprocess_nanopore(files_covid, 5000, 4)
    non_covid = read_preprocess_nanopore(files_non_covid, 5000, 4)
    
    total_sequences_predict =  covid + non_covid
    
    # labes order are ['Cornidovirineae', 'Influenza', 'Metapneumovirus', 'Sars_cov_2', 'human_reads'])
    
    labels_predict = list(np.repeat(3, len(covid))) + \
                     list(np.repeat(4, len(non_covid)))

    sequences_predict = tokenizer.texts_to_sequences(total_sequences_predict)
    
    sequences_preproces_predict = pad_sequences(sequences_predict, maxlen = max_length, padding = 'post')

    predictions = model.predict(sequences_preproces_predict)
    
    # get rid of predictions where the best prediction is lower than 0.5
    labels_high_acc = []
    predictions_high_acc = []
    for i in enumerate(predictions):
        if max(i[1]) > 0.9:
            predictions_high_acc.append(i[1])
            labels_high_acc.append(labels_predict[i[0]])
            
    predictions_high_acc = np.array(predictions_high_acc)
    predictions_argmax = np.argmax(predictions_high_acc, axis=1)

    accuracy_validation = accuracy(np.array(labels_high_acc), predictions_argmax)
    
    results = []
    for i in enumerate(labels_high_acc):
        results.append((i[1], predictions_argmax[i[0]]))
        
    
    























