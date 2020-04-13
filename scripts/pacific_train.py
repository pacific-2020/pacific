#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Apr  5 21:42:12 2020

Grab reads from fasta files and train the model

@author: labuser
"""

from Bio import SeqIO
import random
import os

from sklearn.preprocessing import LabelBinarizer
import numpy as np
from keras.preprocessing.text import Tokenizer

from sklearn.model_selection import train_test_split
from keras.preprocessing.sequence import pad_sequences

from keras.models import Sequential
from keras.layers import Embedding, LSTM, Dense, Bidirectional, Conv1D, CuDNNLSTM
from keras.layers import Dropout, Activation, MaxPooling1D
import tensorflow as tf

from numpy.random import seed
from tensorflow import set_random_seed
import pickle
import time
from sklearn.utils import shuffle
from keras.models import load_model


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


def process_reads(sequences, number_reads, length, kmer):
    '''
    '''
    r_reads = []
    for i in enumerate(sequences):
        # check the reads does not contain weird characters
        if all(c in 'AGCT' for c in i[1].upper()):
            r_reads.append(' '.join(i[1][x:x+kmer].upper() for x in range(len(i[1]) - kmer + 1)))
    return r_reads


def main(directory, number_reads, size_lenght, k_mer_size):
    '''
    '''
    files = os.listdir(directory)
    reads = []
    for file in files:
        all_transcripts = prepare_read(directory+'/'+file)
        reads += process_reads(all_transcripts, 
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
    
    seed_value = 42
    random.seed(seed_value)# 3. Set `numpy` pseudo-random generator at a fixed value
    np.random.seed(seed_value)# 4. Set `tensorflow` pseudo-random generator at a fixed value
    tf.set_random_seed(seed_value)# 5. For layers that introduce randomness like dropout, make sure to set seed values 
       
    config = tf.ConfigProto()
    config.gpu_options.allow_growth = True
    sess = tf.Session(config=config)
    
    # Read lenght
    read_lenght = 150
    
    # make synthetic reads
    Cornidovirineae_reads = main('/media/labuser/Data/COVID-19_classifier/pacific/data/synthetic_reads/Cornidovirineae',
                                         1, read_lenght, 4)
    Influenza_reads = main('/media/labuser/Data/COVID-19_classifier/pacific/data/synthetic_reads/Influenza',
                                         1, read_lenght, 4)
    Metapneumovirus_reads = main('/media/labuser/Data/COVID-19_classifier/pacific/data/synthetic_reads/Metapneumovirus',
                                         1, read_lenght, 4)
    Rhinovirus_reads = main('/media/labuser/Data/COVID-19_classifier/pacific/data/synthetic_reads/Rhinovirus',
                                         1, read_lenght, 4)
    Sars_cov_2_reads = main('/media/labuser/Data/COVID-19_classifier/pacific/data/synthetic_reads/Sars_cov_2',
                                         1, read_lenght, 4)
    Human = main('/media/labuser/Data/COVID-19_classifier/pacific/data/synthetic_reads/Human',
                                         1, read_lenght, 4) 

    total_sequences =  Cornidovirineae_reads + \
                       Influenza_reads +\
                       Metapneumovirus_reads +\
                       Rhinovirus_reads +\
                       Sars_cov_2_reads +\
                       Human
    
    labels_to_fit = ['Cornidovirineae','Influenza',"Metapneumovirus","Rhinovirus","Sars_cov_2", 'Human']
    label_maker = LabelBinarizer()
    transfomed_label = label_maker.fit(labels_to_fit)
    
    # save label_maker
    with open('/media/labuser/Data/COVID-19_classifier/pacific/model/label_maker.pickle', 'wb') as handle:
        pickle.dump(label_maker, handle, protocol=pickle.HIGHEST_PROTOCOL)
    
    labels = list(np.repeat('Cornidovirineae',len(Cornidovirineae_reads))) + \
             list(np.repeat('Influenza',len(Influenza_reads))) + \
             list(np.repeat('Metapneumovirus',len(Metapneumovirus_reads))) + \
             list(np.repeat('Rhinovirus',len(Rhinovirus_reads))) + \
             list(np.repeat('Sars_cov_2',len(Sars_cov_2_reads))) + \
             list(np.repeat('Human',len(Human)))
             
    labels_proces = label_maker.transform(labels)
    
    # Tokenize the vocabulary
    tokenizer = Tokenizer()
    tokenizer.fit_on_texts(total_sequences)
    sequences_preproces = tokenizer.texts_to_sequences(total_sequences)
    
    max_features = len(tokenizer.word_index)+1
    
    max_length = max([len(s.split()) for s in total_sequences])
    # pad sequences
    sequences_preproces = pad_sequences(sequences_preproces, maxlen = max_length, padding = 'post')
    
    sequences_preproces, labels_proces = shuffle(sequences_preproces, labels_proces)
    
    with open('/media/labuser/Data/COVID-19_classifier/pacific/model/tokenizer.pickle', 'wb') as handle:
        pickle.dump(tokenizer, handle, protocol=pickle.HIGHEST_PROTOCOL)
    '''
    # loading
    with open('/media/labuser/Data/COVID-19_classifier/pacific/model/tokenizer.pickle', 'rb') as handle:
        tokenizer = pickle.load(handle)
    '''
    
    # Netweork parameter s
    
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
    model.add(Bidirectional(CuDNNLSTM(lstm_output_size)))
    model.add(Dense(50))
    model.add(Dense(6))

    model.add(Activation('softmax'))
    
    model.compile(loss='categorical_crossentropy',
                  optimizer='adam',
                  metrics=['binary_accuracy', 
                           'categorical_accuracy',
                           ])
    model.summary()
    
    # training time
    start = time.clock()
    
    print('Train...')
    
    for epoch in range(epochs):
        print("epoch %d" %epoch)
        #train in batches of 100k sequences
        for i in range(0, len(sequences_preproces), 200000):
            start, end = i, i+200000
            if end > len(sequences_preproces):
                end = len(sequences_preproces)
            print('chunk: ',start, end)
            training_batch = sequences_preproces[start:end]
            labels_batch = labels_proces[start:end] 
            X_train,X_test,y_train,y_test = train_test_split(training_batch, 
                                                             labels_batch,
                                                             test_size=0.10, 
                                                             random_state=42)
            model.fit(X_train, y_train,
                      batch_size=batch_size,
                      epochs=1,
                      validation_data=(X_test, y_test))
    
    score, binary_acc, categorical_acc = model.evaluate(X_test, y_test, batch_size=batch_size)
    print('Traning time:', time.clock() - start)
    print('Test accuracy:', binary_acc, categorical_acc, score)
    
    # save keras model
    model.save("/media/labuser/Data/COVID-19_classifier/pacific/model/pacific.h5")
    print("Saved model to disk")
    















