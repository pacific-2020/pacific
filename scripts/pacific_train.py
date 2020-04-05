#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  1 15:45:58 2020

RNN classifier for COVID-19 reads

@author: labuser
"""
import os
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
import random 

def get_reads(samfile, kmer):
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
                if len(line) > 50 and 'N' not in line:
                    #make kmer
                    reads.append(' '.join(line[x:x+kmer].upper() for x in range(len(line) - kmer + 1)))
            if counter == 4:
                counter = 0
    return reads


def read_preprocess(files, number, kmer_size):
    '''
    This function take a list of reads and preprocess a number of reads
    with tokenizer
    '''
    reads = []
    for file in enumerate(files):
        reads += get_reads(path+file[1], kmer_size)
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
    
    
    # short reads
    path = '/media/labuser/Data/COVID-19_classifier/pacific/data/illumina/processed/'
    files = os.listdir(path)
    files_covid = [i for i in files if i.endswith('filtered_covid.sam.fastq')]
    files_non_covid = [i for i in files if i.endswith('non_covid.sam.fastq')]
    
    # Lets try 20k reads for training each
    covid = read_preprocess(files_covid, 20000, 4)
    non_covid = read_preprocess(files_non_covid, 20000, 4)
    
    total_sequences =  covid + non_covid
    labels = list(np.ones(len(covid))) + list(np.zeros(len(non_covid)))
    
    # Tokenize the vocabulary
    tokenizer = Tokenizer()
    tokenizer.fit_on_texts(total_sequences)
    sequences_preproces = tokenizer.texts_to_sequences(total_sequences)
    
    max_features = len(tokenizer.word_index) 
    
    max_length = max([len(s.split()) for s in total_sequences])
    sequences_preproces = pad_sequences(sequences_preproces, maxlen = max_length, padding = 'post')
    
    #split dataset into training and testing
    X_train,X_test,y_train,y_test = train_test_split(sequences_preproces, 
                                                     labels,
                                                     test_size=0.10, 
                                                     random_state=42)
    
    X_train = np.asarray(X_train)
    X_test = np.asarray(X_test)
    y_test = np.asarray(y_test)
    y_train = np.asarray(y_train)
    
    ## Define model parameters
    
    # Convolution
    kernel_size = 3
    filters = 64
    pool_size = 4
    
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
    model.add(LSTM(lstm_output_size))
    model.add(Dense(1))
    model.add(Activation('sigmoid'))
    
    model.compile(loss='binary_crossentropy',
                  optimizer='adam',
                  metrics=['accuracy'])
    
    model.summary()
    
    print('Train...')
    model.fit(X_train, y_train,
              batch_size=batch_size,
              epochs=epochs,
              validation_data=(X_test, y_test))
    score, acc = model.evaluate(X_test, y_test, batch_size=batch_size)
    print('Test accuracy:', acc)
    
    
    
    
    
    
    
    
    
    
    






























