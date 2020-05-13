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
from keras.layers import Embedding, LSTM, Dense, Bidirectional, Conv1D, CuDNNLSTM, GRU
from keras.layers import Dropout, Activation, MaxPooling1D, Flatten
import tensorflow as tf

from numpy.random import seed
from tensorflow import set_random_seed
import pickle
import time
from sklearn.utils import shuffle
from keras.models import load_model
import matplotlib.pyplot as plt
import seaborn as sns

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


def process_reads(sequences, length, kmer):
    '''
    '''
    r_reads = []
    for i in enumerate(sequences):
        # check the reads does not contain weird characters
        if all(c in 'AGCT' for c in i[1].upper()):
            r_reads.append(' '.join(i[1][x:x+kmer].upper() for x in range(len(i[1]) - kmer + 1)))
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
    
    
    pacific_9mers = 'pacific_9mers'
    kmers = 9
    
    # Read lenght
    read_lenght = 150
    
    # get synthetic reads
    Cornidovirineae_reads = main('/media/labuser/Data/COVID-19_classifier/pacific/data/synthetictrainingdata_group/Cornidovirineae',
                                 read_lenght, kmers)
    Influenza_reads = main('/media/labuser/Data/COVID-19_classifier/pacific/data/synthetictrainingdata_group/Influenza',
                           read_lenght, kmers)
    Metapneumovirus_reads = main('/media/labuser/Data/COVID-19_classifier/pacific/data/synthetictrainingdata_group/Metapneumovirus',
                                 read_lenght, kmers)
    Rhinovirus_reads = main('/media/labuser/Data/COVID-19_classifier/pacific/data/synthetictrainingdata_group/Rhinovirus',
                            read_lenght, kmers)
    Sars_cov_2_reads = main('/media/labuser/Data/COVID-19_classifier/pacific/data/synthetictrainingdata_group/Sars_Cov-2/',
                            read_lenght, kmers)
    Human = main('/media/labuser/Data/COVID-19_classifier/pacific/data/synthetictrainingdata_group/Human/',
                 read_lenght, kmers)
    
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
    with open('/media/labuser/Data/COVID-19_classifier/pacific/model/label_maker.01.'+pacific_9mers+'.pickle', 'wb') as handle:
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
    
    #np.save('/media/labuser/Data/COVID-19_classifier/pacific/data/training_objects/sequences_'+str(kmers)+'.npy', sequences_preproces)
    #np.save('/media/labuser/Data/COVID-19_classifier/pacific/data/training_objects/labels_'+str(kmers)+'.npy', labels_proces)

    #with open('/media/labuser/Data/COVID-19_classifier/pacific/model/tokenizer.01.'+pacific_9mers+'.pickle', 'wb') as handle:
    #    pickle.dump(tokenizer, handle, protocol=pickle.HIGHEST_PROTOCOL)
    
    # loading tokenizer and label_maker
    with open('/media/labuser/Data/COVID-19_classifier/pacific/model/tokenizer.01.'+pacific_9mers+'.pickle', 'rb') as handle:
        tokenizer = pickle.load(handle)
    with open('/media/labuser/Data/COVID-19_classifier/pacific/model/label_maker.01.'+pacific_9mers+'.pickle', 'rb') as handle:
        label_maker = pickle.load(handle)
    
    # loading data and labels
    sequences_preproces = np.load('/media/labuser/Data/COVID-19_classifier/pacific/data/training_objects/sequences_'+str(kmers)+'.npy') 
    labels_proces = np.load('/media/labuser/Data/COVID-19_classifier/pacific/data/training_objects/labels_'+str(kmers)+'.npy')
    
    sequences_preproces, labels_proces = shuffle(sequences_preproces, labels_proces)

    # Convolution
    kernel_size = 3
    filters = 128
    pool_size = 3
    
    # LSTM
    lstm_output_size = 70
    
    # Training
    batch_size = 30
    epochs = 1
    
    
    # Define the model the model
    model = Sequential()
    model.add(Embedding(max_features, 100, input_length=sequences_preproces.shape[1]))
    model.add(Dropout(0.20))
    model.add(Conv1D(filters,
                     kernel_size,
                     padding='same',
                     activation='relu',
                     strides=1))
    model.add(MaxPooling1D(pool_size=pool_size))
    model.add(Dropout(0.1))
    model.add(Bidirectional(CuDNNLSTM(lstm_output_size)))
    model.add(Dropout(0.1))
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
    start = time.process_time()
    
    histories = []
    print('Train...')
    for epoch in range(epochs):
        print("epoch %d" %epoch)
        #train in batches of 200k sequences
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
            chunk_history = model.fit(X_train, y_train,
                                      batch_size=batch_size,
                                      epochs=1,
                                      validation_data=(X_test, y_test))
            histories.append(chunk_history)
            if i == 15:
                break
            break
    
    end = time.process_time()
    print('Traning time:', start  - end)
    
    # save keras model
    model.save("/media/labuser/Data/COVID-19_classifier/pacific/model/pacific.01."+pacific_9mers+".h5")
    print("Saved model to disk")

    #### plot the accuracies and losses
    bi_acc = []
    cat_acc = []
    loss = []
    val_bi_acc = []
    val_cat_acc = []
    val_loss = []
    for i in histories:
        bi_acc.append(i.history['binary_accuracy'][0])
        cat_acc.append(i.history['categorical_accuracy'][0])
        loss.append(i.history['loss'][0])
        val_bi_acc.append(i.history['val_binary_accuracy'][0])
        val_cat_acc.append(i.history['val_categorical_accuracy'][0])
        val_loss.append(i.history['val_loss'][0])

    f, ax = plt.subplots( figsize=(13,9))
    sns.lineplot(x=np.arange(len(histories)), y=np.array(bi_acc), palette="tab10", linewidth=2.5, label='Binary accuracy')
    sns.lineplot(x=np.arange(len(histories)), y=np.array(cat_acc), palette="tab10", linewidth=2.5, label='Categorical accuracy')
    plt.ylabel('Accuracies')
    plt.ylabel('Accuracies')
    plt.savefig('/media/labuser/Data/COVID-19_classifier/pacific/results/trainning_accuracy_deep'+pacific_9mers+'.pdf',
                format='pdf',
                dpi=1200,
                bbox_inches='tight', pad_inches=0)
    
    f, ax = plt.subplots( figsize=(13,9))
    sns.lineplot(x=np.arange(len(histories)), y=np.array(loss), palette="tab10", linewidth=2.5, label='loss')
    plt.ylabel('Loss')
    plt.savefig('/media/labuser/Data/COVID-19_classifier/pacific/results/training_loss_deep'+pacific_9mers+'.pdf',
                format='pdf',
                dpi=1200,
                bbox_inches='tight', pad_inches=0)
    
    f, ax = plt.subplots( figsize=(13,9))
    sns.lineplot(x=np.arange(len(histories)), y=np.array(val_bi_acc), palette="tab10", linewidth=2.5, label='Validation binary accuracy')
    sns.lineplot(x=np.arange(len(histories)), y=np.array(val_cat_acc), palette="tab10", linewidth=2.5, label='Validation categorical accuracy')
    plt.ylabel('Percentage of predicted reads')
    plt.savefig('/media/labuser/Data/COVID-19_classifier/pacific/results/val_training_accuracy_deep'+pacific_9mers+'.pdf',
                format='pdf',
                dpi=1200,
                bbox_inches='tight', pad_inches=0)
    
    f, ax = plt.subplots( figsize=(13,9))
    sns.lineplot(x=np.arange(len(histories)), y=np.array(val_loss), palette="tab10", linewidth=2.5, label='Validation loss')
    plt.ylabel('Loss')
    plt.savefig('/media/labuser/Data/COVID-19_classifier/pacific/results/val_loss_deep'+pacific_9mers+'.pdf',
                format='pdf',
                dpi=1200,
                bbox_inches='tight', pad_inches=0)
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    




