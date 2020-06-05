#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Apr  5 21:42:12 2020

This script will train PACIFIC and save the model, training and validation
plos in specified folders.

@author: Pablo Acera
"""

import argparse

parser = argparse.ArgumentParser(description=
                                 """
                                 This script train PACIFIC using fasta reads.
                                 The model, training and validation 
                                 plots will be safe in a specific foldere.
                                 """)

OPTIONAL = parser._action_groups.pop()
REQUIRED = parser.add_argument_group('required arguments')

#Inputs
REQUIRED.add_argument("--Coronaviridae_reads",
                      help="file path to folder containing Coronaviridae fasta files to train PACIFIC",
                      required=True)

REQUIRED.add_argument("--Influenza_reads",
                      help="file path to folder containing Influenza fasta files to train PACIFIC",
                      required=True)

REQUIRED.add_argument("--Metapneumovirus_reads",
                      help="file path to folder containing Metapneumovirus fasta files to train PACIFIC",
                      required=True)

REQUIRED.add_argument("--Rhinovirus_reads",
                      help="file path to folder containing Rhinovirus fasta files to train PACIFIC",
                      required=True)

REQUIRED.add_argument("--Sars_cov_2_reads",
                      help="file path to folder containing SARS-CoV-2 fasta files to train PACIFIC",
                      required=True)

REQUIRED.add_argument("--Human_reads",
                      help="file path to folder containing Human fasta files to train PACIFIC",
                      required=True)


#arguments
OPTIONAL.add_argument("--out_folder",
                      help='path to the output folder',
                      default="./")

OPTIONAL.add_argument("--k_mers",
                      help='K-mer number use to train the model',
                      default=9,
                      type=int)

OPTIONAL.add_argument("--model_name",
                      help='Name used to save the model',
                      default="PACIFIC")

OPTIONAL.add_argument("--stop_chunk",
                      help='Chunk number to stop the training',
                      default=15,
                      type=int)

OPTIONAL.add_argument("--GPU",
                      help='If True  PACIFIC will be train using CuDNNLSTM',
                      default=False,
                      )

OPTIONAL.add_argument("--file_type",
                      help='fasta or fastq training files format (all files should have same format)',
                      default='fasta',
                      )


parser._action_groups.append(OPTIONAL)

ARGS = parser.parse_args()

# Inputs
CORONAVIRIDAE_READS = ARGS.Coronaviridae_reads
INFLUENZA_READS = ARGS.Influenza_reads
METAPMEUMOVIRUS_READS = ARGS.Metapneumovirus_reads
RHINOVIRUS_READS = ARGS.Rhinovirus_reads
SARS_COV_2_READS = ARGS.Sars_cov_2_reads
HUMAN_READS = ARGS.Human_reads

# Arguments
OUT_FOLDER = ARGS.out_folder
K_MERS = ARGS.k_mers
MODEL_NAME = ARGS.model_name
STOP_CHUNK = ARGS.stop_chunk
GPU = ARGS.GPU
FILE_TYPE = ARGS.file_type


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
import matplotlib.pyplot as plt
import seaborn as sns


def prepare_read(trancriptome, file_type):
    '''
    function will take tranciprtome and make reads
    '''
    fasta_sequences = SeqIO.parse(open(trancriptome),file_type)
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


def main(directory, size_lenght, k_mer_size, file_type):
    '''
    '''
    files = os.listdir(directory)
    reads = []
    for file in files:
        all_transcripts = prepare_read(directory+'/'+file, file_type)
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
    '''
    config = tf.ConfigProto()
    config.gpu_options.allow_growth = True
    sess = tf.Session(config=config)
    '''
    kmers = K_MERS
    
    # Create output folder if it does not exist
    if os.path.isdir(OUT_FOLDER) is False:
        print('Creating output folder '+OUT_FOLDER)
        os.mkdir(OUT_FOLDER)

    
    # Read lenght
    read_lenght = 150
    
    # get synthetic reads
    print('Loading Coronaviridae reads')
    
    Coronaviridae_reads = main(CORONAVIRIDAE_READS,
                               read_lenght, 
                               kmers,
                               FILE_TYPE,
                              )
    
    print('Loading Influenza reads')
    
    Influenza_reads = main(INFLUENZA_READS,
                           read_lenght,
                           kmers,
                           FILE_TYPE
                          )
    
    print('Loading Metapneumovirus reads')
    
    Metapneumovirus_reads = main(METAPMEUMOVIRUS_READS,
                                 read_lenght, 
                                 kmers,
                                 FILE_TYPE
                                 )
    
    print('Loading Rhinovirus reads')
    
    Rhinovirus_reads = main(RHINOVIRUS_READS,
                            read_lenght, 
                            kmers,
                            FILE_TYPE
                           )
    
    print('Loading SARS-CoV-2 reads')
    
    Sars_cov_2_reads = main(SARS_COV_2_READS,
                            read_lenght, 
                            kmers,
                            FILE_TYPE
                            )
    
    print('Loading Human reads')
    
    Human = main(HUMAN_READS,
                 read_lenght,
                 kmers,
                 FILE_TYPE
                 )
    
    total_sequences =  Coronaviridae_reads + \
                       Influenza_reads +\
                       Metapneumovirus_reads +\
                       Rhinovirus_reads +\
                       Sars_cov_2_reads +\
                       Human
    
    
    labels_to_fit = ['Coronaviridae','Influenza',"Metapneumovirus","Rhinovirus","Sars_cov_2", 'Human']
    label_maker = LabelBinarizer()
    transfomed_label = label_maker.fit(labels_to_fit)
    
    # save label_maker
    print('Saving object to convert output to labels '+ OUT_FOLDER+'/label_maker.'+MODEL_NAME+'.pickle')

    with open(OUT_FOLDER+'/label_maker.'+MODEL_NAME+'.pickle', 'wb') as handle:
        pickle.dump(label_maker, handle, protocol=pickle.HIGHEST_PROTOCOL)
    
    labels = list(np.repeat('Coronaviridae',len(Coronaviridae_reads))) + \
             list(np.repeat('Influenza',len(Influenza_reads))) + \
             list(np.repeat('Metapneumovirus',len(Metapneumovirus_reads))) + \
             list(np.repeat('Rhinovirus',len(Rhinovirus_reads))) + \
             list(np.repeat('Sars_cov_2',len(Sars_cov_2_reads))) + \
             list(np.repeat('Human',len(Human)))
             
    labels_proces = label_maker.transform(labels)
    
    # Tokenize the vocabulary
    tokenizer = Tokenizer()
    tokenizer.fit_on_texts(total_sequences)
    print('Converting reads into k-mers of lenght '+str(K_MERS))
    sequences_preproces = tokenizer.texts_to_sequences(total_sequences)
    
    max_length = max([len(s.split()) for s in total_sequences])
    # pad sequences
    sequences_preproces = pad_sequences(sequences_preproces, maxlen = max_length, padding = 'post')
    
    print('Saving tokenizer object '+ OUT_FOLDER+'/tokenizer.'+MODEL_NAME+'.pickle')
    with open(OUT_FOLDER+'/tokenizer.'+MODEL_NAME+'.pickle', 'wb') as handle:
        pickle.dump(tokenizer, handle, protocol=pickle.HIGHEST_PROTOCOL)
    
    sequences_preproces, labels_proces = shuffle(sequences_preproces, labels_proces)
    
    max_features = len(tokenizer.word_index)+1

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
    if GPU == True:
        model.add(Bidirectional(CuDNNLSTM(lstm_output_size)))
    else:
        model.add(Bidirectional(LSTM(lstm_output_size)))
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
        for chunks in range(0, len(sequences_preproces), 200000):
            start, end = chunks, chunks+200000
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
                                      validation_data=(X_test, y_test)
                                      )

            histories.append(chunk_history)
            if chunks == STOP_CHUNK:
                break
            break
    
    end = time.process_time()
    print('Traning time:', start  - end)
    
    # save keras model
    model.save(OUT_FOLDER+'/'+MODEL_NAME+".h5")
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
    plt.savefig(OUT_FOLDER+'/trainning_accuracy_'+MODEL_NAME+'.pdf',
                format='pdf',
                dpi=1200,
                bbox_inches='tight', pad_inches=0)
    
    f, ax = plt.subplots( figsize=(13,9))
    sns.lineplot(x=np.arange(len(histories)), y=np.array(loss), palette="tab10", linewidth=2.5, label='loss')
    plt.ylabel('Loss')
    plt.savefig(OUT_FOLDER+'/training_loss_'+MODEL_NAME+'.pdf',
                format='pdf',
                dpi=1200,
                bbox_inches='tight', pad_inches=0)
    
    f, ax = plt.subplots( figsize=(13,9))
    sns.lineplot(x=np.arange(len(histories)), y=np.array(val_bi_acc), palette="tab10", linewidth=2.5, label='Validation binary accuracy')
    sns.lineplot(x=np.arange(len(histories)), y=np.array(val_cat_acc), palette="tab10", linewidth=2.5, label='Validation categorical accuracy')
    plt.ylabel('Percentage of predicted reads')
    plt.savefig(OUT_FOLDER+'/val_training_accuracy_'+MODEL_NAME+'.pdf',
                format='pdf',
                dpi=1200,
                bbox_inches='tight', pad_inches=0)
    
    f, ax = plt.subplots( figsize=(13,9))
    sns.lineplot(x=np.arange(len(histories)), y=np.array(val_loss), palette="tab10", linewidth=2.5, label='Validation loss')
    plt.ylabel('Loss')
    plt.savefig(OUT_FOLDER+'/val_loss_'+MODEL_NAME+'.pdf',
                format='pdf',
                dpi=1200,
                bbox_inches='tight', pad_inches=0)
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    




