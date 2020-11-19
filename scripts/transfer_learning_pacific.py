#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Apr  5 21:42:12 2020

This script will add new species to PACIFIC using transfer learning
Notice that at the moment PACIFIC only support 150bp

@author: Pablo Acera
"""

import argparse

parser = argparse.ArgumentParser(description=
                                 """
                                 This script train PACIFIC using fasta reads.
                                 The model, training and validation 
                                 plots will be generated.
                                 Also a model, tokenizer and label_maker will be generated
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

REQUIRED.add_argument("-t", "--tokenizer",
                      help="Tokenizer file path",
                      metavar='\b',
                      required=True)

REQUIRED.add_argument("--new_species_path",
                      help='coma separated path to the folder or folders containing new species of virus to train',
                      required=True,
                      )

REQUIRED.add_argument("--new_species_name",
                      help='coma separated name of the new species to train',
                      required=True,
                      )


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

OPTIONAL.add_argument("--GPU",
                      help='If True  PACIFIC will be train using CuDNNLSTM',
                      default=False,
                      )

OPTIONAL.add_argument("--file_type",
                      help='fasta or fastq training files format (all files should have same format)',
                      default='fasta',
                      )

OPTIONAL.add_argument("--accuracy_limit",
                      help='Stop the training when all individual classes accuracies reaches that level',
                      default=0.999,
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
NEW_SPECIES_PATH = ARGS.new_species_path
NEW_SPECIES_NAME = ARGS.new_species_name
TOKENIZER = ARGS.tokenizer

# Arguments
OUT_FOLDER = ARGS.out_folder
K_MERS = ARGS.k_mers
MODEL_NAME = ARGS.model_name
GPU = ARGS.GPU
FILE_TYPE = ARGS.file_type
ACCURACY_LIMIT = ARGS.accuracy_limit

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
import pickle
from datetime import datetime
from sklearn.utils import shuffle
from keras.models import load_model
import matplotlib.pyplot as plt
import seaborn as sns
import keras

tf.random.set_seed(42)

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
                               K_MERS,
                               FILE_TYPE,
                              )
    
    print('Loading Influenza reads')
    Influenza_reads = main(INFLUENZA_READS,
                           read_lenght,
                           K_MERS,
                           FILE_TYPE
                          )
    
    print('Loading Metapneumovirus reads')
    Metapneumovirus_reads = main(METAPMEUMOVIRUS_READS,
                                 read_lenght, 
                                 K_MERS,
                                 FILE_TYPE
                                 )
    
    print('Loading Rhinovirus reads')
    Rhinovirus_reads = main(RHINOVIRUS_READS,
                            read_lenght, 
                            K_MERS,
                            FILE_TYPE
                           )
    
    print('Loading SARS-CoV-2 reads')
    Sars_cov_2_reads = main(SARS_COV_2_READS,
                            read_lenght, 
                            K_MERS,
                            FILE_TYPE
                            )
    
    print('Loading Human reads')
    Human = main(HUMAN_READS,
                 read_lenght,
                 K_MERS,
                 FILE_TYPE
                 )
    
  
    print('loading new species')
    NEW_SPECIES_PATH = NEW_SPECIES_PATH.split(',')
    NEW_SPECIES_NAME  = NEW_SPECIES_NAME.split(',')
    new_species_sequences = {}
    for i in enumerate(NEW_SPECIES_PATH):
        new_species_sequences[NEW_SPECIES_NAME[i[0]]] = main(i[1],
                                                             read_lenght,
                                                             K_MERS,
                                                             FILE_TYPE
                                                             )
    
    total_sequences =  Coronaviridae_reads + \
                       Influenza_reads +\
                       Metapneumovirus_reads +\
                       Rhinovirus_reads +\
                       Sars_cov_2_reads +\
                       Human
    
    for i in enumerate(new_species_sequences.keys()):
        total_sequences += new_species_sequences[i[1]]
    
    labels_to_fit = ['Coronaviridae','Influenza',"Metapneumovirus","Rhinovirus","Sars_cov_2", 'Human']
    
    for i in enumerate(NEW_SPECIES_NAME):
        labels_to_fit += [i[1]]
    
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
             
    for i in enumerate(NEW_SPECIES_NAME):
        labels += list(np.repeat(i[1], len(new_species_sequences[i[1]])))
    
    labels_proces = label_maker.transform(labels) 
    
    # Import the tokenizer already trainned
    with open(TOKENIZER, 'rb') as handle:
        tokenizer = pickle.load(handle)
    
    tokenizer.fit_on_texts(total_sequences)
    print('Converting reads into k-mers of lenght '+str(K_MERS))
    sequences_preproces = tokenizer.texts_to_sequences(total_sequences)
    
    # pad sequences
    sequences_preproces = pad_sequences(sequences_preproces, maxlen = 142, padding = 'post')
    
    print('Saving tokenizer object '+ OUT_FOLDER+'/tokenizer.'+MODEL_NAME+'.pickle')
    with open(OUT_FOLDER+'/tokenizer.'+MODEL_NAME+'.pickle', 'wb') as handle:
        pickle.dump(tokenizer, handle, protocol=pickle.HIGHEST_PROTOCOL)
    
    sequences_preproces, labels_proces = shuffle(sequences_preproces, labels_proces)
    
    number_labels = 6 + len(NEW_SPECIES_NAME)
    
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
    
    model.load_weights('/media/labuser/Data/pacific/model/pacific.01.pacific_9mers_nonGPU.h5')
    
    #Delete last layers to change the number of neurons
    model.pop()
    model.pop()
    
    model.add(Dense(len(labels_to_fit)))
    
    model.add(Activation('softmax'))
    
    model.compile(loss='categorical_crossentropy',
                  optimizer='Adam',
                  metrics=['binary_accuracy', 
                           'categorical_accuracy',
                           ])
         
    model.summary()

    # training time
    #now = datetime.now().time() # time object

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
            
            inverser = label_maker.inverse_transform(y_test)
            accuracies_per_class = []
            for i in labels_to_fit:
                index = np.where(inverser==i)
                X_test_subselect = X_test[list(index[0])]   
                y_test_subselect = y_test[list(index[0])]
                print('accuracy '+i)
                predictions = np.where(model.predict(X_test_subselect) > 0.5, 1, 0)
                right = 0
                for j in enumerate(predictions):
                    if np.argmax(j[1]) == np.argmax(y_test_subselect[j[0]]):
                        right +=1
                print(right/len(predictions))
                print()
                accuracies_per_class.append(right/len(predictions))
            
            histories.append(chunk_history)

            # check if all classes have equal or more than 0.99 accuracy
            if all([i>=ACCURACY_LIMIT for i in accuracies_per_class]):
                #now = datetime.now().time() # time object
                print()
                break
    
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
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    




