# PACIFIC 

PACIFIC: A lightweight alignment-free deep-learning classifier of SARS-CoV-2 and co-infecting viral sequences  
#Add reference

PACIFIC implements deep learning to classify RNA sequencing reads into human, SARS-CoV-2 or additional respiratory viruses.

PACIFIC an input FASTA/FASTQ file and predicts the presence of the following viruses:
	    SARS-CoV-2
	    128 taxonomic units from Influenza
	    5 species from Metapneumovirus
	    130 species from Rhinovirus 
	    11 species from Coronaviridae (non-SARS-CoV-2).

## Table of Contents

1. Quick start
2. System requirements
3. Input and output
4. Usage
5. Test data

## Quick start

### Install and test PACIFIC
```
git clone https://github.com/pabloacera/pacific.git
cd pacific;
python ./PACIFIC.py --FILE_IN testdata etc.
#Add test data to check if it works
```
#Describe expected output


### System requirements
Python 3.X+ (python.org/) with the following libraries:
  Bio (v.)
  sklearn (v.)
  numPy (v.)
  keras (v.)
  tensorflow (v.)

## Input 
PACIFIC expects FASTA or FASTQ RNA-seq files as input. 

#Comment on training files

## Output
PACIFIC will output the following files:

#Describe output files

## Usage

Run PACIFIC
```
Usage: python ./PACIFIC.py [options]
```

Required arguments:
```
  --FILE_IN FILE_IN     FASTA/FASTQ file path to use PACIFIC
  --model MODEL         PACIFIC model path PACIFIC
  --tokenizer TOKENIZER
                        Tokenizer file path
  --label_maker LABEL_MAKER
                        Label maker object file path
  --file_type FILE_TYPE
                        fasta or fastq training files format (all files should
                        have same format)
```

Optional arguments:
```
  -h, --help            show this help message and exit
  --FILE_OUT FILE_OUT   path to the output file
  --k_mers K_MERS       K-mer number use to train the model
  --prediction_threshold PREDICTION_THRESHOLD
                        Threshold to use for the prediction
```

# Test data  


