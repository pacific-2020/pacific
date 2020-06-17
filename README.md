# PACIFIC 

PACIFIC: A lightweight alignment-free deep-learning classifier of SARS-CoV-2 and co-infecting viral sequences  
(Add preprint link once submitted)

PACIFIC implements deep learning to classify RNA sequencing reads into human, SARS-CoV-2 or additional respiratory viruses. PACIFIC takes an input FASTA/FASTQ file and predicts the presence of the following viruses and their relative proportions within a sample:
- SARS-CoV-2
- 128 taxonomic units from Influenza
- 5 species from Metapneumovirus
- 130 species from Rhinovirus 
- 11 species from Coronaviridae (non-SARS-CoV-2)

## Table of Contents

1. [Quick start](#Quick-start)
1. [System requirements](#System-requirements)
1. [Usage](#Usage)
1. [Input](#Input)
1. [Output](#Output)

## Quick start

### Install and test PACIFIC
```
git clone https://github.com/pabloacera/pacific.git
cd pacific;
python ./scripts/PACIFIC.py -i ./testdata/testdata.fa -m ./model/pacific.01.pacific_9mers_nonGPU.h5  -t ./model/tokenizer.01.pacific_9mers.pickle -l ./model/label_maker.01.pacific_9mers.pickle
!!!!!!Add test data to check if it works for a user - maybe in a "testdata" directory as above
```
!!!!!!!Describe expected output

## System requirements
- Python 3.X+ (python.org/) with the following libraries:
    - Bio 1.74
    - sklearn 0.20.3
    - numPy 1.16.4
    - keras 2.2.4
    - pandas 0.25.1
    - tensorflow 1.14.0
    - scikit-learn 0.21.3
    - cudatoolkit 10.1.168
    - cudnn 7.6.0
    
  (for more packages versions, look at pacific_versions.txt file)

## Usage

**Run PACIFIC**
```
usage: python PACIFIC.py [options] -i <in.fa>|<in.fq> -m <model> -t <tokenizer> -l <label-maker>
```

**Required arguments:**
```
  -i, --input_file  FASTA/FASTQ input file path
  -m, --model       PACIFIC model file path
  -t, --tokenizer   Tokenizer file path
  -l, --label_maker Label maker object file path
```

**Optional arguments:**
```
  -h, --help            show this help message and exit
  -f <fasta/fastq>, --file_type <fasta/fastq>
                        FASTA or FASTQ training file format [fasta]
  -o <dir>, --outputdir <dir>
                        Path to output directory [.]
  -T <float>, --prediction_threshold <float>
                        Threshold/cutoff for predictions [0.95]
  -c <int>, --chunk_size <int>
                        Number of reads per chunk [10000]
  -O, --output_fasta    If this option is "True", a FASTA file containing
                        predictions for each read will be provided [False]
  -v, --version         show program's version number and exit
```

## Input 
PACIFIC expects four arguments as input: 
 - FASTA or FASTQ RNA-seq file # Multiple files accepted?
 - Training model file (recommended: ./model/pacific.01.pacific_9mers_nonGPU.h5)
 - Tokenizer file (recommended: ./model/tokenizer.01.pacific_9mers.pickle)
 - Label maker file (recommended: ./model/label_maker.01.pacific_9mers.pickle)

PACIFIC allows users to use their own custom training model, tokenizer and label maker files. However, we recommend the use of default parameters and the following files above as input into the program.

## Output
PACIFIC will output the following files:

!!!!!Describe output files

