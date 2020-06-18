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
1. [Test and model data](#Test-and-model-data)

## Quick start

### Install and test PACIFIC

Note: As the model file is ~300MB, we hosted this file Cloudstor.
We recommend that users download this file and place it in the model directory as in the code below:

```
git clone https://github.com/pacific-2020/pacific.git
cd pacific/model
wget MODEL FILE .
cd ../test
python ../PACIFIC.py \
  -i ./test.fa \
  -m ../model/pacific.01.pacific_9mers_nonGPU.h5 \
  -t ../model/tokenizer.01.pacific_9mers.pickle \
  -l ../model/label_maker.01.pacific_9mers.pickle
```

If installed correctly, PACIFIC should create output_PACIFIC.fasta and output_PACIFIC.txt in the test directory, and should provide the following results in the terminal:

```
From a total of 5000 reads, 0 were discarded (e.g. non-ACGT nucleotides/characters or short reads (<150bp))

             Class  # predicted reads  # predicted reads (%)  # predicted reads above 0.95  # predicted reads above 0.95 (%)
0       SARS-CoV-2                  0                   0.00                             0                               0.0
1    Coronaviridae               4998                  99.96                          4998                             100.0
2        Influenza                  0                   0.00                             0                               0.0
3  Metapneumovirus                  0                   0.00                             0                               0.0
4       Rhinovirus                  1                   0.02                             0                               0.0
5            Human                  1                   0.02                             0                               0.0

Thank you for using PACIFIC =^)
```

## System requirements
- Python 3.X+ (python.org/) with the following libraries:
    - Bio 1.74
    - sklearn 0.20.3
    - numPy 1.16.4
    - keras 2.2.4
    - pandas 0.25.1
    - tensorflow 2.2.0
    - scikit-learn 0.21.3
    - cudatoolkit 10.1.168
    - cudnn 7.6.0
    
  (for a full list of package versions, view metadata/pacific_versions.txt)

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
PACIFIC will output the following files (using default parameters):

- output_PACIFIC.fasta
A fasta file with modified sequence headers from the input fasta file. PACIFIC includes the prediction class and score in the header. For example, the following describes a sequence predicted to be of the Coronaviridae class with a prediction score of 0.97:

```
>fastaheader:0.97:Coronaviridae
```

- output_PACIFIC.txt
A csv file which summarises the number of predicted reads for each class and their predicted proportions (%)both in the entire dataset, as well as in predicted reads with a score above 0.95. This is the same information that is provided as output into the terminal when PACIFIC is run.

## Test and model data

1. Model and test data are available [here](https://cloudstor.aarnet.edu.au/plus/s/sRLwF3IJQ12pNGQ)
2. PACIFIC model is available [here](https://cloudstor.aarnet.edu.au/plus/s/Hwg20YRlua9a2OH)