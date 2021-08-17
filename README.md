![picture](msc/pacific_brand.png)

## PACIFIC: A lightweight deep-learning classifier of SARS-CoV-2 and co-infecting viral sequences  
Read more about PACIFIC here:
> [__Acera Mateos P., Balboa R.F., Easteal S., Eyras E., and Patel, H.R.__ PACIFIC: A lightweight deep-learning classifier of SARS-CoV-2 and co-infecting RNA viruses. *Scientific Reports*, 2021.](https://www.nature.com/articles/s41598-021-82043-4)

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

Note: As the model file is ~300MB, we hosted this file on Cloudstor. The model file can be found [here](https://cloudstor.aarnet.edu.au/plus/s/Hwg20YRlua9a2OH). We recommend that users download this file and place it in the model directory as in the code below:

```
##Create required virtual environment for running PACIFIC
conda create --name pacific python=3.7.6

##Activate the environment
conda activate pacific

##Install required dependencies
conda install numpy=1.18.1 tensorflow=2.2.0 keras=2.3.1 pandas=1.0.1 scikit-learn=0.21.3 biopython=1.76

#Clone PACIFIC repository
git clone https://github.com/pacific-2020/pacific.git

#Download model file
cd pacific
wget -O model/pacific.01.pacific_9mers_nonGPU.h5 https://cloudstor.aarnet.edu.au/plus/s/Hwg20YRlua9a2OH/download

#Change to test directory and run PACIFIC
cd test
python ../PACIFIC.py \
  -i test.fa.gz \
  -m ../model/pacific.01.pacific_9mers_nonGPU.h5 \
  -t ../model/tokenizer.01.pacific_9mers.pickle \
  -l ../model/label_maker.01.pacific_9mers.pickle \ 
  -f fasta
```

If installed correctly, PACIFIC should generate pacificoutput_test.fa.gz and test.fa.gz_summary.txt in the test directory, and should provide the following results in the terminal:

```
From a total of 5000 reads, 0 were discarded (e.g. non-ACGT nucleotides/characters or short reads (<150bp))

   filename            class  # predicted reads  predicted reads (%)  # predicted reads above 0.95  predicted reads above 0.95 (%)
 test.fa.gz       SARS-CoV-2                  0                 0.00                             0                            0.00
 test.fa.gz    Coronaviridae               4997                99.94                          4997                           99.94
 test.fa.gz        Influenza                  0                 0.00                             0                            0.00
 test.fa.gz  Metapneumovirus                  0                 0.00                             0                            0.00
 test.fa.gz       Rhinovirus                  0                 0.00                             0                            0.00
 test.fa.gz            Human                  1                 0.02                             0                            0.00
 test.fa.gz        Discarded                  0                 0.00                             0                            0.00
 test.fa.gz     rc_discarded                  2                 0.04                             1                            0.02

Thank you for using PACIFIC =^)
```

## System requirements
- Python 3.X+ (python.org/) with the following libraries:
    - Bio 1.76
    - numPy 1.18.1
    - keras 2.3.1
    - pandas 1.0.1
    - tensorflow 2.2.0
    - scikit-learn 0.21.3
    - cudatoolkit 10.1.168 (only for GPU mode)
    - cudnn 7.6.0 (only for GPU mode)
    
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
                        FASTA or FASTQ training file format [fastq]
  -o <dir>, --outputdir <dir>
                        Path to output directory [.]
  -d <dir>, --tmpdir <dir>
                        Path to tmp directory [outputdir]
  -T <float>, --prediction_threshold <float>
                        Threshold/cutoff for predictions [0.95]
  -c <int>, --chunk_size <int>
                        Number of reads per chunk [100000]
  -v, --version         show program's version number and exit
```

## Input 
PACIFIC expects four arguments as input: 
 - FASTA or FASTQ RNA-seq file (can handle either gzipped or non-gzipped files)
 - Training model file (recommended: model/pacific.01.pacific_9mers_nonGPU.h5)
 - Tokenizer file (recommended: model/tokenizer.01.pacific_9mers.pickle)
 - Label maker file (recommended: model/label_maker.01.pacific_9mers.pickle)

PACIFIC allows users to use their own custom training model, tokenizer and label maker files. However, we recommend the use of default parameters and the following files above as input into the program.

## Output
PACIFIC will output the following files (using default parameters):

*1. pacificoutput_$input.gz:*
A gzipped fasta file with modified sequence headers from the input fasta file. PACIFIC includes the prediction class and score in the header, as well as whether the sequences are discarded when run through the program e.g. non-ACGT nucleotides/characters or short reads (<150bp) - "discarded" or whether a read fails to meet the reverse complement prediction test "rc_discarded"). 

For example, the following describes a sequence predicted to be of the Coronaviridae class with a prediction score of 0.97:

```
>fastaheader:0.97:Coronaviridae
```

*2. $input_summary.txt:*
A text file which summarises the number of predicted reads for each class and their predicted proportions (%) both in the entire dataset, as well as in predicted reads with a score above 0.95. The table also provides information on discarded and rc_discarded reads. This is the same information that is provided as output in the terminal when PACIFIC is run.

## Test and model data

1. Model and test data are available [here](https://cloudstor.aarnet.edu.au/plus/s/sRLwF3IJQ12pNGQ)
2. PACIFIC model can be downloaded from [here](https://cloudstor.aarnet.edu.au/plus/s/Hwg20YRlua9a2OH)
