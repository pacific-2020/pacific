# PACIFIC 

PACIFIC: A lightweight alignment-free deep-learning classifier of SARS-CoV-2 and co-infecting viral sequences  
#Add reference

PACIFIC implements deep learning to classify RNA sequencing reads into human, SARS-CoV-2 or additional respiratory viruses.

## Table of Contents

1. Quick start
2. Usage

## Quick start

### Install and test PACIFIC
```
git clone https://github.com/pabloacera/pacific.git
cd pacific;
python ./PACIFIC.py --FILE_IN testdata etc.
#Add test data to check if it works
```

### System requirements
Python 3.X+ (python.org/) with the following libraries:
  Bio (v.)
  sklearn (v.)
  numPy (v.)
  keras (v.)
  tensorflow (v.)

## Usage

Run PACIFIC
```
python ./PACIFIC.py \
  --FILE_IN  my.--model MODEL --tokenizer TOKENIZER
                  --label_maker LABEL_MAKER [--file_type FILE_TYPE]
                  [--FILE_OUT FILE_OUT] [--k_mers K_MERS]
                  [--prediction_threshold PREDICTION_THRESHOLD]
###Set defaults here
```

Output
PACIFIC will output the following files:

  





# SRA downloads

Reads were downloaded using the following command:

```
fasterq-dump --threads 4 -p --split-files --outfile PLATFORM.ASSAYTYPE.LIBRARYLAYOUT/BIOPROJECTID_RUNID.fastq RunID
```

## NOTES
1. PRJNA605907: Reads were split into paired and unpaired reads. There were three files for each run, RunID_1.fastq.gz, RunID_2.fastq.gz and RunID.fastq.gz 
2. PRJNA614995: SRR11410536, SRR11410538, SRR11410540 are labelled as SINGLE library layout in metadata file. They are actually paired data. These data are treated as paired data.
3. PRJNA612578: SRR11314339 was split into paired and unpaired reads. There were three files for each run, RunID_1.fastq.gz, RunID_2.fastq.gz and RunID.fastq.gz

