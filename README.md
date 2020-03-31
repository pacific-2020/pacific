# SRA downloads

Reads were downloaded using the following command:

```
fasterq-dump --threads 4 -p --split-files RunID
```

## NOTES
1. PRJNA605907: Reads were split into paired and unpaired reads. There were three files for each run, RunID_1.fastq.gz, RunID_2.fastq.gz and RunID.fastq.gz 
