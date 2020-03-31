# SRA downloads

Reads were downloaded using the following command:

```
fasterq-dump --threads 4 -p --split-files RunID
```

## NOTES
1. PRJNA605907: Reads were split into paired and unpaired reads. There were three files for each run, RunID_1.fastq.gz, RunID_2.fastq.gz and RunID.fastq.gz 
2. PRJNA614995: SRR11410540 is labelled as SINGLE library layout in metadata file, it is actually paired data. This data is treated as paired data.


