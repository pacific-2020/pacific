# SRA downloads

Reads were downloaded using the following command:

```
fasterq-dump --threads 4 -p --split-files --outfile PLATFORM.ASSAYTYPE.LIBRARYLAYOUT/BIOPROJECTID_RUNID.fastq RunID
```

## NOTES
1. PRJNA605907: Reads were split into paired and unpaired reads. There were three files for each run, RunID_1.fastq.gz, RunID_2.fastq.gz and RunID.fastq.gz 
2. PRJNA614995: SRR11410536, SRR11410538, SRR11410540 are labelled as SINGLE library layout in metadata file. They are actually paired data. These data are treated as paired data.
3. PRJNA612578: SRR11314339 was split into paired and unpaired reads. There were three files for each run, RunID_1.fastq.gz, RunID_2.fastq.gz and RunID.fastq.gz

