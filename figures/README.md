## Data availability

Virus genome assembly data and human cDNA data used in this study are available at https://cloudstor.aarnet.edu.au/plus/s/sRLwF3IJQ12pNGQ.
1. `referencedata` directory contains merged assembly files per class
2. `benchmarkdata` directory contains synthetic data used for benchmarking PACIFIC (precision, recall, accuracy, etc)

## Benchmark tests

1. Benchmark data were generated using [generatebenchmarkdata.pl](https://github.com/pacific-2020/pacific/blob/master/figures/generatebenchmarkdata.pl) script.
2. Summary counts for benchmark tests is available in [benchmarkresults.txt](https://github.com/pacific-2020/pacific/blob/master/figures/benchmarkresults.txt).
3. Results of benchmarking were processed and analysed using [benchmark.R](https://github.com/pacific-2020/pacific/blob/master/figures/benchmark.R).

```
#testdata files for this experiment are in benchmarkdata directory
cd benchmarkdata
gunzip *.fasta.gz
printf %s\\n {1..100} | xargs -n 1 -P 4 -I{} \
python ~/pacific/PACIFIC.py \
-i testdata{}.fasta \
--model pacific.01.pacific_9mers_nonGPU.h5 \
--tokenizer tokenizer.01.pacific_9mers.pickle \
--label_maker label_maker.01.pacific_9mers.pickle \
--filte_type fasta
```

