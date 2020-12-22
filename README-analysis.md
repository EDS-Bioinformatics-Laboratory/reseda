# Data analysis of a test 10X dataset

There are fastq files: R1/R2 and I1/I2 on L001/L002 lanes for the same sample (at least, that is what I guess from the name).
I1/I2 seem to be indices and R1/R2 contain T or B cell data. I ran PEAR on the R1/R2 files and there is not much overlap. Around 2.5%.

## Ran the CDR3 extraction script on R1, R2 and the assembled reads

In this example for the B-cell heavy chain (IGH_HUMAN).

```
mv reftables/ref.table.heavy.csv .
nohup python2 TranslateAndExtractCdr3.py -c IGH_HUMAN sc5p_v2_hs_B_1k_b_S1_L001_R* sc5p_v2_hs_B_1k_b_S1_L002_R* *.assembled.fastq.gz > nohup-cdr3.out 2> nohup-cdr3.err < /dev/null &
```

Result: not many CDR3's found.
I analyzed a few sequences with NCBI's IgBlast. 
The fragments are relatively short and contain different parts of the BCR.
From the data analysis results provided by 10X it turns out that sequences need to be assembled per cell/umi first.
I will install and use CellRanger, at least for the first assembly step.
I removed the results of the step above (TranslateAndExtractCdr3.py)
