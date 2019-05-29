# AMP via LSK109

Move to Amplicon Library data
```
cd AMP_demo
```

Set your team name (e.g. 'TEAM=T.Alaska')
```
TEAM=T.YOURNAME
```


## Basecalling

Already completed
```
guppy_basecaller --input_path rawdata --save_path $TEAM.basecalled --flowcell FLO-MIN106 --kit SQK-LSK109  -q 0 –r --device cuda:0
```

OUTPUT
```
ONT Guppy basecalling software version 3.1.5+781ed57
config file:        /home/up9demo/ont-guppy/data/dna_r9.4.1_450bps_hac.cfg
model file:         /home/up9demo/ont-guppy/data/template_r9.4.1_450bps_hac.jsn
input path:         rawdata
save path:          basecalled
chunk size:         1000
chunks per runner:  1000
records per file:   0
num basecallers:    4
gpu device:         cuda:0
kernel path:
runners per device: 2

Found 111 fast5 files to process.
Init time: 1348 ms

0%   10   20   30   40   50   60   70   80   90   100%
|----|----|----|----|----|----|----|----|----|----|
***************************************************
Caller time: 527764 ms, Samples called: 3924332994, samples/s: 7.43577e+06
Finishing up any open output files.
Basecalling completed successfully.
```
## Demultiplexing

Split data by barcode first with `guppy_barcoder`
```
guppy_barcoder -i basecalled/ -s $TEAM.demultiplexed_guppy -t 16 -q 0 --barcode_kits "EXP-NBD104" --device cuda:0
```
Split again with `porechop` and trim adapters
```
porechop -i demultiplexed_guppy -b $TEAM.demultiplexed --format fastq --discard_middle --threads 40 --check_reads 1000
```

## Filtering

Create a directory for your FASTQ/FASTA
```
mkdir $TEAM.fasta
````

Filter low quality reads (q score < 10)
```
filtlong --min_mean_q 90 demultiplexed/BC01.fastq > $TEAM.fasta/$TEAM.BC01.reads.qcreads.fastq
filtlong --min_mean_q 90 demultiplexed/BC02.fastq > $TEAM.fasta/$TEAM.BC02.reads.qcreads.fastq
filtlong --min_mean_q 90 demultiplexed/BC03.fastq > $TEAM.fasta/$TEAM.BC03.reads.qcreads.fastq
filtlong --min_mean_q 90 demultiplexed/BC04.fastq > $TEAM.fasta/$TEAM.BC04.reads.qcreads.fastq
filtlong --min_mean_q 90 demultiplexed/BC05.fastq > $TEAM.fasta/$TEAM.BC05.reads.qcreads.fastq
```


## Variant Calling

Below script `map2variant` will use `nanopolish` to map reads to **READS*** to reference **GENOME**
```
map2variant.sh -i $TEAM.fasta/$TEAM.BC01.reads.qcreads.fastq -r rawdata -s basecalled/sequencing_summary.txt -g $HOME/genomes/AfricanSwineFever_Georgia2007_1_complete_genome.fasta -d 100 -f 0.75 -o $TEAM.BC01.variants -p Georgia2007_1 -t 48
```
