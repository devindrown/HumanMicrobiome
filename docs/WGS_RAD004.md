# WGS via RAD004

Set your team name (e.g. 'TEAM=T.Alaska')
```
TEAM=T.YOURNAME
```


## Basecalling

Already completed
```
guppy_basecaller --input_path rawdata --save_path $TEAM.basecalled --flowcell FLO-MIN106 --kit SQK-RAD004  -q 0 –r --device cuda:0
```

Previous OUTPUT
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

Found 38 fast5 files to process.
Init time: 3426 ms

0%   10   20   30   40   50   60   70   80   90   100%
|----|----|----|----|----|----|----|----|----|----|
***************************************************
Caller time: 285000 ms, Samples called: 2185773886, samples/s: 7.66938e+06
Finishing up any open output files.
Basecalling completed successfully.
```
## Quality Control

Create a directory for your FASTQ/FASTA and combine basecalled data
```
mkdir $TEAM.fasta
cat basecalled/*.fastq > $TEAM.fasta/$TEAM.reads.fastq
```

Trim adapters with `porechop`
```
porechop -i $TEAM.fasta/$TEAM.reads.fastq -o $TEAM.fasta/$TEAM.reads.trim.fastq --format fastq -t 40 --discard_middle --check_reads 1000
```

Filter low quality reads (q score < 10)
```
filtlong --min_mean_q 90 $TEAM.fasta/$TEAM.reads.trim.fastq > $TEAM.fasta/$TEAM.reads.qcreads.fastq
```

# Variant Calling

Below script `map2variant` will use `nanopolish` to map reads to **READS*** to reference **GENOME**

```
map2variant.sh -i $TEAM.fasta/$TEAM.reads.qcreads.fastq -r rawdata -s basecalled/sequencing_summary.txt -g $HOME/genomes/AfricanSwineFever_Georgia2007_1_complete_genome.fasta -d 5 -f 0.75 -o $TEAM.variants -p Georgia2007_1 -t 48
```
