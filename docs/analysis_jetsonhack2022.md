# Jetson Hackathon 2022 Analysis

## Software list
* Filtlong (https://github.com/rrwick/Filtlong)
* Seqtk (https://github.com/lh3/seqtk)
* Samtools (http://www.htslib.org/download/)
* Assembly-stats (https://github.com/sanger-pathogens/assembly-stats)


### Filtlong
Build from source
````
git clone https://github.com/rrwick/Filtlong.git
cd Filtlong; make -j
bin/filtlong -h
````
If this works, then copy to PATH
````
sudo cp bin/filtlong /usr/local/bin
````

### Seqtk
Build from source
````
git clone https://github.com/lh3/seqtk.git;
cd seqtk; make
````
If this works, then copy to PATH
````
sudo cp seqtk /usr/local/bin
````
## Samtools
````
sudo apt install samtools
````

### Assembly-stats
Install cmake first
````
sudo apt install cmake
````
Then build from source
````
git clone https://github.com/sanger-pathogens/assembly-stats.git
cd assembly-stats
mkdir build
cd build
cmake ..
make
make test
sudo make install
````

# Analysis

Change to analysis directory
````
cd /data/minknow/
sudo chown -R minit.minit jetsonhack2022
cd jetsonhack2022
ls -l
````
Pick run1 and enter that below for example, if you color is `RAINBOW`
````
SAMPLEID=jetsonhack_RAINBOW_run1_20220207
````

## Filter out length, and separate

### Setup paramters
````
MAXCHAN=256
THREADS=6
STORAGE=/data/minknow/jetsonhack2022
READS="$STORAGE"/"$SAMPLEID"/fasta/"$SAMPLEID".reads.pass.fastq
FILTERED="$STORAGE"/"$SAMPLEID"/fasta/"$SAMPLEID".reads.pass.l1000.fastq
````
### Collect all reads and filter by length
````
mkdir "$STORAGE"/"$SAMPLEID"/fasta
cat "$STORAGE"/"$SAMPLEID"/basecalling/pass/*.fastq > "$READS"
# filter out short (<1000 bases)
filtlong --min_length 1000 "$READS" > "$FILTERED"
````
#### Split by treatment

Enriched Batch
````
TRMT=enrich
OUTPUTFASTQ="$STORAGE"/"$SAMPLEID"/fasta/"$SAMPLEID".reads."$TRMT".fastq
OUTPUT_SUM="$STORAGE"/"$SAMPLEID"/fasta/"$SAMPLEID".reads."$TRMT".seqsum.txt
LIST_OUT="$STORAGE"/"$SAMPLEID"/fasta/"$SAMPLEID".reads."$TRMT".lst
````
````
head -n1 "$SAMPLEID"/basecalling/sequencing_summary.txt > "$OUTPUT_SUM"
awk -v c="$MAXCHAN" '$5 <= c && $10 == "TRUE" && $14 >= 1000 {print}' \
  "$STORAGE"/"$SAMPLEID"/basecalling/sequencing_summary.txt >> "$OUTPUT_SUM"
awk -v c="$MAXCHAN" '$5 <= c && $10 == "TRUE" && $14 >= 1000 {print $2}' \
  "$STORAGE"/"$SAMPLEID"/basecalling/sequencing_summary.txt > "$LIST_OUT"
seqtk subseq "$FILTERED" "$LIST_OUT" > "$OUTPUTFASTQ"
````
Control Batch
````
TRMT=control
OUTPUTFASTQ="$STORAGE"/"$SAMPLEID"/fasta/"$SAMPLEID".reads."$TRMT".fastq
OUTPUT_SUM="$STORAGE"/"$SAMPLEID"/fasta/"$SAMPLEID".reads."$TRMT".seqsum.txt
LIST_OUT="$STORAGE"/"$SAMPLEID"/fasta/"$SAMPLEID".reads."$TRMT".lst
````
````
head -n1 "$SAMPLEID"/basecalling/sequencing_summary.txt > "$OUTPUT_SUM"
awk -v c="$MAXCHAN" '$5 > c && $10 == "TRUE" && $14 >= 1000 {print}' \
  "$STORAGE"/"$SAMPLEID"/basecalling/sequencing_summary.txt >> "$OUTPUT_SUM"
awk -v c="$MAXCHAN" '$5 > c && $10 == "TRUE" && $14 >= 1000 {print $2}' \
  "$STORAGE"/"$SAMPLEID"/basecalling/sequencing_summary.txt > "$LIST_OUT"
seqtk subseq "$FILTERED" "$LIST_OUT" > "$OUTPUTFASTQ"
````

## Map output to reference genomes

### Select appropriate Reference Genome. 
Choose only **one** and remove the `#` in front
````
# REFERENCE=/data/genomes/Bacillus_subtilis_complete_genome.fasta
# REFERENCE=/data/genomes/Enterococcus_faecalis_complete_genome.fasta
# REFERENCE=/data/genomes/Escherichia_coli_complete_genome.fasta
# REFERENCE=/data/genomes/Listeria_monocytogenes_complete_genome.fasta
# REFERENCE=/data/genomes/Pseudomonas_aeruginosa_complete_genome.fasta
# REFERENCE=/data/genomes/Salmonella_enterica_complete_genome.fasta
````
### Create location for output
````
mkdir "$SAMPLEID"/mapped
````

### Map control and enriched reads to reference and output those reads
````
for TRMT in control enrich; do
  FASTQ="$SAMPLEID"/fasta/"$SAMPLEID".reads."$TRMT".fastq
  MAPPED_OUT="$SAMPLEID"/mapped/"$SAMPLEID"."$TRMT".bam
  LIST_OUT="$SAMPLEID"/mapped/"$SAMPLEID"."$TRMT".lst
  SUM_OUT="$SAMPLEID"/fasta/"$SAMPLEID".mapped."$TRMT".seqsum.txt
  FASTQ_OUT="$SAMPLEID"/fasta/"$SAMPLEID".mapped."$TRMT".fastq
  ## Map reads to Reference
  minimap2 -ax map-ont -t "$THREADS" "$REFERENCE" "$FASTQ" | \
    samtools sort -T reads.tmp | samtools view -F 0x04 -b > "$MAPPED_OUT"
  ## Index BAM and Extract mapped reads
  samtools index "$MAPPED_OUT"
  samtools view -F 0x04 "$MAPPED_OUT" | cut -f1 | sort | uniq > "$LIST_OUT"
  seqtk subseq "$SAMPLEID"/fasta/"$SAMPLEID".reads.pass.fastq "$LIST_OUT" > "$FASTQ_OUT"
  head -n1 "$SAMPLEID"/basecalling/sequencing_summary.txt > "$SUM_OUT"
  grep -f "$LIST_OUT" "$SAMPLEID"/basecalling/sequencing_summary.txt >> "$SUM_OUT"
done
````

## Grab the results
````
assembly-stats  "$SAMPLEID"/fasta/*.fastq | grep -B 1 sum
````

For the 4 files that end in
* `mapped.enrich.fastq`
* `mapped.control.fastq`
* `reads.enrich.fastq`
* `reads.control.fastq`

Copy the values from `sum = ` and `n =` into the shared spreadsheet for further analysis.

## Explore sequencing summary
You can use NanoPlot Online (http://nanoplot.bioinf.be/) to explore any of the files that end in `seqsum.txt`

