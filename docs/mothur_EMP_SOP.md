# Overview
Here we will analyze our own data. We'll start with a small set.

# Getting started

* Open a File Explorer window and navigate to `Documents/mothur/`
* Confirm data files and databases

You should verify that you have the databases we used in the previous lab. If you are in doubt visit those instructions:

[UAF mothur MiSeq SOP](https://devindrown.github.io/HumanMicrobiome/mothur_MiSeq_SOP)

**Check the raw data folder**

 `rawdata_neg` confirm that you have **12 original data files**

```
ExtractionNEG-B-2017001057_S96_L001_R1_001.fastq
ExtractionNEG-B-2017001057_S96_L001_R2_001.fastq
ExtractionNEG-A-2017000993_S32_L001_R1_001.fastq
ExtractionNEG-A-2017000993_S32_L001_R2_001.fastq
PCR-NEG-SetA2-2017000969_S8_L001_R1_001.fastq
PCR-NEG-SetA2-2017000969_S8_L001_R2_001.fastq
PCR-NEG-SetA1-2017000977_S16_L001_R1_001.fastq
PCR-NEG-SetA1-2017000977_S16_L001_R2_001.fastq
PCR-NEG-SetB2-2017001033_S72_L001_R1_001.fastq
PCR-NEG-SetB2-2017001033_S72_L001_R2_001.fastq
PCR-NEG-SetB1-2017001041_S80_L001_R1_001.fastq
PCR-NEG-SetB1-2017001041_S80_L001_R2_001.fastq
```

## Starting the program

Now you are ready to start the local version of mothur, you can double click on the `mothur` file

**Setup a directory to put all of your output.**
```
set.dir(output=./negativeQC)
```

# Prepare the custom database

Like last time, we'll be aligning our sequencings to a currated set of 16s sequences to the Silva database. We need to trim our samples to just the v4 region that we actually amplified. We used the 515f and 806rb primers from the Earth Microbiome Project. I've already aligned those to the E. coli genome and found their positions (start=13861, end=23444).

**Customize database to our region of interest**
```
pcr.seqs(fasta=./silva.bacteria/silva.bacteria.fasta, start=13861, end=23444, keepdots=F)
```

**Rename the reference database**

```
rename.file(input=INPUT, new=silva.EMP.fasta)
```

# Make contigs

First we need to take the paired-end read data from the MiSeq and join the read 1 and read 2 into a single sequences or contig. For example, you will see two files: ExtractionNEG-B-2017001057_S96_L001_R1_001.fastq and ExtractionNEG-B-2017001057_S96_L001_R2_001.fastq.. The first and all those with R1 correspond to read 1 while the second and all those with R2 correspond to the second or reverse read. These sequences are 250 bp and overlap in the V4 region of the 16S rRNA gene; this region is about 253 bp long.

Last time, you had a file preapred for you that contained all of the raw data files. This time around, you'll be making and editing that file from scratch. Mothur helps to get you started. You give the program a location (directory) and type of file (fastq).

```
make.file(inputdir=./rawdata_neg, type=fastq, prefix=neg.stability)
```

**Fix the stability file**

Unfortunately, the current version of `mothur` has a bug and the file isn't quite right. You'll need to open the file in `notepad` and change the first column. You shoud give each sample a short name without spaces or dashes. 

Remeber, the first column is the name othe sample. The second column is the name of the forward read for that sample and the third columns in the name of the reverse read for that sample.

**Combine our two sets of reads for each sample and then to combine the data from all of the samples**

```
make.contigs(inputdir=./rawdata_neg, file=./negativeQC/neg.stability.files, processors=8)
```

**Report how many sequences are in each group**

**What do these sequences look like?**
```
summary.seqs(# Reducing sequencing errorsfasta=current)
```
**How many total sequences do you have? What's the distribution of the length?**

# Reducing sequencing errors

**Next, we want to get rid of some of the bad reads**

Note how long most of the reads are. We'll use that to trim out the reads that assembled backwards.

```
screen.seqs(fasta=current, group=current, maxambig=0, maxlength=311)
```

**How many sequences do you have now? What command do you have to use?**


Many of our sequences are duplicates of each other. It's computationally wasteful to align the same thing a bazillion times.

**Select only the unique sequences**

```
unique.seqs()
```

**How many sequences do you have now?**

**Generate a table where the rows are the names of the unique sequences and the columns are the names of the groups.**
```
count.seqs(name=current,group=current)
```

**Let's see where you are in terms of reducing your dataset**
```
summary.seqs(count=current)
```

# Aligning to the database

Now we need to align our sequences to the reference alignment. You’ve already trimmed the silva database to just the right region of 16s. If you were using a different section, you’d need to do this from scratch. See the MiSeq SOP online for details.

**Align your dataset to the database**
```
align.seqs(fasta=current, reference=./negativeQC/silva.EMP.fasta)
```

**Check the results**

```
summary.seqs(count=current)
```
 You'll see that the bulk of the sequences start at position 1 and end at position 9583. Deviants from the mode positions are likely due to an insertion or deletion at the terminal ends of the alignments. 

To make sure that everything overlaps the same region we'll re-run `screen.seqs`  
```
screen.seqs(fasta=current,count=current,summary=current,start=1, end=9583, maxhomop=8)
```

Check to see if our sequences overlap the same alignment coordinates
```
summary.seqs(fasta=current, count=current)
```


**Filter the sequences** to remove the overhangs at both ends and any columns in the alignment that only contain gap characters.

```
filter.seqs(fasta=current, vertical=T, trump=.)
```

Perhaps we’ve created some redundancy across our sequences by trimming the ends

**Re-run `unique.seqs`**
```
unique.seqs(fasta=current,count=current)
```

**How many unique sequences are left?**

**Further de-noise our sequences** using the `pre.cluster` command allowing for up to 2 differences between sequences.

```
pre.cluster(fasta=current, count=current, diffs=2)
```

**How many unique sequences are left?**

# We have removed as much sequencing error as we can!

**Check for chimeras**

```
chimera.vsearch(fasta=current, count=current, dereplicate=t)
```

**Remove those sequences from the fasta file.**

```
remove.seqs(fasta=current, accnos=current)
```

**Check your results**
```
summary.seqs(fasta=current, count=current)
```

**How many unique and total sequences are left now?**

Sometimes when we pick a primer set they will amplify other stuff such as 16S rRNA from chloroplasts, and mitochondria.

**Classify your sequences**

```
classify.seqs(fasta=current, count=current, reference=./reference/trainset9_032012.pds.fasta, taxonomy=./reference/trainset9_032012.pds.tax, cutoff=80)
```

Now that everything is classified we want to **remove our undesirables**

```
remove.lineage(fasta=current, count=current, taxonomy=current, taxon=Chloroplast-Mitochondria-unknown-Archaea-Eukaryota)
```

**Create an updated taxonomy summary file that reflects these removals**
```
summary.tax(taxonomy=current, count=current)
```
This creates a new summary file with the undesirables removed. At this point we have curated our data as far as possible

# Get some OTUs

**Clustering sequences into OTUs.**

We use the taxonomic information to bin our sequences. In this command we use taxlevel=4, which corresponds to the level of Order.

**Split the sequences into bins and then cluster within each bin**
```
cluster.split(fasta=current, count=current, taxonomy=current, splitmethod=classify, taxlevel=4, cutoff=0.03)
```

**Combine your data across different samplee**
```
make.shared(list=current, count=current, label=0.03)
```

**Get the consensus taxonomy for each OTU**
```
classify.otu(list=current, count=current, taxonomy=current, label=0.03)
```

# Quick Visualize

**Create a `.biom` file**
```
make.biom(shared=current,constaxonomy=current)
```

**Rename the file to something simpler**
```
rename.file(biom=current, new=TEAMNAME.biom)
```
