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
make.file(inputdir=./rawdata_neg/, type=fastq, prefix=neg.stability)
```

**Fix the stability file**

Unfortunately, the current version of `mothur` has a bug and the file isn't quite right. You'll need to open the file in `notepad` and change the first column. You shoud give each sample a short name without spaces or dashes. 

Remeber, the first column is the name of the sample. The second column is the name of the forward read for that sample and the third columns in the name of the reverse read for that sample.

# Reducing sequencing and PCR errors

**Combine our two sets of reads for each sample and then to combine the data from all of the samples**

```
make.contigs(inputdir=./rawdata_neg, file=./negativeQC/neg.stability.files, processors=8)
```

**Report how many sequences are in each group**

**What do these sequences look like?**
```
summary.seqs(fasta=current)
```
**How many totoal sequences do you have? What's the distribution of the length?**

**Next, we want to get rid of some of the bad reads**

Note how long most of the reads are. We'll use that to trim out the reads that assembled backwards.

```
screen.seqs(fasta=current, group=current, maxambig=0, maxlength=311)
```

**It is a lot of work to keep typing in `fasta=example.fasta` etc. Try**
```
get.current()
```
mothur remembers your latest fasta file and group file as well as the number of processors you have, so you don’t have to type them in each time. You should see something like this list:

```
Current files saved by mothur:
fasta=/home/microbiome/mothur/testrun/stability.trim.contigs.good.fasta
group=/home/microbiome/mothur/testrun/stability.contigs.good.groups
qfile=/home/microbiome/mothur/testrun/stability.trim.contigs.qual
contigsreport=/home/microbiome/mothur/testrun/stability.contigs.report
processors=8
summary=/home/microbiome/mothur/testrun/stability.trim.contigs.summary

Current input directory saved by mothur: /home/microbiome/mothur/MiSeq_SOP/
Current output directory saved by mothur: /home/microbiome/mothur/testrun/
Current default directory saved by mothur: /home/microbiome/mothur/mothur/
Current working directory: /home/microbiome/mothur/

Output File Names:
/home/microbiome/mothur/testrun/current_files.summary
```

# Processing improved sequences

Many of our sequences are duplicates of each other. It's computationally wasteful to align the same thing a bazillion times.

**Select only the unique sequences**

```
unique.seqs()
```

**Compare your output with the expected below**
```
Output File Names:
stability.trim.contigs.good.names
stability.trim.contigs.good.unique.fasta
```

**Generate a table where the rows are the names of the unique sequences and the columns are the names of the groups.**
```
count.seqs(name=current,group=current)
```

**Compare your output with the expected below**
```
Total number of sequences: 128872

Output File Names: stability.trim.contigs.good.count_table
```

**Let's see where you are in terms of reducing your dataset**
```
summary.seqs(count=current)
```

**Compare your output with the expected below**
```
                Start   End     NBases  Ambigs  Polymer NumSeqs
Minimum:        1       250     250     0       3       1
2.5%-tile:      1       252     252     0       3       411
25%-tile:       1       252     252     0       4       4107
Median:         1       253     253     0       4       8214
75%-tile:       1       253     253     0       5       12320
97.5%-tile:     1       254     254     0       6       16016
Maximum:        1       270     270     0       12      16426
Mean:   1       252.594 252.594 0       4.44277
# of Seqs:      16426
total # of seqs:        128872

Output File Names:
/home/microbiome/mothur/testrun/stability.trim.contigs.good.unique.summary
```

# Aligning to the database

Now we need to align our sequences to the reference alignment. You’ve already trimmed the silva database to just the right region of 16s. If you were using a different section, you’d need to do this from scratch. See the MiSeq SOP online for details.

**Align your dataset to the database**
```
align.seqs(fasta=current, reference=./testrun/silva.v4.fasta)
```

**Check the results**
```
summary.seqs(count=current)
```
So what does this mean? You'll see that the bulk of the sequences start at position 1968 and end at position 11550. Deviants from the mode positions are likely due to an insertion or deletion at the terminal ends of the alignments. 

To make sure that everything overlaps the same region we'll re-run `screen.seqs` to get sequences that start at or before position 1 and end at or after position 11550. We'll also set the maximum homopolymer length to 8 since there's nothing in the database with a stretch of 9 or more of the same base in a row. 
```
screen.seqs(fasta=current,count=current,summary=current,start=1968, end=11550, maxhomop=8)
```

Check to see if our sequences overlap the same alignment coordinates
```
summary.seqs(fasta=current, count=current)
```

**Compare your output with the expected below**
```
 		Start	End	NBases	Ambigs	Polymer	NumSeqs
Minimum:	1965	11550	250	0	3	1
2.5%-tile:	1968	11550	252	0	3	3217
25%-tile:	1968	11550	252	0	4	32164
Median: 	1968	11550	252	0	4	64328
75%-tile:	1968	11550	253	0	5	96492
97.5%-tile:	1968	11550	253	0	6	125439
Maximum:	1968	13400	270	0	8	128655
Mean:	1968	11550	252.463	0	4.36666
# of unique seqs:	16298
total # of seqs:        128655
```


**Filter the sequences** to remove the overhangs at both ends. In addition, there are many columns in the alignment that only contain gap characters.
```
filter.seqs(fasta=current, vertical=T, trump=.)
```

**Compare your output with the expected below**
```
Length of filtered alignment: 376
Number of columns removed: 13049
Length of the original alignment: 13425
Number of sequences used to construct filter: 16298
```


Perhaps we’ve created some redundancy across our sequences by trimming the ends

**Re-run `unique.seqs`**
```
unique.seqs(fasta=current,count=current)
```

**Further de-noise our sequences** using the `pre.cluster` command allowing for up to 2 differences between sequences. This command will split the sequences by group and then sort them by abundance and go from most abundant to least and identify sequences that are within 2 nt of each other. If they are then they get merged. We generally favor allowing 1 difference for every 100 bp of sequence
```
pre.cluster(fasta=current, count=current, diffs=2)
```

# We have removed as much sequencing error as we can!

**Check for chimeras**

```
chimera.vsearch(fasta=current, count=current, dereplicate=t)
```

You still need to remove those sequences from the fasta file. Depending on the outcome of the chimera check, the following command may give an error. If there are no bad sequences to remove, then it will produce an error message.
```
remove.seqs(fasta=current, accnos=current)
```

**Check your results**
```
summary.seqs(fasta=current, count=current)
```

**Compare your output with the expected below**
```
                Start   End     NBases  Ambigs  Polymer NumSeqs
Minimum:        1       376     249     0       3       1
2.5%-tile:      1       376     252     0       3       2954
25%-tile:       1       376     252     0       4       29539
Median:         1       376     252     0       4       59077
75%-tile:       1       376     253     0       5       88615
97.5%-tile:     1       376     253     0       6       115200
Maximum:        1       376     256     0       8       118153
Mean:   1       376     252.464 0       4.37545
# of unique seqs:       2283
total # of seqs:        118153
```

Sometimes when we pick a primer set they will amplify other stuff such as 16S rRNA from chloroplasts, and mitochondria.

**Classify your sequences**
```
classify.seqs(fasta=current, count=current, reference=./reference/trainset9_032012.pds.fasta, taxonomy=./reference/trainset9_032012.pds.tax, cutoff=80)
```

Now that everything is classified we want to **remove our undesirables**
```
remove.lineage(fasta=current, count=current, taxonomy=current, taxon=Chloroplast-Mitochondria-unknown-Archaea-Eukaryota)
```

Also of note is that "unknown" only pops up as a classification if the classifier cannot classify your sequence to one of the domains. If you run summary.seqs you'll see that we now have 2281 unique sequences and a total of 118150 total sequences. This means about 350 of our sequences were in these various groups.

**Create an updated taxonomy summary file that reflects these removals**
```
summary.tax(taxonomy=current, count=current)
```
This creates a new summary file with the undesirables removed. At this point we have curated our data as far as possible

# Let's get some OTUs

**Remove a Mock community group**
```
remove.groups(count=current, fasta=current, taxonomy=current, groups=Mock)
```

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
rename.file(biom=current, new=testrun.biom)
```
