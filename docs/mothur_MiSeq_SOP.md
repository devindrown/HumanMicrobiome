# Overview

What follows is a modified version of the mothur MiSeq SOP. You can find the full details at:

[mothur MiSeq SOP](https://www.mothur.org/wiki/MiSeq_SOP)

The author of that manual encourages you to cite the following paper if you use this SOP:

> Kozich JJ, Westcott SL, Baxter NT, Highlander SK, Schloss PD. (2013): Development of a dual-index sequencing strategy and curation pipeline for analyzing amplicon sequence data on the MiSeq Illumina sequencing platform. Applied and Environmental Microbiology. 79(17):5112-20.

# Getting started with mothur

* In the Files window and make sure that you are in `Home`

## Confirm data files and databases

Already included on your hard drive are several data files, along with some databases. Let’s have a quick look at the files on the system at the start before we dig any deeper. 

**Check the following folders**

 `MiSeq_SOP` confirm that you have original data files
```
F3D0_S188_L001_R1_001.fastq
F3D0_S188_L001_R2_001.fastq
...
stability.files
```
`silva.bacteria` confirm seven files with different extensions such as:
```
silva.bacteria.fasta
silva.bacteria.rdp.tax
...
```
`reference`, confirm you have database files:
```
trainset9_032012.pds.fasta
trainset9_032012.pds.tax
```

## Starting the program

Now you are ready to start the local version of mothur, Go to the Terminal window and it should look like
```
mb2021@hikita:~$
```
At the `$` prompt type `mothur`

**Verify with the instructor that you are indeed running running: mothur v.1.42.3.**

From here on out you’ll be entering commands directly to the mothur program. The command prompt should look like this:
```
mothur>
```

You'll find that `mothur` generates a lot of files. It is helpful to a single location to help organize.

**Setup a directory to put all of your output.**
```
set.dir(output=./testrun)
```

# Prepare the custom database

We'll be aligning our sequencings to a currated set of 16S sequences called the Silva database. We'll go more into this later, but for now we want to trim the database down to just the v4 region.

**Customize database to our region of interest**
```
pcr.seqs(fasta=./silva.bacteria/silva.bacteria.fasta, start=11894, end=25319, keepdots=F, processors=4)
```
The output should look something like
```
It took 2 secs to screen 14956 sequences.

Output File Names:
/home/mb2021/testrun/silva.bacteria.pcr.fasta
```

**Rename the output**
Paying cafeul attention to the output file, use this information in the below command substituting **INPUT** with the file location
```
rename.file(input=INPUT, new=silva.v4.fasta)
```

# Reducing sequencing and PCR errors

**Combine our two sets of reads for each sample and then to combine the data from all of the samples**
```
make.contigs(inputdir=./MiSeq_SOP, file=./MiSeq_SOP/stability.files, processors=4)
```

You’ll might get a [WARNING] message, but don’t worry.

**Compare your output with the expected below**
Output
```
Group count:
F3D0    7793
F3D1    5869
F3D141  5958
...
F3D8    5294
F3D9    7070
Mock    4779
Total of all groups is 152360

Output File Names:
/home/mb2021/testrun/stability.trim.contigs.fasta
/home/mb2021/testrun/stability.trim.contigs.qual
/home/mb2021/testrun/stability.scrap.contigs.fasta
/home/mb2021/testrun/stability.scrap.contigs.qual
/home/mb2021/testrun/stability.contigs.report
/home/mb2021/testrun/stability.contigs.groups

```

**Go back to the Files window, click on the `MiSeq_SOP` folder, and then open `stability.files`**

What do you see?  The first column is the name of the sample. The second column is the name of the forward read for that sample and the third columns in the name of the reverse read for that sample.

**Let's see what these sequences look like**, try getting a summary by typing
```
summary.seqs(fasta=current)
```

**Compare your output with the expected below**

```
Using /home/mb2021/testrun/stability.trim.contigs.fasta as input file for the fasta parameter.
Using 4 processors.
                Start   End     NBases  Ambigs  Polymer NumSeqs
Minimum:        1       248     248     0       3       1
2.5%-tile:      1       252     252     0       3       3810
25%-tile:       1       252     252     0       4       38091
Median:         1       252     252     0       4       76181
75%-tile:       1       253     253     0       5       114271
97.5%-tile:     1       253     253     6       6       148552
Maximum:        1       502     502     249     243     152360
Mean:   1       252     252     0       4
# of Seqs:      152360

It took 0 secs to summarize 152360 sequences.
Output File Names:
/home/mb2021/testrun/stability.trim.contigs.summary
```

**Verify that you have the same number of reads `# of Seqs:      152360`**

**Next, we want to get rid of some of the bad reads**

```
screen.seqs(fasta=current, group=current, maxambig=0, maxlength=275)
```

**It is a lot of work to keep typing in `fasta=example.fasta` etc. Try**
```
get.current()
```
`mothur` remembers your latest fasta file and group file as well as the number of processors you have, so you don’t have to type them in each time. You should see something like this list:

```
Current RAM usage: 0.0926437 Gigabytes. Total Ram: 125.848 Gigabytes.
Current files saved by mothur:
accnos=/home/mb2021/testrun/stability.trim.contigs.bad.accnos
fasta=/home/mb2021/testrun/stability.trim.contigs.good.fasta
group=/home/mb2021/testrun/stability.contigs.good.groups
qfile=/home/mb2021/testrun/stability.trim.contigs.qual
contigsreport=/home/mb2021/testrun/stability.contigs.report
processors=4
summary=/home/mb2021/testrun/stability.trim.contigs.summary

Current input directory saved by mothur: /home/mb2021/MiSeq_SOP/
Current output directory saved by mothur: /home/mb2021/testrun/
Current default directory saved by mothur: /usr/local/bin/
Current working directory: /home/mb2021/

Output File Names:
/home/mb2021/testrun/current_files.summary
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
count.seqs(name=current, group=current)
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
Mean:   1       252     252     0       4
# of Seqs:      16426
total # of seqs:        128872

Output File Names:
/home/mb2021/testrun/stability.trim.contigs.good.unique.summary
```

# Aligning to the database

Now we need to align our sequences to the reference alignment. You’ve already trimmed the silva database to just the right region of 16S. If you were using a different section, you’d need to do this from scratch. See the MiSeq SOP online for details.

**Align your dataset to the database**
```
align.seqs(fasta=current, reference=./testrun/silva.v4.fasta)
```

**Check the results**
```
summary.seqs(count=current)
```
So what does this mean? You'll see that the bulk of the sequences start at position 1968 and end at position 11550. Deviants from the mode positions are likely due to an insertion or deletion at the terminal ends of the alignments. 

To make sure that everything overlaps the same region we'll re-run `screen.seqs` to get sequences that start at or before position 1968 and end at or after position 11550. We'll also set the maximum homopolymer length to 8 since there's nothing in the database with a stretch of 9 or more of the same base in a row. 
```
screen.seqs(fasta=current, count=current, summary=current, start=1968, end=11550, maxhomop=8)
```

Check to see if our sequences overlap the same alignment coordinates
```
summary.seqs(fasta=current, count=current)
```

**Compare your output with the expected below**
```
                Start   End     NBases  Ambigs  Polymer NumSeqs
Minimum:        1965    11550   250     0       3       1
2.5%-tile:      1968    11550   252     0       3       3217
25%-tile:       1968    11550   252     0       4       32164
Median:         1968    11550   252     0       4       64328
75%-tile:       1968    11550   253     0       5       96492
97.5%-tile:     1968    11550   253     0       6       125439
Maximum:        1968    13400   270     0       8       128655
Mean:   1967    11550   252     0       4
# of unique seqs:       16298
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
unique.seqs(fasta=current, count=current)
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
2.5%-tile:      1       376     252     0       3       2955
25%-tile:       1       376     252     0       4       29544
Median:         1       376     252     0       4       59087
75%-tile:       1       376     253     0       5       88630
97.5%-tile:     1       376     253     0       6       115219
Maximum:        1       376     256     0       8       118173
Mean:   1       376     252     0       4
# of unique seqs:       2489
total # of seqs:        118173
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

Also of note is that "unknown" only pops up as a classification if the classifier cannot classify your sequence to one of the domains. If you run `summary.seqs(fasta=current, count=current)`, you'll see that we now have 2469 unique sequences and a total of 118011 total sequences. This means about 20 of our sequences were in these various groups.

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
make.biom(shared=current, constaxonomy=current, metadata=./MiSeq_SOP/stability.metadata)
```

**Rename the file to something simpler**
```
rename.file(biom=current, new=TEAMNAME.biom)
```
