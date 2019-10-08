# Overview
Here we will analyze our own data. We'll start with a small set.

# Getting started

* Confirm data files and databases

You should verify that you have the databases we used in the previous lab. If you are in doubt visit those instructions:

[UAF mothur MiSeq SOP](https://devindrown.github.io/HumanMicrobiome/mothur_MiSeq_SOP)

**Check the raw data folder**

 `rawdata_neg` confirm that you have **12 original data files**

```
2017a1PCRneg_S16_L001_R1_001.fastq
2017a1PCRneg_S16_L001_R2_001.fastq
2017a2PCRneg_S8_L001_R1_001.fastq
2017a2PCRneg_S8_L001_R2_001.fastq
2017aDNAneg_S32_L001_R1_001.fastq
2017aDNAneg_S32_L001_R2_001.fastq
2017b1PCRneg_S80_L001_R1_001.fastq
2017b1PCRneg_S80_L001_R2_001.fastq
2017b2PCRneg_S72_L001_R1_001.fastq
2017b2PCRneg_S72_L001_R2_001.fastq
2017bDNAneg_S96_L001_R1_001.fastq
2017bDNAneg_S96_L001_R2_001.fastq
2018aDNAneg_PlateA_G12_R1.fastq
2018aDNAneg_PlateA_G12_R2.fastq
2018aPCRneg_PlateA_H12_R1.fastq
2018aPCRneg_PlateA_H12_R2.fastq
2018bDNAneg_PlateB_G12_R1.fastq
2018bDNAneg_PlateB_G12_R2.fastq
2018bHOUSEDNAneg_PlateB_G10_R1.fastq
2018bHOUSEDNAneg_PlateB_G10_R2.fastq
2018bHOUSEPCRneg_PlateB_H10_R1.fastq
2018bHOUSEPCRneg_PlateB_H10_R2.fastq
2018bPCRneg_PlateB_H12_R1.fastq
2018bPCRneg_PlateB_H12_R2.fastq
```

## Starting the program

Now you are ready to start the local version of mothur, you can double click on the `mothur` file

**Setup a directory to put all of your output.**
```
set.dir(output=./negativeQC)
```
<!-- 
# Prepare the custom database

DO NOT COMPLETE THIS STEP

Like last time, we'll be aligning our sequencings to a currated set of 16s sequences to the Silva database. We need to trim our samples to just the v4 region that we actually amplified. We used the 515f and 806rb primers from the Earth Microbiome Project. We'll use the actual primers that we used for PCR to perform a digitial PCR on the SILVA database.

**Customize database to our region of interest**
```
pcr.seqs(fasta=./reference/silva.seed_v132.align, oligos=./reference/EMP_primers.oligos, keepdots=F)
```

**Rename the reference database**

Us the output from the previous command (the .pcr.align file) as input for the next command.
```
rename.file(input=LASTOUTPUT, new=silva.seed_v132.EMP.fasta)
```
**Trim the taxonomy database**
First list all the remaining sequences left in the FASTA database, then select only those in the Taxonomy database (~7816 taxa)
```
list.seqs(fasta=./negativeQC/silva.seed_v132.EMP.fasta)
get.seqs(accnos=./negativeQC\silva.seed_v132.EMP.accnos, taxonomy=./reference/silva.seed_v132.tax)
```

**Rename the taxonomy database**
```
rename.file(input=./negativeQC/silva.seed_v132.pick.tax, new=silva.seed_v132.EMP.tax)
```
-->
# Make contigs

First we need to take the paired-end read data from the MiSeq and join the read 1 and read 2 into a single sequences or contig. For example, you will see two files: ExtractionNEG-B-2017001057_S96_L001_R1_001.fastq and ExtractionNEG-B-2017001057_S96_L001_R2_001.fastq.. The first and all those with R1 correspond to read 1 while the second and all those with R2 correspond to the second or reverse read. These sequences are 250 bp and overlap in the V4 region of the 16S rRNA gene; this region is about 253 bp long.

Last time, you had a file preapred for you that contained all of the raw data files. This time around, you'll be making and editing that file from scratch. Mothur helps to get you started. You give the program a location (directory) and type of file (fastq).

```
make.file(inputdir=./rawdata_neg, type=fastq, prefix=neg.stability)
```

**Check the stability file**

You'll need to open the file in `notepad` and check the first column. Make sure each sample a short name without spaces or dashes. 

Remeber, the first column is the name othe sample. The second column is the name of the forward read for that sample and the third columns in the name of the reverse read for that sample.

**Combine our two sets of reads for each sample and then to combine the data from all of the samples**

```
make.contigs(inputdir=./rawdata_neg, file=./negativeQC/neg.stability.files, processors=8)
```

**Q1: Report how many sequences are in each group**

**Determine the distributions of these sequences**
```
summary.seqs(fasta=current)
```
**Q2: How many total sequences do you have? What's the median length?**

# Reducing sequencing errors

**Next, we want to get rid of some of the bad reads**

Note how long most of the reads are. We'll use that to trim out the reads that assembled backwards.

```
screen.seqs(fasta=current, group=current, maxambig=0, maxlength=311)
```

**Q3 How many sequences do you have now?** 

What command do you have to use?

Many of our sequences are duplicates of each other. It's computationally wasteful to align the same thing a bazillion times.

**Select only the unique sequences**

```
unique.seqs()
```

**Q4 How many unique sequences do you have?**

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
align.seqs(fasta=current, reference=./reference/silva.seed_v132.EMP.fasta)
```

**Check the results**

```
summary.seqs(count=current)
```
 You'll see that the bulk of the sequences start at position 27 and end at position 9609. Deviants from the mode positions are likely due to an insertion or deletion at the terminal ends of the alignments. 
 
 **Q5: What's the median sequence length after aligning to the reference?**

To make sure that everything overlaps the same region we'll re-run `screen.seqs`  

```
screen.seqs(fasta=current,count=current,summary=current,start=27, end=9609, maxhomop=8)
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

**Q6: How many unique sequences are left?**

**Further de-noise our sequences** using the `pre.cluster` command allowing for up to 2 differences between sequences.

```
pre.cluster(fasta=current, count=current, diffs=2)
```

**Q7: How many sequences are left? How many unique sequences?**

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

**Q8: How many sequences are left? How many unique sequences?**

Sometimes when we pick a primer set they will amplify other stuff such as 16S rRNA from chloroplasts, and mitochondria.

**Classify your sequences**

```
classify.seqs(fasta=current, count=current, reference=./reference/silva.seed_v132.EMP.fasta, taxonomy=./reference/silva.seed_v132.EMP.tax, cutoff=80)
```

Now that everything is classified we want to **remove our undesirables**

```
remove.lineage(fasta=current, count=current, taxonomy=current, taxon=Chloroplast-Mitochondria-unknown-Archaea-Eukaryota)
```

**Q9: How many sequences are left? How many unique sequences?**

**Create an updated taxonomy summary file that reflects these removals**
```
summary.tax(taxonomy=current, count=current)
```

This creates a new summary file with the undesirables removed. At this point we have curated our data as far as possible

# Generate OTUs

**Clustering sequences into OTUs**

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

# Generate a BIOM to Visualize

**Create a `.biom` file**
```
make.biom(shared=current,constaxonomy=current, metadata=./rawdata_neg/neg.stability.metadata)
```

**Rename the file to include a Team Name**
```
rename.file(biom=current, new=TEAMNAME.biom)
```
