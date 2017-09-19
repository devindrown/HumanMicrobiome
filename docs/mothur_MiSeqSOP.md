Welcome class

From here on out you’ll be entering commands directly to the mothur program. The command prompt should look like this:
```
mothur>
```

# Prepare database

```
set.dir(output=./testrun)
```


Customize database to our region of interes
```
pcr.seqs(fasta=./silva.bacteria/silva.bacteria.fasta, start=11894, end=25319, keepdots=F)
```
Output File Names:
/home/microbiome/mothur/testrun/silva.bacteria.pcr.fasta

Paying careful attention to the output file, use this in the below command <INPUT>

```
rename.file(input=/home/microbiome/mothur/testrun/silva.bacteria.pcr.fasta, new=silva.v4.fasta)
rename.file(input=<INPUT>, new=silva.v4.fasta)
```

# Reducing sequencing and PCR errors

Combine our two sets of reads for each sample and then to combine the data from all of the samples

```
make.contigs(inputdir=./MiSeq_SOP, file=./MiSeq_SOP/stability.files, processors=8)
```

You’ll probably get a [WARNING] message, but don’t worry.

Go back to the file explorer and open `stability.files`

What do you see?  The first column is the name of the sample. The second column is the name of the forward read for that sample and the third columns in the name of the reverse read for that sample.

Output
```
Group count:
F3D0    7793
F3D1    5869
F3D141  5958
F3D142  3183
F3D143  3178
F3D144  4827
F3D145  7377
F3D146  5021
F3D147  17070
F3D148  12405
F3D149  13083
F3D150  5509
F3D2    19620
F3D3    6758
F3D5    4448
F3D6    7989
F3D7    5129
F3D8    5294
F3D9    7070
Mock    4779

Total of all groups is 152360
```

```
Output File Names:
/home/microbiome/mothur/testrun/stability.trim.contigs.fasta
/home/microbiome/mothur/testrun/stability.trim.contigs.qual
/home/microbiome/mothur/testrun/stability.contigs.report
/home/microbiome/mothur/testrun/stability.scrap.contigs.fasta
/home/microbiome/mothur/testrun/stability.scrap.contigs.qual
/home/microbiome/mothur/testrun/stability.contigs.groups
```

Let's see what these sequences look like
```
summary.seqs(fasta=current)
```

OUTPUT
```
Using /home/microbiome/mothur/testrun/stability.trim.contigs.fasta as input file for the fasta parameter.

Using 8 processors.

                Start   End     NBases  Ambigs  Polymer NumSeqs
Minimum:        1       248     248     0       3       1
2.5%-tile:      1       252     252     0       3       3810
25%-tile:       1       252     252     0       4       38091
Median:         1       252     252     0       4       76181
75%-tile:       1       253     253     0       5       114271
97.5%-tile:     1       253     253     6       6       148552
Maximum:        1       502     502     249     243     152360
Mean:   1       252.811 252.811 0.70063 4.44854
# of Seqs:      152360

Output File Names:
/home/microbiome/mothur/testrun/stability.trim.contigs.summary
```

Verify that you have the same number of reads `# of Seqs:      152360`

Next, we want to get rid of some of the bad reads

```
screen.seqs(fasta=current, group=current, maxambig=0, maxlength=275)
```

It is a lot of work to keep typing in `fasta=example.fasta` etc. Try this: 
```
get.current()
```
mothur remembers your latest fasta file and group file as well as the number of processors you have, so you don’t have to type them in each time

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

Many of our sequences are duplicates of each other. Because it's computationally wasteful to align the same thing a bazillion times, we'll unique our sequences using

```
unique.seqs()
```

OUTPUT
```
Output File Names:
/home/microbiome/mothur/testrun/stability.trim.contigs.good.names
/home/microbiome/mothur/testrun/stability.trim.contigs.good.unique.fasta
```


Generate a table where the rows are the names of the unique sequences and the columns are the names of the groups. In this case, we’ll only have one group.

```
count.seqs()
```

OUTPUT
```
Total number of sequences: 128872

Output File Names:
/home/microbiome/mothur/testrun/stability.trim.contigs.good.count_table
```

```
summary.seqs(count=current)
```

OUTPUT
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

Now we need to align our sequences to the reference alignment. I’ve already trimmed the silva database to just the v4 region of 16s, so if you were using a different section, you’d need to do this from scratch. See the MiSeq SOP online for details. Here we’ll just use the shortcut.

```
align.seqs(fasta=current, reference=./testrun/silva.v4.fasta)
```

Check the results. Enter the following command all on one line:
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

OUTPUT
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


Filter the sequences to remove the overhangs at both ends. In addition, there are many columns in the alignment that only contain gap characters.
```
filter.seqs(fasta=current, vertical=T, trump=.)
```

OUTPUT
```
Length of filtered alignment: 376
Number of columns removed: 13049
Length of the original alignment: 13425
Number of sequences used to construct filter: 16298
```


Perhaps we’ve created some redundancy across our sequences by trimming the ends, we can re-run unique.seqs:
```
unique.seqs(fasta=current,count=current)
```

Further de-noise our sequences using the pre.cluster command allowing for up to 2 differences between sequences. This command will split the sequences by group and then sort them by abundance and go from most abundant to least and identify sequences that are within 2 nt of each other. If they are then they get merged. We generally favor allowing 1 difference for every 100 bp of sequence
```
pre.cluster(fasta=current, count=current, diffs=2)
```

# We have removed as much sequencing error as we can!
Check for chimeras
chimera.vsearch(fasta=current, count=current, dereplicate=t)

You still need to remove those sequences from the fasta file. Depending on the outcome of the chimera check, the following command may give an error. If there are no bad sequences to remove, then it will produce an error message.
```
remove.seqs(fasta=current, accnos=current, count=current)
```
Sometimes when we pick a primer set they will amplify other stuff such as 16S rRNA from chloroplasts, and mitochondria. 
```
classify.seqs(fasta=current, count=current, reference=./reference/trainset9_032012.pds.fasta, taxonomy=./reference/trainset9_032012.pds.tax, cutoff=80)
```
Now that everything is classified we want to remove our undesirables
```
remove.lineage(fasta=current, count=current, taxonomy=current, taxon=Chloroplast-Mitochondria-unknown-Archaea-Eukaryota)
```



