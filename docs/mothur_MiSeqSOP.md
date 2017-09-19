Welcome class

From here on out you’ll be entering commands directly to the mothur program. The command prompt should look like this:
```
mothur>
```




```
set.dir(input=./MiSeq_SOP, output=./testrun)
```
# Reducing sequencing and PCR errors

Combine our two sets of reads for each sample and then to combine the data from all of the samples
```
make.contigs(file=stability.files, processors=8)
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
summary.seqs()
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

summary.seqs()

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

Output File Names:
/home/microbiome/mothur/testrun/stability.trim.contigs.good.unique.summary

```
## Preparing for aligning to database

Customize database to our region of intereste.
```
pcr.seqs(fasta=./silva.bacteria/silva.bacteria.fasta, start=11894, end=25319, keepdots=F)
```
OUTPUT
```
Output File Names:
/home/microbiome/mothur/testrun/silva.bacteria.pcr.fasta
```
rename.file(input=/home/microbiome/mothur/testrun/silva.bacteria.pcr.fasta, new=silva.v4.fasta) 


```



