# Overview
Here you will start to explore your own dataset.

You should have most of the basics firgured out, but to get you started

# Getting your workspace ready

```
#Load libraries
library(ggplot2)
library(vegan)
library(dplyr)
library(scales)
library(grid)
library(reshape2)
library(phyloseq)
library(colorspace)
library(ape)
```

Here is some code to create and change your working directory. Before you run, be sure to chang **MYNAME** to something specific.

```
# Set working directory
mainDir <- "C:/Users/microbiome/documents"
subDir <- "MYNAME"

ifelse(!dir.exists(file.path(mainDir, subDir)), dir.create(file.path(mainDir, subDir)), FALSE)
setwd(file.path(mainDir, subDir))

# Set plotting theme
theme_set(theme_bw())
```

# Importing your data

First, we will import the mothur shared file, consensus taxonomy file, and our sample metadata and store them in one phyloseq object. By storing all of our data structures together in one object we can easily interface between each of the structures. For example, as we will see later, we can use criteria in the sample metadata to select certain samples from the OTU table

```
# Assign variables for imported data
sharedfile = "demodataR.shared"
taxfile = "demodataR.taxonomy"

# Import mothur data
mothur_data <- import_mothur(mothur_shared_file = sharedfile,
                             mothur_constaxonomy_file = taxfile)
```

The sample metadata is just a basic csv with columns for sample attributes. Here is a preview of what the sample metadata looks like.

```
SampleID	House	Site	Type
ExtractionNEGA	NONE	NEG	QC
ExtractionNEGB	NONE	NEG	QC
PCRNEGSetA1	NONE	NEG	QC
PCRNEGSetA2	NONE	NEG	QC
PCRNEGSetB1	NONE	NEG	QC
PCRNEGSetB2	NONE	NEG	QC
```

As you can see, there is one column called SampleID with the names of each of the samples. The remaining columns contain information on the sampling conditions related to each sample. The only formatting required to merge the sample data into a phyloseq object is that the rownames must match the sample names in your shared and taxonomy files.

```
mapfile = "demodataR.csv"
# Import sample metadata
map <- read.csv(mapfile)

# Convert this dataframe into phyloseq forma
map <- sample_data(map)

# Assign rownames to be Sample ID's
rownames(map) <- map$SampleID
```

We need to merge our metadata into our phyloseq object.
```
# Merge mothurdata object with sample metadata
moth_merge <- merge_phyloseq(mothur_data, map)
```

If you type `moth_merge` you should see the following output

```
> moth_merge
phyloseq-class experiment-level object
otu_table()   OTU Table:         [ 213 taxa and 6 samples ]
sample_data() Sample Data:       [ 6 samples by 4 sample variables ]
tax_table()   Taxonomy Table:    [ 213 taxa by 6 taxonomic ranks ]
```

Now we have a phyloseq object called moth.merge. 

Before we move on with analysis, we need to do some basic reformatting and filtering.

What are the column names of our taxonomy file?

```
colnames(tax_table(moth_merge))
## [1] "Rank1" "Rank2" "Rank3" "Rank4" "Rank5" "Rank6"
```
These taxonomy names are not helpful, so let’s rename them
```
colnames(tax_table(moth_merge)) <- c("Kingdom", "Phylum", "Class", 
  "Order", "Family", "Genus")
```

**Let's make a tree of our OTUs in case we want to use that later**
```
random_tree = rtree(ntaxa(moth_merge), rooted=TRUE, tip.label=taxa_names(moth_merge))
```

**Finally put your data in a new container**
```
mydata <- merge_phyloseq(moth_merge, random_tree)
```
