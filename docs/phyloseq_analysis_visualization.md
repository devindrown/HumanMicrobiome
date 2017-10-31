# Overview
Here we will analyze out small dataset in R and produce some visualizations


Most of these instructions are modified from:
[Denef lab howto](http://deneflab.github.io/MicrobeMiseq/demos/mothur_2_phyloseq.html)

# Installing lots of Packages

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
```

```
# Set working directory
setwd("C:/Users/Microbiome/documents/demodataR")

# Source code files
# miseqR.R can be found in this repository
source("setwd("C:/Users/Microbiome/documents/demodataR/miseqR.R")

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
**Finally put your data in a new container**
```
mydata <- moth_merge
```

# Basic info and plots

**Calculate the number of reads per sample**

```
sample_sums(mydata)
```

Output should look like:
```
ExtractionNEGA ExtractionNEGB    PCRNEGSetA1    PCRNEGSetA2    PCRNEGSetB1    PCRNEGSetB2
          5108          71128          13710           4192          53655          12605
```
**Plot it**
```
plot_bar(mydata)
```

# Alpha Diversity

Estimating alpha diversity of microbial communities is problematic no matter what you do. My best stab at it is to subsample the libraries with replacement to estimate the species abundance of the real population while standardizing sampling effort.

```
min_lib <- min(sample_sums(mydata))
```
We will subsample to 4192, the minimum number of reads. We will repeat this 100 times and average the diversity estimates from each trial.

**Initialize matrices to store richness and evenness estimates**
```
nsamp = nsamples(mydata)
trials = 100

richness <- matrix(nrow = nsamp, ncol = trials)
row.names(richness) <- sample_names(mydata)

evenness <- matrix(nrow = nsamp, ncol = trials)
row.names(evenness) <- sample_names(mydata)

diversity <- matrix(nrow = nsamp, ncol = trials)
row.names(diversity) <- sample_names(mydata)
```

**It is always important to set a seed when you subsample so your result is replicable**
```
set.seed(3)
```
The create a loop to do all the subsampling
```
for (i in 1:100) {
  # Subsample
  r <- rarefy_even_depth(mydata, sample.size = min_lib, verbose = FALSE, replace = TRUE)
  
  # Calculate richness
  rich <- as.numeric(as.matrix(estimate_richness(r, measures = "Observed")))
  richness[ ,i] <- rich
  
  # Calculate evenness
  even <- as.numeric(as.matrix(estimate_richness(r, measures = "InvSimpson")))
  evenness[ ,i] <- even
  
}
```

Let’s calculate the mean and standard deviation per sample for observed richness and inverse simpson’s index and store those values in a dataframe.
```
# Create a new dataframe to hold the means and standard deviations of richness estimates
SampleID <- row.names(richness)
mean <- apply(richness, 1, mean)
sd <- apply(richness, 1, sd)
measure <- rep("Richness", nsamp)
rich_stats <- data.frame(SampleID, mean, sd, measure)

# Create a new dataframe to hold the means and standard deviations of evenness estimates
SampleID <- row.names(evenness)
mean <- apply(evenness, 1, mean)
sd <- apply(evenness, 1, sd)
measure <- rep("Inverse Simpson", nsamp)
even_stats <- data.frame(SampleID, mean, sd, measure)
```

**Now we will combine our estimates for richness and evenness into one dataframe**
```
alpha <- rbind(rich_stats, even_stats)
```

Let’s add the sample metadata into this dataframe using the `merge()` command
```
s <- data.frame(sample_data(mydata))
alphadiv <- merge(alpha, s, by = "SampleID") 
```

**Finally, we will plot the two alpha diversity measures using a facet**

```
ggplot(alphadiv, aes(x = SampleID, y = mean)) +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), colour="black", width=.1) +
  geom_point(size = 2) +
  facet_wrap(~measure, ncol = 1, scales = "free")
```
**Phyloseq also can do this for you**
```
plot_richness(mydata, measures="Chao1")
```
