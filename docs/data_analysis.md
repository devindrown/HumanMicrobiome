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
library(colorspace)
library(phyloseq)
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

I've put four four data files on the Lab/Project Analysis section of Blackboard. You can download the mothur output of our class data there. You should end up with

```
microbiome2017.biom
microbiome2017.csv
microbiome2017.shared
microbiome2017.taxonomy
```
You want to make sure these end up in the directory you created in the previous steps.

Next, as you learned before you want to import these for use for phyloseq

```
# Assign variables for imported data
sharedfile = "microbiome2017.shared"
taxfile = "microbiome2017.taxonomy"

# Import mothur data
mothur_data <- import_mothur(mothur_shared_file = sharedfile,
                             mothur_constaxonomy_file = taxfile)
```

We collected the sample metadata too, next we need to import this.

```
# Import sample metadata
mapfile = "microbiome2017.csv"
map <- read.csv(mapfile)

# Convert this dataframe into phyloseq forma
map <- sample_data(map)

# Assign rownames to be Sample ID's
rownames(map) <- map$SampleID

# Merge mothurdata object with sample metadata
mb2017 <- merge_phyloseq(mothur_data, map)
```

If you type `mb2017` you should see the following output

```
phyloseq-class experiment-level object
otu_table()   OTU Table:         [ 12951 taxa and 128 samples ]
sample_data() Sample Data:       [ 128 samples by 5 sample variables ]
tax_table()   Taxonomy Table:    [ 12951 taxa by 6 taxonomic ranks ]
```

Now we have a phyloseq object called mb2017. 

Next reset the column names of the taxonomy
```
colnames(tax_table(mb2017)) <- c("Kingdom", "Phylum", "Class", 
  "Order", "Family", "Genus")
```

# Create your samples

The complete dataset is too big to really look at all at once.

Pick a Site that you want to work with can create a smaller dataset. You can look in the metadata file to see how the labels are formatted (e.g. SiteA, SiteB)
```
mysiteA <- subset_samples(mb2017, Site=="SiteC")
```

While you're making datasets, pick a house that you want to work with can create a smaller dataset. You houses are identified by the last 4 digits of the ID

```
myhouse <- subset_samples(mb2017, House=="5676")
```

Now you have two new datasets, go back to the [previous lab](phyloseq_analysis_visualization) and run the following to explore you data. **Show your instructor the set of figures on your house data before you move on to the site data.**

1. Calculate number of reads
2. Alpha Diversity
3. Bar Plots of diversity at different scales
4. Ordination plot


