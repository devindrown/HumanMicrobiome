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
```
# Importing your data
I've put three data files in a folder in your drive. These are two of the 'mothur' output files and the start of the class metadata file. Here is some code to set your working directory. 
```
# Set working directory
setwd("~/BIOL491.2019.microbe")

# Set plotting theme
theme_set(theme_bw())
```


 You should see the three files as below in that directory

```
BIOL491.2019.metadata.csv
BIOL491.2019.shared
BIOL491.2019.taxonomy
```

Next, as you learned before you want to import these for use for phyloseq

```
# Assign variables for imported data
sharedfile = "BIOL491.2019.shared"
taxfile = "BIOL491.2019.taxonomy"

# Import mothur data
mothur_data <- import_mothur(mothur_shared_file = sharedfile,
                             mothur_constaxonomy_file = taxfile)
```
<!-- COMMENT OUT THIS STEP
## Some simple QC

We know that there are a couple of OTUs that are contaminating our samples, so let's remove those before we proceed any further

```
# Remove bad OTUs
badOTU = c("Otu00001","Otu00002")
allTaxa = taxa_names(mothur_data)
allTaxa <- allTaxa[!(allTaxa %in% badOTU)]
mothur_data_clean = prune_taxa(allTaxa, mothur_data)

```
END OF COMMENT
-->
## Metadata

We collected the sample metadata too, next we need to import this.

```
# Import sample metadata
mapfile = "BIOL491.2019.metadata.csv"
map <- read.csv(mapfile)

# Convert this dataframe into phyloseq forma
map <- sample_data(map)

# Assign rownames to be Sample ID's
rownames(map) <- map$SampleID

# Merge mothurdata object with sample metadata
mb <- merge_phyloseq(mothur_data, map)
```

If you type `mb` you should see the following output

```
phyloseq-class experiment-level object
otu_table()   OTU Table:         [ 8208 taxa and 192 samples ]
sample_data() Sample Data:       [ 192 samples by 18 sample variables ]
tax_table()   Taxonomy Table:    [ 8208 taxa by 6 taxonomic ranks ]
```

Now we have a phyloseq object called mb. 

Next reset the column names of the taxonomy
```
colnames(tax_table(mb)) <- c("Kingdom", "Phylum", "Class", 
  "Order", "Family", "Genus")
```


# Create your samples

The complete dataset is too big to really look at all at once.

Let's start by looking at all the sites within a single home. Pick a house that you want to work with can create a smaller dataset. You houses are identified by the last 4 digits of the ID

```
myhouse <- subset_samples(mb, House=="ab8a")
```

While you're making datasets, pick a Site that you want to work with and create a smaller dataset. You can look in the metadata file to see how the labels are formatted (e.g. SiteA, SiteB)
```
mysiteC <- subset_samples(mb, Site=="SiteC")
```

# Explore diversity

Now you have two new datasets, go back to the previous lab, [Phyloseq and R for analysis and visualization](phyloseq_analysis_visualization), and run the following to explore you data.

1. Calculate number of reads
2. Alpha Diversity
3. Bar Plots of diversity at different scales
4. Ordination plot

**Show your instructor the set of figures on your house data before you move on to the site data.**

When looking at your site data, you might want to consider these options for plotting your diversity

We can specify a sample variable on which to group/organize samples along the horizontal (x) axis. An experimentally meaningful categorical variable is usually a good choice (e.g. Plate, House, Source, Type)

For a single house, think about groups of sites
```
plot_richness(myhouse, x = "Source", measures="Chao1")
plot_richness(myhouse, x = "Type", measures="Chao1")
```
For a single Site, think about other ways of showing all of the data
```
plot_richness(mysiteC, x = "Plate", measures="Chao1")
plot_richness(mysiteC, x = "House", measures="Chao1")
```


# Sampling complex subsets

If you'd like to sample more than one site at a time or more than one house you can do that in the following way

Create some lists, each time in enclosed in double quotes `"` and separated by a comma `,`
```
mysitelist = c("SiteB","SiteC","SiteD")
myhouselist = c("3a4c","4226","3f92","415e")
```

Next, create a subset as before, but with some masking
```
mysiteBCD <- subset_samples(mb, ((Site %in% mysitelist) & (House %in% myhouselist)))
```

## Ordination with two variables

With this more complete dataset, you can create an ordination plot, here we'll use an NMDS plot

```
#Ordinate
mysiteBCD_nmds <- ordinate(
  physeq = mysiteBCD, 
  method = "NMDS", 
  distance = "bray"
)
```

Next, we want to plot our results, but we'll use symbols for the different sites and solors for the various houses

```
# Plot 

house_colors <- rainbow_hcl(length(unique(myhouselist)))

plot_ordination(
  physeq = mysiteBCD,
  ordination = mysiteBCD_nmds,
  color = "House",
  shape = "Site",
  title = "NMDS of mysiteBCD bacterial Communities"
) + 
  scale_color_manual(values = house_colors
  ) +
  #  geom_line() +
  geom_point(aes(color = House), alpha = 0.7, size = 6) +
  geom_point(colour = "grey90", size = 1.5) 

```


## Plot alpha diversity

Heere we'll divide up the X axis by Site and then color each house differently

```
plot_richness(mysiteBCD, x = "Site", color = "House", measures="Chao1")
```

## Complex bar plots with this kind of data

Transform to relative abundances

```
relmydata = transform_sample_counts(mysiteBCD,function(x) 100 * x / sum(x))
```

You probably don't want to look at all of your data at once. Here we are looking at the Phylum level and filtering out anything less than 1%. You might want to do something else for your own dataset.

```
relmydata_phylum <- relmydata %>%
  tax_glom(taxrank = "Phylum") %>%                     # group at Phylum level
  psmelt() %>%                                         # Melt to long format
  filter(Abundance > 1) %>%                         # Filter out low abundance taxa
  arrange(Phylum)                                   # Sort data frame alphabetically by Phylum
```

Pick some colors based on the phylum data (you can do deeper if you choose)

```
phylum_colors <- diverge_hcl(length(unique(relmydata_phylum$Phylum)))
```

Plot **Sites** across the **X axis** and make a separate **Row** for each **house**

```
ggplot(relmydata_phylum, aes(x = Site, y = Abundance, fill = Phylum)) + 
  facet_grid(House~.) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = phylum_colors) +
  # Remove x axis title
  theme(axis.title.x = element_blank()) + 
  ylab("Relative Abundance (Phylum > 1%) \n") +
  theme(axis.text.x=element_text(angle=90,hjust=1)) +
  ggtitle("Composition, Phylum")
```

What if you want to compare in the other dimension? Try this:

```
ggplot(relmydata_phylum, aes(x = House, y = Abundance, fill = Phylum)) + 
  facet_grid(Site~.) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = phylum_colors) +
  # Remove x axis title
  theme(axis.title.x = element_blank()) + 
  ylab("Relative Abundance (Phylum > 1%) \n") +
  theme(axis.text.x=element_text(angle=90,hjust=1)) +
  ggtitle("Composition, Phylum")
```

Finally, explore differences at the Phylum level with this plot (note the `x = Phylum`)
```
ggplot(relmydata_phylum, aes(x = Phylum, y = Abundance, fill = Phylum)) + 
  facet_grid(House ~ Site) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = phylum_colors) +
  # Remove x axis title
  theme(axis.title.x = element_blank()) + 
  ylab("Relative Abundance (Phylum > 1%) \n") +
  theme(axis.text.x=element_text(angle=90,hjust=1)) +
  ggtitle("Composition, Phylum")
```

The **`facet_grid`** function controls the formatting as `facet_grid(ROW_variable ~ COLUMN_variable)`

# Advanced QC

First we want to remove the negative controls from our dataset. You can subsample your data and select only those samples that are not of the QC type
```
mydata <- subset_samples(mb, Type!="QC")
```

The phyloseq package includes functions for filtering, subsetting, and merging abundance data. In the following example, the data is first transformed to relative abundance, creating the new GPr object, which is then filtered such that only OTUs with a mean greater than 10^-5 are kept.
```
mbr  = transform_sample_counts(mydata, function(x) x / sum(x) )
mbfr = filter_taxa(mbr, function(x) mean(x) > 1e-5, TRUE)
```

This results in a highly-subsetted object, mbfr, removing the really rare OTUs.

Another method: Remove taxa not seen more than 3 times in at least 20% of the samples. This protects against an OTU with small mean & trivially large C.V.

```
GP = filter_taxa(mb, function(x) sum(x > 3) > (0.2*length(x)), TRUE)
```

Standardize abundances to the median sequencing depth
```
total = median(sample_sums(GP))
standf = function(x, t=total) round(t * (x / sum(x)))
gps = transform_sample_counts(GP, standf)
```

Filter the taxa using a cutoff of 3.0 for the Coefficient of Variation
```
gpsf = filter_taxa(gps, function(x) sd(x)/mean(x) > 3.0, TRUE)
```

# What's next

Now you have lots of code and your head should be full of lots of ideas. Next week, we'll move on to testing.
