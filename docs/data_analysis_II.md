# Overview
Exploring Metadata with your OTU data.

You should have most of the basics firgured out, but to get you started, one last time...

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
```
setwd("~/BIOL491.2019.microbe")
theme_set(theme_bw())
```
# Import data for use for phyloseq
```
sharedfile = "BIOL491.2019.shared"
taxfile = "BIOL491.2019.taxonomy"
mothur_data <- import_mothur(mothur_shared_file = sharedfile, mothur_constaxonomy_file = taxfile)
```
# Import sample metadata
```
mapfile = "BIOL491.2019.metadata.csv"
map <- read.csv(mapfile)
map <- sample_data(map)
rownames(map) <- map$SampleID
mb <- merge_phyloseq(mothur_data, map)
colnames(tax_table(mb)) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")
```
Now we have a phyloseq object called mb. 

# To start today, focus on your Single Site Sample, below is an example for Site C

```
mysiteC <- subset_samples(mb, Site=="SiteC")
```



# Introductory stats

```
# Calculate
mydata_pcoa_bray <- ordinate(
  physeq = mydata, 
  method = "PCoA",
  distance = "bray"
)
# Plot
plot_ordination(
  physeq = mydata,
  ordination = mydata_pcoa_bray,
  title = "PCoA of mydata (bray)",
  color = "Plate"
) + 
  geom_point(aes(color = Plate), alpha = 0.7, size = 4)


## Permanova
Here is an example of how to run a permanova test using the adonis function in vegan. In this example we are testing the hypothesis that samples from the two different plates have different centroids

```
# Calculate bray curtis distance matrix
mysiteC_bray <- phyloseq::distance(mysiteC, method = "bray")

# make a data frame from the sample_data
sampledf <- data.frame(sample_data(mysiteC))

# Adonis test
adonis(mysiteC_bray ~ Plate, data = sampledf)
```

Example output
```
Call:
adonis(formula = mysiteA_bray ~ Plate, data = sampledf) 
Permutation: free
Number of permutations: 999

Terms added sequentially (first to last)
          Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)
Plate      1    0.2244 0.22440 0.78306 0.04169  0.665
Residuals 18    5.1583 0.28657         0.95831       
Total     19    5.3827                 1.00000
```
This output tells us that our adonis test is not significant (p = 0.665) so we cannot reject the null hypothesis that our samples from different plates have same centroid.

If we had a significant test, then it would be worth running a **Homogeneity of dispersion test**

```
beta <- betadisper(mysiteC_bray, sampledf$Plate)
permutest(beta)
```

Example  output
```
Permutation test for homogeneity of multivariate dispersions
Permutation: free
Number of permutations: 999

Response: Distances
          Df  Sum Sq   Mean Sq      F N.Perm Pr(>F)
Groups     1 0.00067 0.0006744 0.0257    999  0.866
Residuals 18 0.47307 0.0262816  
```

Additionally, our betadisper results are not significant, meaning we cannot reject the null hypothesis that our groups have the same dispersions. This means we can be more confident that our adonis result is a real result, and not due to differences in group dispersions

There is a lot more analysis that can be done here. We could test different grouping variables, or we could create a more complex permanova by testing a model that combines multiple variables. We'll get back to that later.

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

## Testing signifcance with two variables

Prepare data

```
mysiteBCD_bray <- phyloseq::distance(mysiteBCD, method = "bray")
sampledf <- data.frame(sample_data(mysiteBCD))
```

We can write a more complex formula as below (typical model formula such as Y ~ A + B*C)

```
adonis(mysiteBCD_bray ~ House + Site, data = sampledf)
```

Example output

```
adonis(formula = mysiteBCD_bray ~ House + Site, data = sampledf) 
          Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
House      3    1.4630 0.48765  2.4831 0.45750  0.001 ***
Site       2    0.5564 0.27821  1.4166 0.17400  0.143    
Residuals  6    1.1783 0.19639         0.36849           
Total     11    3.1977                 1.00000 
```

It appears that we can reject the null hypothesis that samples from different house are from the same centroid.

adonis adds the terms of formula sequentially, so it is worth comparing the two orders so that you can be more confident of your results.
```
adonis(mysiteBCD_bray ~ Site + House, data = sampledf)
```

Example output
```
adonis(formula = mysiteBCD_bray ~ Site + House, data = sampledf)
          Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)   
House      3    1.4630 0.48765  2.4831 0.45750  0.002 **
```

Again, House is significant, so we should move on the final test of homogeneity of dispersions

```
beta <- betadisper(mysiteBCD_bray, sampledf$House)
permutest(beta)
```

Example output

```
Permutation test for homogeneity of multivariate dispersions
Response: Distances
          Df  Sum Sq  Mean Sq      F N.Perm Pr(>F)
Groups     3 0.04076 0.013586 0.3013    999  0.878
Residuals  8 0.36068 0.045085
```

Not significant, so we can be more confident of our earlier results.

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

Now you have lots of code and your head should be full of lots of ideas. Go ahead and start testing.
