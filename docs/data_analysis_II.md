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

As a shortcut, you can copy your new dataset into a generic container. This will allow you to reuse your visualization code for multiple datasets.
```
mydata <- mysiteC
```


# Introductory stats

## Plot
We'll look at all of our particular site, but color by Sequencing `Plate`
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
```

## Permanova
Here is an example of how to run a permanova test using the adonis function in vegan. In this example we are testing the hypothesis that samples from the two different plates have different centroids

```
# Calculate bray curtis distance matrix
mydata_bray <- phyloseq::distance(mydata, method = "bray")

# make a data frame from the sample_data
sampledf <- data.frame(sample_data(mydata))

# Adonis test
adonis(mydata_bray ~ Plate, data = sampledf)
```

Example output
```
Call:
adonis(formula = mydata_bray ~ Plate, data = sampledf) 
Permutation: free
Number of permutations: 999

Terms added sequentially (first to last)
          Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)
Plate      1    0.2964 0.29644 0.97378 0.05132   0.43
Residuals 18    5.4795 0.30442         0.94868       
Total     19    5.7760                 1.00000
```
This output tells us that our adonis test is not significant (p = 0.43) so we cannot reject the null hypothesis that our samples from different plates have same centroid.

If we had a significant test, then it would be worth running a **Homogeneity of dispersion test** Go ahead and run it now

```
beta <- betadisper(mydata_bray, sampledf$Plate)
permutest(beta)
```

Example  output
```
Permutation test for homogeneity of multivariate dispersions
Permutation: free
Number of permutations: 999

Response: Distances
          Df  Sum Sq   Mean Sq      F N.Perm Pr(>F)
Groups     1 0.00291 0.0029126 0.1134    999  0.711
Residuals 18 0.46214 0.0256745  
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

If you need to find a list of house IDs, then you can use the following command
```
print(unique(map$House))
```
The same command works for $Site


Next, create a subset as before, but with some masking
```
mycomplexdata <- subset_samples(mb, ((Site %in% mysitelist) & (House %in% myhouselist)))
```

## Ordination with two variables

With this more complete dataset, you can create an ordination plot, here we'll use an NMDS plot

```
#Ordinate
mycomplexdata_nmds <- ordinate(
  physeq = mycomplexdata, 
  method = "NMDS", 
  distance = "bray"
)

```

Next, we want to plot our results, but we'll use symbols for the different sites and solors for the various houses

```
# Plot 

house_colors <- rainbow_hcl(length(unique(myhouselist)))

plot_ordination(
  physeq = mycomplexdata,
  ordination = mycomplexdata_nmds,
  color = "House",
  shape = "Site",
  title = "NMDS of mycomplexdata bacterial Communities"
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
mycomplexdata_bray <- phyloseq::distance(mycomplexdata, method = "bray")
sampledf <- data.frame(sample_data(mycomplexdata))
```

We can write a more complex formula as below (typical model formula such as Y ~ A + B*C)

```
adonis(mycomplexdata_bray ~ House + Site, data = sampledf)
```

Example output

```
adonis(formula = mycomplexdata_bray ~ House + Site, data = sampledf) 
          Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)  
House      3    1.9064 0.63545  2.5398 0.47278  0.014 *
Site       2    0.6247 0.31234  1.2484 0.15492  0.249  
Residuals  6    1.5012 0.25019         0.37229         
Total     11    4.0322                 1.00000
```

It appears that we can reject the null hypothesis that samples from different house are from the same centroid.

adonis adds the terms of formula sequentially, so it is worth comparing the two orders so that you can be more confident of your results.
```
adonis(mycomplexdata_bray ~ Site + House, data = sampledf)
```

Example output
```
adonis(formula = mycomplexdata_bray ~ Site + House, data = sampledf)
          Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)   
House      3    1.9064 0.63545  2.5398 0.47278  0.016 *
```

Again, House is significant, so we should move on the final test of homogeneity of dispersions

```
beta <- betadisper(mycomplexdata_bray, sampledf$House)
permutest(beta)
```

Example output

```
Permutation test for homogeneity of multivariate dispersions
Response: Distances
          Df  Sum Sq  Mean Sq      F N.Perm Pr(>F)
Groups     3 0.01769 0.005896 0.0946    999  0.852
Residuals  8 0.49889 0.062361
```

Not significant, so we can be more confident of our earlier results.

## Plot alpha diversity

Here we'll divide up the X axis by Site and then color each house differently

```
plot_richness(mycomplexdata, x = "Site", color = "House", measures="Chao1")
```

## Complex bar plots with this kind of data

Transform to relative abundances

```
relmydata = transform_sample_counts(mycomplexdata,function(x) 100 * x / sum(x))
```

You probably don't want to look at all of your data at once. Here we are looking at the Order level and filtering out anything less than 5%. You might want to do something else for your own dataset.

```
relmydata_order <- relmydata %>%
  tax_glom(taxrank = "Order") %>%
  psmelt() %>%
  filter(Abundance > 5) %>%
  arrange(Order)
```

Pick some colors based on the Order data (you can do deeper if you choose). You can use your old code or if you'd like to explore some other colors, try this code, and then look at the names.
```
hcl_palettes(plot = TRUE)
```
If you like a palette under `Qualitative`, then you can pick a new color palette and swap out `diverge_hcl` with the alternative name `qualitative_hcl` and specify the palette
```
order_colors <- qualitative_hcl(length(unique(relmydata_order$Order)), palette = "Dark 3")
```
Plot **Sites** across the **X axis** and make a separate **Row** for each **house**

```
ggplot(relmydata_order, aes(x = Site, y = Abundance, fill = Order)) + 
  facet_grid(House~.) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = order_colors) +
  # Remove x axis title
  theme(axis.title.x = element_blank()) + 
  ylab("Relative Abundance (Order > 5%) \n") +
  theme(axis.text.x=element_text(angle=90,hjust=1)) +
  ggtitle("Composition, Order")
```

What if you want to compare in the other dimension? Try this:

```
ggplot(relmydata_order, aes(x = House, y = Abundance, fill = Order)) + 
  facet_grid(Site~.) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = order_colors) +
  # Remove x axis title
  theme(axis.title.x = element_blank()) + 
  ylab("Relative Abundance (Order > 5%) \n") +
  theme(axis.text.x=element_text(angle=90,hjust=1)) +
  ggtitle("Composition, Order")
```

Finally, explore differences at the Order level with this plot (note the `x = Order`)
```
ggplot(relmydata_order, aes(x = Order, y = Abundance, fill = Order)) + 
  facet_grid(House ~ Site) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = order_colors) +
  # Remove x axis title
  theme(axis.title.x = element_blank()) + 
  ylab("Relative Abundance (Order > 5%) \n") +
  theme(axis.text.x=element_text(angle=90,hjust=1)) +
  ggtitle("Composition, Order")
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
