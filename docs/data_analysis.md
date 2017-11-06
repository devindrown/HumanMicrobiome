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

# Explore diversity

Now you have two new datasets, go back to the [previous lab](phyloseq_analysis_visualization) and run the following to explore you data.

1. Calculate number of reads
2. Alpha Diversity
3. Bar Plots of diversity at different scales
4. Ordination plot

**Show your instructor the set of figures on your house data before you move on to the site data.**

# Introductory stats

## Permanova
Here is an example of how to run a permanova test using the adonis function in vegan. In this example we are testing the hypothesis that samples from the two different plates have different centroids

```
# Calculate bray curtis distance matrix
mysiteA_bray <- phyloseq::distance(mysiteA, method = "bray")

# make a data frame from the sample_data
sampledf <- data.frame(sample_data(mysiteA))

# Adonis test
adonis(mysiteA_bray ~ Plate, data = sampledf)
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
beta <- betadisper(mysiteA_bray, sampledf$Plate)
permutest(beta)
``

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

# Sampling other sets

If you'd like to sample more than one site at a time or more than one house you can do that in the following way

Create some lists, each time in enclosed in double quotes `"` and separated by a comma `,`
```
mysitelist = c("SiteA","SiteC","SiteD")
myhouselist = c("3a4c","4226","3f92","415e")
```

Next, create a subset as before, but with some masking
```
mysiteACD <- subset_samples(mb2017, ((Site %in% mysitelist) & (House %in% myhouselist)))
```

## Ordination with two variables

With this more complete dataset, you can create an ordination plot, here we'll use an NMDS plot

```
#Ordinate
mysiteACD_pcoa <- ordinate(
  physeq = mysiteACD, 
  method = "NMDS", 
  distance = "bray"
)

```

Next, we want to plot, but we'll use symbols for the different sites and solors for the various houses

```
# Plot 

house_colors <- rainbow_hcl(length(unique(myhouselist)))

plot_ordination(
  physeq = mysiteACD,
  ordination = mysiteACD_pcoa,
  color = "House",
  shape = "Site",
  title = "NMDS of mysiteACD bacterial Communities"
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
mysiteACD_bray <- phyloseq::distance(mysiteACD, method = "bray")
sampledf <- data.frame(sample_data(mysiteACD))
```

We can write a more complex formula as below (typical model formula such as Y ~ A + B*C)

```
adonis(mysiteACD_bray ~ House + Site, data = sampledf)
```

adonis adds the terms of formula sequentially
```
adonis(mysiteACD_bray ~ Site + House, data = sampledf)
```



## Complex bar plots with this kind of data


Transform to relative abundances

```
relmydata = transform_sample_counts(mysiteACD,function(x) 100 * x / sum(x))
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
