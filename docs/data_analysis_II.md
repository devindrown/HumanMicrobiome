# Data Analysis Part II

Exploring Metadata with your OTU data.

# Loading the complete data set
To get you started, I have provided a script to load a larger data set including 383 samples spanning two years of data collection

In most cases, your data should still be available in RStudio. However, if you need to reload your project, you can navigate to the `BIOL491.combined.microbe` folder and then click on the *project file* `BIOL491.combined.microbe.Rproj`.

If you want to get a clean start, you can use the broom tool to clear objects from the workspace. After that you can enter `source('~/BIOL491.combined.microbe/BIOL491.combined.LoadData.R')`.

## Output
Now we have a number of phyloseq objects:

* `mb` contains the entire dataset with 383 samples including negative controls
* `mbQC` excludes the negative controls for a reduced 239 samples
* `mb_dirty` includes the entire dataset along with some contaminating OTUs

# Create your data sets, a REVIEW from last week

The complete data set is too big to really look at all at once.

Pick a house that you want to work with can create a smaller data set. You houses are identified by the last 4 digits of the ID. You can get a list of the included houses by typing `unique(map$House)`. Let's start by looking at all the sites within a single home. You can use the code below to put all of the samples from a single house (e.g. `ab8a`) into a container (`myhouse`)

```
myhouse <- subset_samples(mbQC, House=="ab8a")
```

While you're making data sets, pick a Site that you want to work with and create a smaller dataset. You can look in the metadata file to see how the labels are formatted (e.g. SiteA, SiteB)
```
mysite <- subset_samples(mb, Site=="SiteC")
```

**HINT** The code you used last week relied on your dataset being in a container called `mydata`. You can copy your own dataset into that temporary container with this short command `mydata <- myhouse` or `mydata <- mysite`.

# Introductory stats

## Plot the data
We'll look at all of our particular site, but color by sequencing `Plate`. In this case, Plate was just an arbritrary assignment during the library preparation. We do not expect Plate to have a signficiant impact on our sequencing output.
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

## Test for differences using Permanova
Here is an example of how to run a permanova test using the `adonis` function in vegan. In this example we are testing the hypothesis that samples from the two different plates have different centroids

```
# Calculate a distance matrix using Bray Curtis distances
mydata_distance <- phyloseq::distance(mydata, method = "bray")

# Create a data frame from the sample_data
sampledf <- data.frame(sample_data(mydata))

# Run the Adonis test
adonis(mydata_distance ~ Plate, data = sampledf)
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
This output tells us that our adonis test is not significant (`p = 0.43`). We cannot reject the null hypothesis that our samples from different plates have same centroid.

If we had a significant test, then it would be worth running a **Homogeneity of dispersion** test. Go ahead and run it now

```
beta <- betadisper(mydata_distance, sampledf$Plate)
permutest(beta)
```

Example output
```
Permutation test for homogeneity of multivariate dispersions
Permutation: free
Number of permutations: 999

Response: Distances
          Df  Sum Sq   Mean Sq      F N.Perm Pr(>F)
Groups     1 0.00291 0.0029126 0.1134    999  0.711
Residuals 18 0.46214 0.0256745  
```

Additionally, our betadispersion results are not significant, meaning we cannot reject the null hypothesis that our groups have the same dispersions. We can be more confident that our adonis result is a real result, and not due to differences in group dispersions.

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
The same command works for `$Site`


Next, create a subset as before, but with some masking
```
mycomplexdata <- subset_samples(mb, ((Site %in% mysitelist) & (House %in% myhouselist)))
```

If you copy this container into `mydata`, then you can reuse some of your previous graphing code.
```
mydata <- mycomplexdata
```

## Ordination with two variables

With this more complete dataset, you can create an ordination plot, here we'll use an PCoA plot

```
# Calculate distances
mydata_ord <- ordinate(
  physeq = mydata, 
  method = "PCoA",
  distance = "bray"
)
```

Next, we want to plot our results, but we'll use symbols for the different sites and colors for the various houses

```
# Get a list of colors
house_colors <- rainbow_hcl(length(unique(myhouselist)))

# Plot 
plot_ordination(
  physeq = mydata,
  ordination = mydata_ord,
  color = "House",
  shape = "Site",
  title = "PCoA of mycomplexdata bacterial Communities"
) + 
  scale_color_manual(values = house_colors) +
  geom_point(aes(color = House), alpha = 0.7, size = 6)
```

## Testing signifcance with two variables

Prepare data, this is the same as you did previously
```
# Calculate a distance matrix using Bray Curtis distances
mydata_distance <- phyloseq::distance(mydata, method = "bray")
# Create a data frame from the sample_data
sampledf <- data.frame(sample_data(mydata))
```

We can write a more complex formula as below (typical model formula such as `Y ~ A + B`)

```
adonis(mydata_distance ~ House + Site, data = sampledf)
```

Example output

```
adonis(formula = mycomplexdata_bray ~ House + Site, data = sampledf) 
          Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)  
House      3    1.9108 0.63695  1.9422 0.40882  0.006 **
Site       2    0.7955 0.39775  1.2128 0.17019  0.187   
Residuals  6    1.9677 0.32796         0.42099         
Total     11    4.0322                 1.00000
```

It appears that we can reject the null hypothesis that samples from different houses are from the same centroid (`p = 0.006`)

**adonis** adds the terms of formula sequentially, so it is worth comparing the two orders so that you can be more confident of your results.
```
adonis(mydata_distance ~ Site + House, data = sampledf)
```

Example output
```
adonis(formula = mydata_distance ~ Site + House, data = sampledf)
          Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)   
Site       2    0.7955 0.39775  1.2128 0.17019  0.198   
House      3    1.9108 0.63695  1.9422 0.40882  0.004 **
Residuals  6    1.9677 0.32796         0.42099          
Total     11    4.6741                 1.00000       
```

Again, House is significant (`p = 0.004`), so we should move on the final test of homogeneity of dispersions and specify `House` in the dataframe.

```
beta <- betadisper(mydata_distance, sampledf$House)
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

## Alpha diversity with two variables

We can group the data along the X axis by Site and then use color to distinguish houses.

```
plot_richness(mycomplexdata, x = "Site", color = "House", measures="Chao1")
```
*If you wanted to look at differences among house, how might you plot the data?*

## Community composition bar plots with two variables

You probably don't want to look at all of your data at once. Here we are looking at the Order level and filtering out anything less than 5%. You might want to do something else for your own dataset.

```
myTaxLevel <- "Order"
myFilter <- 5
myYaxis <- paste("Relative Abundance (", myTaxLevel, " > ", myFilter, "%) \n")
```
Transform to relative abundances
```
relmydata = transform_sample_counts(mydata,function(x) 100 * x / sum(x))

relmydata_grouped <- relmydata %>%
  tax_glom(taxrank = myTaxLevel) %>%        # group at your Taxonomic level
  psmelt() %>%                              # Melt to long format
  filter(Abundance > myFilter) %>%          # Filter out low abundance taxa
  arrange(myTaxLevel)                       # Sort data frame alphabetically by your Taxonomic level
relmydata_grouped_clean <- subset(relmydata_grouped, relmydata_grouped[[myTaxLevel]] != "Bacteria_unclassified")
```

Pick some colors based on the Order data (you can do deeper if you choose). You can use your old code or if you'd like to explore some other colors, try this code, and then look at the names.
```
hcl_palettes(plot = TRUE)
```
If you like a palette under `Qualitative`, then you can pick a new color palette and swap out `diverge_hcl` with the alternative name `qualitative_hcl` and specify the palette
```
my_colors <- qualitative_hcl(length(unique(relmydata_grouped_clean[[myTaxLevel]])), palette = "Dark 3")
```
Plot **Sites** across the **X axis** and make a separate **Row** for each **house**

```
ggplot(relmydata_grouped_clean, aes_string(x = "Site", y = "Abundance", fill = myTaxLevel)) + 
  facet_grid(House~.) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = my_colors) +
  theme(axis.title.x = element_blank()) + 
  ylab(myYaxis) +
  theme_bw()+
  theme(axis.text.x=element_text(angle=90,hjust=1)) +
  ggtitle("Community Composition")
```

What if you want to compare in the other dimension? Try this:

```
ggplot(relmydata_grouped_clean, aes_string(x = "House", y = "Abundance", fill = myTaxLevel)) + 
  facet_grid(Site~.) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = my_colors) +
  theme(axis.title.x = element_blank()) + 
  ylab(myYaxis) +
  theme_bw()+
  theme(axis.text.x=element_text(angle=90,hjust=1)) +
  ggtitle("Community Composition")
```

Finally, explore differences at the Order level with this plot (note the `x = Order`)
```
ggplot(relmydata_grouped_clean, aes_string(x = myTaxLevel, y = "Abundance", fill = myTaxLevel)) + 
  facet_grid(House ~ Site) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = my_colors) +
  theme(axis.title.x = element_blank()) + 
  ylab(myYaxis) +
  theme_bw()+
  theme(axis.text.x=element_text(angle=90,hjust=1)) +
  ggtitle("Community Composition")
```

The **`facet_grid`** function controls the formatting as `facet_grid(ROW_variable ~ COLUMN_variable)`. You can explore your data in many different way.

# Advanced QC, redux

In case you did not comlete this in the previous lab. Here are the instructions again for advance quality control. 

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

**Now what?** You could use this new data set and explore the impact of reducing some of the noise in your sequence data.
