# Data Analysis Part II

Exploring Statistics and Complex Datasets

# Loading the current year data set, a review from last week

1. You should still be in the `/BIOL491_2024` RStuido Project.
2. Get a clean start, Use the broom tool to clear objects from the workspace. After that, i the Console, you can enter
```source('treasurechest/LoadData.R')```
3. For today’s analysis, you want to create a new R Script (File>New File>R Script) to hold all of the code you’re writing. This will create an empty document in a new panel. You should go ahead and save this document (File>Save). You can name the file anything you want, but keep the title informative and without space (e.g. `house_stats.R`). It’s important to end the file in `.R` so that Rstudio knows it’s an R script.

## Output
Now we have a number of phyloseq objects:

* `mb` contains the entire dataset with 187 samples including negative controls
* `mbQC` excludes the negative controls for a reduced 120 samples, also removes some low abundance ASVs
* `mb_dirty` excludes the negative controls, but includes some contaminating ASVs

# Create your data sets, a REVIEW from last week

The complete data set is too big to really look at all at once.

Pick a house that you want to work with can create a smaller data set. You houses are identified by the last 4 digits of the ID. You can get a list of the included houses by typing `levels(mb@sam_data$House)`. Let's start by looking at all the sites within a single home. You can use the code below to put all of the samples from a single house (e.g. `ab8a`) into a container (`myhouse`)

```
myhouse <- subset_samples(mbQC, House=="ab8a")
```

While you're making data sets, pick a Site that you want to work with and create a smaller dataset. You can look in the metadata file `view(mb@sam_data)` to see how the labels are formatted (e.g. SiteA, SiteB)
```
mysite <- subset_samples(mb, Site=="SiteZ")
```

**HINT** The code you used last week relied on your dataset being in a container called `mydata`. You can copy your own dataset into that temporary container with this short command `mydata <- myhouse` or `mydata <- mysite`.

# Introductory stats

## Plot the data
For this example, we'll look at our site data set, but color by sequencing `Plate`. In this case, Plate was just an arbritrary assignment during the library preparation. We do not expect Plate to have a signficiant impact on our sequencing output.

Start by copying the your site data set into a generic container `mydata`.
```
mydata <- mysite`
```
Now let's make an ordination to visualize the information
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

## Test for differences in community composition using Permanova
Here is an example of how to run a permanova test using the `adonis2` function in vegan. In this example we are testing the hypothesis that samples from the two different plates have different centroids

```
# Perform PERMANOVA analysis using the adonis2 function
# Formula: Bray-Curtis distance ~ Plate
# Data: Sample data from the phyloseq object
adonis2(
  phyloseq::distance(mydata, method = "bray") ~ Plate,  
  data = data.frame(mydata@sam_data)                  
)
```
What does the above code do?
* `adonis2()`: This function performs permutational multivariate analysis of variance (PERMANOVA) using the `vegan` package in R. It is used to test for significant differences in community composition between groups.
* `phyloseq::distance(mydata, method = "bray")`: This calculates the Bray-Curtis dissimilarity matrix between samples in your mydata object (a phyloseq object). Bray-Curtis is a common distance metric used for ecological data.
* `~ Plate`: This is the formula specifying the model. It means you are testing for differences in community composition explained by the variable `Plate` (a factor representing different groups or treatments).
* data = data.frame(mydata@sam_data): This specifies the data frame containing the explanatory variable (`Plate`). It extracts the sample data from your phyloseq object (`mydata@sam_data`) and converts it into a data frame.

Example output
```
Permutation test for adonis under reduced model
Permutation: free
Number of permutations: 999

adonis2(formula = phyloseq::distance(mydata, method = "bray") ~ Plate, data = data.frame(mydata@sam_data))
         Df SumOfSqs      R2      F Pr(>F)
Model     1   0.3032 0.04541 0.8564  0.712
Residual 18   6.3725 0.95459              
Total    19   6.6757 1.00000
```
This output tells us that our adonis test is not significant (`p > 0.05`). We cannot reject the null hypothesis that our samples from different plates have same centroid.

## Test for differences in dispersion

If we had a significant test, then it would be worth running a **Homogeneity of dispersion** test. Go ahead and run it now.

This code performs a permutation test to assess the homogeneity of multivariate dispersions between groups defined by a variable. This is important to check before conducting a PERMANOVA analysis, as significant differences in dispersion can influence the results of PERMANOVA.
```
# Perform a permutation test for homogeneity of multivariate dispersions
# Calculate Bray-Curtis distances
# Grouping factor: "Plate" variable
permutest(
  betadisper(phyloseq::distance(mydata, method = "bray"),
             mydata@sam_data$Plate
             )
  )
```
What does the above code do?
* `permutest()`: This function performs a permutation test on an object of class `betadisper`. It tests the null hypothesis that the dispersions (variances) of the groups are equal.
* `betadisper()`: This function from the `vegan` package calculates the multivariate dispersions for each group. It takes two arguments:
    * A distance matrix: `phyloseq::distance(mydata, method = "bray")` calculates the Bray-Curtis distance matrix.
    * A grouping factor: `mydata@sam_data$Plate` specifies the groups based on the `Plate` variable.

Example output
```
Permutation test for homogeneity of multivariate dispersions
Permutation: free
Number of permutations: 999

Response: Distances
          Df   Sum Sq   Mean Sq      F N.Perm Pr(>F)
Groups     1 0.000755 0.0007545 0.0923    999  0.759
Residuals 18 0.147117 0.0081732 
```

Additionally, our betadispersion results are not significant, meaning we cannot reject the null hypothesis that our groups have the same dispersions. We can be more confident that our adonis result is a real result, and not due to differences in group dispersions.

There is a lot more analysis that can be done here. We could test different grouping variables, or we could create a more complex permanova by testing a model that combines multiple variables. We'll get back to that later.

# Sampling complex subsets

If you'd like to sample more than one site at a time or more than one house you can do that in the following way

Create some lists, each item is enclosed in double quotes `"` and separated by a comma `,`
```
mysitelist = c("SiteX","SiteY","SiteZ")
myhouselist = c("3a4c","4226","3f92","415e")
```
**Note, the above houses do not exist. You will have to pick some from the full data set**
If you need to find a list of house IDs, then you can use the following command
```
print(levels(mb@sam_data$House))
```
The same command works for `$Site`

Next, create a subset as before, but with some masking
```
mycomplexdata <- subset_samples(mb, ((Site %in% mysitelist) & (House %in% myhouselist)))
```
Copy this into a new container so that you code doesn't rely on this specific name
```
mydata <- mycomplexdata
```

## Ordination with two variables

With this more complete dataset, you can create an ordination plot, here we'll use an PCoA plot

```
# Calculate distances
mydata_pcoa_bray <- ordinate(
  physeq = mydata, 
  method = "PCoA",
  distance = "bray"
)
```

Next, we want to plot our results, but we'll use symbols for the different houses and colors for the various sites

```
# Get a list of colors
site_colors <- rainbow_hcl(length(unique(mysitelist)))

# Plot 
plot_ordination(
  physeq = mydata,
  ordination = mydata_pcoa_bray,
  color = "Site",
  shape = "House",
  title = "PCoA of mycomplexdata bacterial Communities"
) + 
  scale_color_manual(values = site_colors) +
  geom_point(aes(color = Site), alpha = 0.7, size = 6)
```

## Testing signifcance with two variables

We can write a more complex formula as below (typical model formula such as `Y ~ A + B`)
```
# Perform PERMANOVA analysis using the adonis2 function
# Formula: Bray-Curtis distance ~ House + Site
# Data: Sample data from the phyloseq object
  adonis2(
    formula  phyloseq::distance(mydata, method = "bray") ~ House + Site,  
    data = data.frame(mydata@sam_data)                  
  )
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

What is the impact of the quality control? Can you measure it? Remember you have been working with the quality controlled data. However, you have multiple completed data sets

* `mbQC` excludes the negative controls for a reduced 120 samples, also removes some low abundance ASVs
* `mb_dirty` excludes the negative controls, but includes some contaminating ASVs

**Now what?** You could compare the `mbQC` and `mb_dirty` versions of your house or site data set to impact of reducing some of the noise in your sequence data.
