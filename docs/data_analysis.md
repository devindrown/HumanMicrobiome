# Data Analysis
Here you will start to explore your own dataset. You have already aquired some basic skills. You'll apply that knowledge to this new dataset.

# Loading the complete data set
To get you started, I have provided a script to load our class data set.

1. If you haven't already, create a new RStudio Project (File>New Project). If RStudio asks you to “Save Current Workspace”, you should select “Don’t Save”.  You then want to select Existing Directory. On the next window, set the Project working directory to: `~/BIOL491_2024`
2. In the Console you can enter `source('treasurechest/LoadData.R')`. You can open this file in your explorer to see all of the data you've loaded.
3. For today’s analysis, you want to create a new R Script (File>New File>R Script) to hold all of the code you’re writing. This will create an empty document in a new panel. You should go ahead and save this document (File>Save). You can name the file anything you want, but keep the title informative and without space (e.g. `house_analysis.R`). It’s important to end the file in `.R` so that Rstudio knows it’s an R script.

## Output
Now we have a number of phyloseq objects:

* `mb` contains the entire dataset with 187 samples including negative controls
* `mbQC` excludes the negative controls for a reduced 120 samples, also removes some low abundance ASVs
* `mb_dirty` excludes the negative controls, but includes some contaminating ASVs

# Create your data sets

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

# Explore diversity

Now you have two new datasets Run the following to explore you data.

1. Calculate number of reads You may back to the previous lab, [Phyloseq and R for analysis and visualization](phyloseq_analysis_visualization)
2. [Alpha Diversity](alpha_diversity_plot)
3. [Bar Plots of diversity at different scales](community_composition_plot)
4. Ordination plot (only complete for the Site data set). Hint: Check out the treasurechest.

**Show your instructor the set of figures on your house data before you move on to the site data.**

## Optional idea
When looking at your site data, you might want to consider these options for plotting your diversity

We can specify a sample variable on which to group/organize samples along the horizontal (x) axis. An experimentally meaningful categorical variable is usually a good choice (e.g. Plate, House, Source, Type)

For a single house, think about groups of sites
```
plot_richness(myhouse, x = "Source", measures="InvSimpson")
```
For a single Site, think about other ways of showing all of the data
```
plot_richness(mysite, x = "Plate", measures="InvSimpson")
```


# Sampling complex subsets

If you'd like to sample more than one site at a time or more than one house you can do that in the following way

Create some lists, each item is enclosed in double quotes `"` and separated by a comma `,`
```
mysitelist = c("SiteX","SiteY","SiteZ")
myhouselist = c("3a4c","4226","3f92","415e")
```
**Note, the above houses do not exist. You will have to pick some from the full data set**

Next, create a subset as before, but with some masking
```
mysiteL <- subset_samples(mb, ((Site %in% mysitelist) & (House %in% myhouselist)))
```
Copy this into a new container so that you code doesn't rely on this specific name
```
mydata <- mysiteL
```

## Ordination with two variables

With this more complete dataset, you can create an ordination plot, here we'll use an NMDS plot

```
#Ordinate
mydata_nmds_bray <- ordinate(
  physeq = mydata, 
  method = "NMDS", 
  distance = "bray"
)
```

Next, we want to plot our results, but we'll use symbols for the different sites and solors for the various houses

```
# Plot 

house_colors <- rainbow_hcl(length(unique(myhouselist)))

plot_ordination(
  physeq = mydata,
  ordination = mydata_nmds_bray,
  color = "House",
  shape = "Site",
  title = "NMDS of mydata bacterial Communities"
) + 
  scale_color_manual(values = house_colors
  ) +
  geom_point(aes(color = House), alpha = 0.7, size = 6) +
  geom_point(colour = "grey90", size = 1.5) 

```


## Plot alpha diversity

Heere we'll divide up the X axis by Site and then color each house differently

```
plot_richness(mydata, x = "Site", color = "House", measures="InvSimpson")
```

## Complex bar plots with this kind of data

Transform to relative abundances

```
relmydata = transform_sample_counts(mydata,function(x) 100 * x / sum(x))
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

# What's next

Now you have lots of code and your head should be full of lots of ideas. Next week, we'll move on to testing.
