# Overview
Here we will visualize our small dataset in R.


Most of these instructions are modified from:
[Denef lab howto](http://deneflab.github.io/MicrobeMiseq/demos/mothur_2_phyloseq.html)

# Getting your workspace ready

```
#Load libraries
library(tidyverse)
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

```
# Set working directory
setwd("~/demodataR")

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

The sample metadata is just a basic `.csv` with columns for sample attributes. Here is a preview of what the sample metadata looks like.

| SampleID     | year | type |
|--------------|------|------|
| 2017a1PCRneg | 2017 | PCR  |
| 2017a2PCRneg | 2017 | PCR  |
| 2017aDNAneg  | 2017 | DNA  |

As you can see, there is one column called SampleID with the names of each of the samples. The remaining columns contain information on the sampling conditions related to each sample. The only formatting required to merge the sample data into a phyloseq object is that the rownames must match the sample names in your shared and taxonomy files.

```
mapfile = "demodataR.csv"
# Import sample metadata
map <- read.csv(mapfile)

# Convert year into categorical factor
map$year <- as.factor(map$year)

# Convert this dataframe into phyloseq format
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
otu_table()   OTU Table:         [ 395 taxa and 12 samples ]
sample_data() Sample Data:       [ 12 samples by 3 sample variables ]
tax_table()   Taxonomy Table:    [ 395 taxa by 6 taxonomic ranks ]
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
mydata <- merge_phyloseq(moth_merge)
```

At this point you should have a Data object called `mydata` in the **Environment** panel. This object should be a Formal *class phyloseq*.

# Basic info and plots

**Calculate the number of reads per sample**

```
sample_sums(mydata)
```

Output should look like (but contain more values):
```
    2017a1PCRneg     2017a2PCRneg      2017aDNAneg     2017b1PCRneg 
           13712             4205             5113            53672
```
**Plot it**
```
plot_bar(mydata)
```
**Before you move on**, show your plot to your neighbor or the instructor. What does each bar represent? How is each bar divided? 
Does your plot look like this?
![Raw Bar Plot](demodata.rawbarplot.png)


# Bar plots

Of course, you know that your samples had different numbers of reads after the QC, so you should convert your dataset to relative abundances. Use this command and create a new dataset
```
relmydata = transform_sample_counts(mydata,function(x) 100 * x / sum(x))
```
Here you’ll divide all the OTU counts by the total sample counts and then multiple by 100. Now your bars will sum to 100% and represent the relative abundance within a sample. 

You can use Phylseq's built in function to color your bar plot
```
plot_bar(relmydata,fill="Class")
```
The above command is plotting all the OTUs colored by Class. This can get pretty confusing pretty quickly. You can use the below code to pull together some of the OTUs by whatever taxonomic level you're intereseted in.

```
relmydata_phylum <- relmydata %>%
  tax_glom(taxrank = "Phylum") %>%                     # group at Phylum level
  psmelt() %>%                                         # Melt to long format
  filter(Abundance > 1) %>%                         # Filter out low abundance taxa
  arrange(Phylum)                                   # Sort data frame alphabetically by Phylum
```
Here, we are combining at the **Phylum** level and filtering out any phylum that is represented less than 1%

Next we plot the results
```
phylum_colors <- diverge_hcl(length(unique(relmydata_phylum$Phylum)))
ggplot(relmydata_phylum, aes(x = SampleID, y = Abundance, fill = Phylum)) + 
  geom_bar(stat = "identity") +
  scale_fill_manual(values = phylum_colors) +
  # Remove x axis title
  theme(axis.title.x = element_blank()) + 
  ylab("Relative Abundance (Phylum > 1%) \n") +
  theme(axis.text.x=element_text(angle=90,hjust=1)) +
  ggtitle("Composition, Phylum")
```
Does your plot look like this?
![Phylum Bar Plot](demodata.phylumbarplot.png)

**Keep digging deeper into the data**

Now, let's look at class and family level. Below is code to combine the OTUs at each of those levels.
```
relmydata_class <- relmydata %>%
  tax_glom(taxrank = "Class") %>% 
  psmelt() %>% 
  filter(Abundance > 1) %>% 
  arrange(Class)

relmydata_family <- relmydata %>%
  tax_glom(taxrank = "Family") %>% 
  psmelt() %>% 
  filter(Abundance > 1) %>% 
  arrange(Family)            
```
Here is another exmaples to create some pretty color pallets for your categories.  There are lots more!
```
class_colors <- rainbow_hcl(length(unique(relmydata_class$Class)))
```
Finally, you can plot each level separately. First **Class**
```
ggplot(relmydata_class, aes(x = SampleID, y = Abundance, fill = Class)) + 
  geom_bar(stat = "identity") +
  scale_fill_manual(values = class_colors) +
  # Remove x axis title
  theme(axis.title.x = element_blank()) + 
  ylab("Relative Abundance (Class > 10%) \n") +
  theme(axis.text.x=element_text(angle=90,hjust=1)) +
  ggtitle("Composition, Class") 
```
and now by **Family**
In this case, we've defined the color directly in the `ggplot` command (see the `scale_fill_discrete_qualitative`). This is another alternative method.
```
ggplot(relmydata_family, aes(x = SampleID, y = Abundance, fill = Family)) + 
  geom_bar(stat = "identity") +
  scale_fill_discrete_qualitative(palette = "Dark 3") +
  # Remove x axis title
  theme(axis.title.x = element_blank()) + 
  ylab("Relative Abundance (Family > 1%) \n") +
  theme(axis.text.x=element_text(angle=90,hjust=1)) +
  ggtitle("Composition, Family")
```

**Before you move on**, show your plot to your neighbor or the instructor. What does each bar represent? How is each bar divided? 

# Ordinations
One of the best exploratory analyses for amplicon data is an ordinations to visualize the beta diversity. Phyloseq can compute these in two simple steps. You calculate the distances between all your points, then you plot that data.

Caculation
```
# Calculate
mydata_nmds_bray <- ordinate(
  physeq = mydata, 
  method = "NMDS",
  distance = "bray"
)
````
Plotting
````
# Plot
plot_ordination(
  physeq = mydata,
  ordination = mydata_nmds_bray,
  title = "NMDS of mydata (bray)",
  color = "year",
  shape = "type"
) + 
  geom_point(aes(color = year), alpha = 0.7, size = 4)
```

# Alpha Diversity

Estimating alpha diversity of microbial communities is problematic no matter what you do. My best stab at it is to subsample the libraries with replacement to estimate the species abundance of the real population while standardizing sampling effort.

```
min_lib <- min(sample_sums(mydata))
```
We will subsample to 4205, the minimum number of reads. We will repeat this 100 times and average the diversity estimates from each trial.

**Initialize matrices to store richness and evenness estimates**
```
nsamp = nsamples(mydata)
trials = 100

richness <- matrix(nrow = nsamp, ncol = trials)
row.names(richness) <- sample_names(mydata)

evenness <- matrix(nrow = nsamp, ncol = trials)
row.names(evenness) <- sample_names(mydata)

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
*Phyloseq can also do the hard work for you with a one line code*
```
plot_richness(mydata, measures = "InvSimpson")
```
There are many metrics built into Phyloseq
```
# Shannon Diveristy
plot_richness(mydata,measures = "Shannon")
# Chao1 measure
plot_richness(mydata, measures="Chao1")
```
