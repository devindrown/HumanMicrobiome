#Load libraries
library(ggplot2)
library(vegan)
library(dplyr)
library(scales)
library(grid)
library(reshape2)
library(phyloseq)
library(colorspace)
library(ape)

# Set working directory
setwd("MYHOMEDIRECTORY/demodataR")

# Set plotting theme
theme_set(theme_bw())

# Assign variables for imported data
sharedfile = "demodataR.shared"
taxfile = "demodataR.taxonomy"

# Import mothur data
mothur_data <- import_mothur(mothur_shared_file = sharedfile,
                             mothur_constaxonomy_file = taxfile)

mapfile = "demodataR.csv"
# Import sample metadata
map <- read.csv(mapfile)

# Convert this dataframe into phyloseq forma
map <- sample_data(map)

# Assign rownames to be Sample ID's
rownames(map) <- map$SampleID

# Merge mothurdata object with sample metadata
moth_merge <- merge_phyloseq(mothur_data, map)

colnames(tax_table(moth_merge)) <- c("Kingdom", "Phylum", "Class", 
                                     "Order", "Family", "Genus")
# Create a tree
random_tree = rtree(ntaxa(moth_merge), rooted=TRUE, tip.label=taxa_names(moth_merge))
mydata <- merge_phyloseq(moth_merge, random_tree)

# Alpha diversity
min_lib <- min(sample_sums(mydata))

nsamp = nsamples(mydata)
trials = 100

richness <- matrix(nrow = nsamp, ncol = trials)
row.names(richness) <- sample_names(mydata)

evenness <- matrix(nrow = nsamp, ncol = trials)
row.names(evenness) <- sample_names(mydata)

set.seed(42)

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

alpha <- rbind(rich_stats, even_stats)
s <- data.frame(sample_data(mydata))
alphadiv <- merge(alpha, s, by = "SampleID") 

ggplot(alphadiv, aes(x = SampleID, y = mean)) +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), colour="black", width=.1) +
  geom_point(size = 2) +
  facet_wrap(~measure, ncol = 1, scales = "free")

# BAR PLOTS
relmydata = transform_sample_counts(mydata,function(x) 100 * x / sum(x))

plot_bar(relmydata,fill="Class")

relmydata_phylum <- relmydata %>%
  tax_glom(taxrank = "Phylum") %>%                     # agglomerate at Phylum level
  psmelt() %>%                                         # Melt to long format
  filter(Abundance > 1) %>%                         # Filter out low abundance taxa
  arrange(Phylum)                                   # Sort data frame alphabetically by Phylum

phylum_colors <- diverge_hcl(length(unique(relmydata_phylum$Phylum)))

ggplot(relmydata_phylum, aes(x = SampleID, y = Abundance, fill = Phylum)) + 
  geom_bar(stat = "identity") +
  scale_fill_manual(values = phylum_colors) +
  # Remove x axis title
  theme(axis.title.x = element_blank()) + 
  ylab("Relative Abundance (Phylum > 10%) \n") +
  theme(axis.text.x=element_text(angle=90,hjust=1)) +
  ggtitle("Composition, Phylum") 

relmydata_class <- relmydata %>%
  tax_glom(taxrank = "Class") %>%                     # agglomerate at class level
  psmelt() %>%                                         # Melt to long format
  filter(Abundance > 1) %>%                         # Filter out low abundance taxa
  arrange(Class)                                      # Sort data frame alphabetically by class

relmydata_family <- relmydata %>%
  tax_glom(taxrank = "Family") %>%                     # agglomerate at family level
  psmelt() %>%                                         # Melt to long format
  filter(Abundance > 1) %>%                         # Filter out low abundance taxa
  arrange(Family)            

class_colors <- diverge_hcl(length(unique(relmydata_class$Class)))
family_colors <- rainbow_hcl(length(unique(relmydata_family$Family)))

ggplot(relmydata_class, aes(x = SampleID, y = Abundance, fill = Class)) + 
  geom_bar(stat = "identity") +
  scale_fill_manual(values = class_colors) +
  # Remove x axis title
  theme(axis.title.x = element_blank()) + 
  ylab("Relative Abundance (Class > 10%) \n") +
  theme(axis.text.x=element_text(angle=90,hjust=1)) +
  ggtitle("Composition, Class") 

ggplot(relmydata_family, aes(x = SampleID, y = Abundance, fill = Family)) + 
  geom_bar(stat = "identity") +
  scale_fill_manual(values = family_colors) +
  # Remove x axis title
  theme(axis.title.x = element_blank()) + 
  ylab("Relative Abundance (Family > 1%) \n") +
  theme(axis.text.x=element_text(angle=90,hjust=1)) +
  ggtitle("Composition, Phylum") 

# Plot
mydata_pcoa_bray <- ordinate(
  physeq = mydata, 
  method = "PCoA",
#  weight=TRUE,
  distance = "bray"
)

mydata_pcoa_Wunifrac = ordinate(mydata,"PCoA","unifrac",weighted=TRUE)
mydata_pcoa_unifrac = ordinate(mydata,"PCoA","unifrac",weighted=FALSE)

mydata_NMDS_bray <- ordinate(
  physeq = mydata, 
  method = "NMDS", 
  distance = "bray"
)
mydata_NMDS_Wunifrac = ordinate(mydata,"NMDS","unifrac",weighted=TRUE)
mydata_NMDS_unifrac = ordinate(mydata,"NMDS","unifrac",weighted=FALSE)

plot_ordination(
  physeq = mydata,
  ordination = mydata_pcoa_bray,
  title = "PCoA of mydata (bray)",
  color = "SampleID"
) + 
  geom_point(aes(color = SampleID), alpha = 0.7, size = 4)

plot_ordination(mydata, mydata_pcoa_Wunifrac,
                title = "PCoA of mydata (Weighted Unifrac)",
                color = "SampleID") + 
  geom_point(aes(color = SampleID), alpha = 0.7, size = 4)

plot_ordination(mydata, mydata_pcoa_unifrac,
                title = "PCoA of mydata (Unweighted Unifrac)",
                color = "SampleID") + 
  geom_point(aes(color = SampleID), alpha = 0.7, size = 4)

plot_ordination(
  physeq = mydata,
  ordination = mydata_NMDS_bray,              
  title = "NMDS of mydata (bray)",
  color = "SampleID"
) + 
  geom_point(aes(color = SampleID), alpha = 0.7, size = 4)

 plot_ordination(mydata, mydata_NMDS_Wunifrac,
                title = "NMDS of mydata (Weighted Unifrac)",
                color = "SampleID") + 
  geom_point(aes(color = SampleID), alpha = 0.7, size = 4)

plot_ordination(mydata, mydata_NMDS_unifrac,
                title = "NMDS of mydata (Unweighted Unifrac)",
                color = "SampleID") + 
  geom_point(aes(color = SampleID), alpha = 0.7, size = 4)
