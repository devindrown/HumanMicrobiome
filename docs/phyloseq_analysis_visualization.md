# Phyloseq and R for Analysis and Visualization

These instructions guide you through the analysis and visualization of microbiome data using the **phyloseq** package in R. You load data that was processed in QIIME2, create plots of taxonomic composition, and calculate alpha and beta diversity.

The lab is a notebook that runs in Google Colab:

* **[View the lab notebook on GitHub](https://github.com/devindrown/demodataR/blob/main/phyloseq_lab.ipynb)**: read the instructions and code, for review or reference.
* **[Open the lab notebook in Google Colab](https://colab.research.google.com/github/devindrown/demodataR/blob/main/phyloseq_lab.ipynb)**: run it yourself. Save a copy to your Drive first (**File → Save a copy in Drive**), and check that the runtime is set to **R**.

## What the lab covers

1. **Getting your workspace ready:** installing and loading phyloseq and the plotting packages.
2. **Importing your data:** loading a phyloseq object and identifying its parts (feature table, taxonomy, tree, sample metadata).
3. **Taxonomic composition plots:** reads per sample, relative abundance, and stacked bar plots at the Phylum, Class, Family, and Order levels.
4. **Beta diversity:** NMDS ordination with Bray-Curtis dissimilarity (and checking stress), PCoA with weighted and unweighted UniFrac, and PERMANOVA to test for differences between groups.
5. **Alpha diversity:** repeated rarefaction, richness and Inverse Simpson diversity, and comparing groups with boxplots.

The data in the lab are all **negative controls** from a 16S rRNA gene sequencing study, so along the way you also consider where the microbes in "blank" samples come from and what they mean for real samples.

## Using these skills with the class data

In Data Analysis I and II you apply the same tools to our class dataset, in Quarto notebooks in RStudio. See [Data Analysis I](data_analysis), [Alpha Diversity](alpha_diversity_plot), and [Bar Plots of diversity at different scales](community_composition_plot).
