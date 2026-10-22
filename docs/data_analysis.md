# Data Analysis I: Exploring the Class Microbiome Dataset
In this part of the lab, you will apply the analysis skills you learned last week to a new, larger dataset. The goal is to explore the microbial diversity within different houses and sites from our class project.

This year you will work in **Quarto notebooks** inside RStudio. A notebook mixes instructions and code cells: you run each cell and its output (tables, plots) appears right below it. This page is an overview. The notebooks have the step-by-step instructions.

# Step 1: Open the class project
1. Log in to RStudio.
2. Open the class project: `File > Open Project...`, go to the `2026.only` folder in your home directory, and open `2026.only.Rproj`. If RStudio asks you to "Save Current Workspace", select "Don't Save".
3. In the **Files** pane (bottom right) you should see:
   * `LoadData_walkthrough.qmd`: explains how the class data is loaded and prepared
   * `data_analysis_I.qmd`: this week's notebook
   * `data_analysis_II.qmd`: next week's notebook
   * `data/`: the class data, including the metadata spreadsheet `BIOL491_2026.metadata.tsv`
   * `treasurechest/`: helper scripts, including `LoadData.R`

**How to run a notebook:** click the green ▶ (*Run Current Chunk*) at the top right of each gray code cell, from top to bottom. If you restart R or something stops working, click the ▼ button (*Run All Chunks Above*) on the cell you are working on to catch up.

# Step 2: Understand the class data
Open `LoadData_walkthrough.qmd` and work through it first. Every analysis starts with this line, which runs the script `treasurechest/LoadData.R`:
```
source("treasurechest/LoadData.R")
```
The walkthrough contains the same code as that script, one step at a time, with an explanation of each step. It creates three phyloseq objects, each a slightly different version of the class dataset:

* `mb`: the complete dataset, including the negative controls.
* `mb_dirty`: negative controls and flagged samples (samples with too few reads to use) removed.
* `mbQC`: `mb_dirty` with rare ASVs (Amplicon Sequence Variants) removed. **You will use this for most of your analysis.**

The walkthrough also shows a histogram of sequencing depth (reads per sample) for the whole class, so you can see how the real samples compare with the negative controls.

# Step 3: Create your personal datasets
Open `data_analysis_I.qmd`. The full dataset is too large to analyze all at once, so your first task is to create two smaller subsets.

1. **Subset by house.** The notebook lists the house IDs. Choose your house and type its ID into this cell, in upper case exactly as listed (if the ID is not in the list, the cell stops and tells you so):
```
my_house <- "AB8A"   # <-- change to your house ID

myhouse <- subset_samples(mbQC, House == my_house)
```

2. **Subset by site.** The notebook also lists the sites. Choose one site to compare across all houses:
```
my_site <- "SiteZ"   # <-- change to the site you chose

mysite <- subset_samples(mbQC, Site == my_site)
```

The analysis code works on an object named `mydata`. To point it at one of your subsets, the notebook copies the subset into `mydata`: `mydata <- myhouse` or `mydata <- mysite`.

# Step 4: Analyze your subsets
The notebook walks you through each analysis, first for your house and then for your site:

1. **Sequencing depth**: a bar plot of the number of reads in each sample.
2. **Alpha diversity**: the Inverse Simpson index (and richness), built up one step at a time. See [Alpha Diversity](alpha_diversity_plot) for the same steps with explanations.
3. **Community composition**: stacked bar plots at the Phylum, Class, and Family levels, also built up step by step, ending with the [ggsci](https://nanx.me/ggsci/) color palettes. See [Bar Plots of diversity at different scales](community_composition_plot).
4. **Beta diversity (ordination)**: an NMDS ordination plot showing how communities relate to each other. This is most informative for your `mysite` dataset, which compares multiple houses at the same site.

**Checkpoint: Please show your instructor the set of figures you've generated for your house dataset before you proceed with analyzing the site data.**

To review the basics from last week's lab, see [Phyloseq and R for analysis and visualization](phyloseq_analysis_visualization).

# What's next
You now have a complete workflow for exploring microbiome data and a notebook full of your own analysis code. Next week, in `data_analysis_II.qmd`, we will build on this foundation to perform statistical testing and formally test hypotheses.
