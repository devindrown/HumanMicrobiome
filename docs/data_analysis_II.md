# Data Analysis II: Statistics and Complex Datasets
Last week you explored your house and your site with plots. This week you test whether the patterns you saw are real: do groups of samples have different microbial communities? You will also work with larger subsets that combine several houses and sites.

As before, you work in a **Quarto notebook**, `data_analysis_II.qmd`. This page is an overview with explanations of the statistics. The notebook has the step-by-step instructions and all of the code.

# Step 1: Open the notebook and load the data
1. Open the class project (`2026.only.Rproj`) in RStudio, as in [Data Analysis I](data_analysis).
2. Open `data_analysis_II.qmd` from the **Files** pane.
3. Run the first cell. It loads the class data, exactly as last week:

   ```
   source("treasurechest/LoadData.R")
   ```

This gives you the same three phyloseq objects as last week: `mb` (everything, including negative controls), `mb_dirty` (controls and flagged samples removed) and `mbQC` (rare ASVs also removed). **Use `mbQC` for your analyses.**

**How to run a notebook:** click the green ▶ (*Run Current Chunk*) at the top right of each gray code cell, from top to bottom. If you restart R, click the ▼ button (*Run All Chunks Above*) on the cell you are working on to catch up.

# Step 2: Make your datasets (a review)
The notebook lists the house IDs and the sites. Type your house and your site into this cell. House IDs are in **upper case**, exactly as listed. If you mistype one, the cell stops and tells you so.
```
my_house <- "AB8A"    # <-- your house ID (upper case, as listed above)
my_site  <- "SiteZ"   # <-- the site you chose

myhouse <- subset_samples(mbQC, House == my_house)   # one home, all six sites
mysite  <- subset_samples(mbQC, Site == my_site)     # one site, all homes
```

The analysis code works on an object named `mydata`, so the notebook copies a subset into it: `mydata <- mysite`.

# Step 3: Introductory statistics
We start with your **site** dataset: one site, sampled in every house. Does the community at that site differ between groups of houses?

## Ordination (PCoA)
An ordination places samples so that similar communities sit close together. This week we use **PCoA** (principal coordinates analysis) with **Bray-Curtis** dissimilarity, which compares which taxa are present and how abundant they are.

The first example colors the samples by sequencing `Plate`. Plate was an arbitrary assignment during library preparation, so we do *not* expect it to matter. That makes it a useful check: if plate explained the differences between communities, it would point to a technical problem rather than biology.

## PERMANOVA: do the groups differ?
The notebook tests the plate grouping with `adonis2()` from the vegan package:
```
adonis2(phyloseq::distance(mydata, method = "bray") ~ Plate,
        data = data.frame(sample_data(mydata)), by = "terms")
```

What each part does:

* **`adonis2()`** runs a PERMANOVA (permutational multivariate analysis of variance). It tests whether groups of samples have different community **centers** (centroids).
* **`phyloseq::distance(mydata, method = "bray")`** calculates the Bray-Curtis dissimilarity between every pair of samples.
* **`~ Plate`** is the model: the grouping you are testing.
* **`data = data.frame(sample_data(mydata))`** is the table of metadata that contains `Plate`.
* **`set.seed(491)`**, on the line before, makes the result repeatable. PERMANOVA shuffles the samples at random (999 times) to build its p-value, so without a seed the p-value changes slightly each run.

Example output (SiteA, the toilet seat, 2026). Your numbers will differ.
```
         Df SumOfSqs     R2      F Pr(>F)
Plate     1   0.5187 0.0729 1.4154   0.12
Residual 18   6.5970 0.9271
Total    19   7.1157 1.0000
```

How to read it:

* **R²** is the share of the variation in community composition explained by the grouping. Here plate explains about 7%.
* **Pr(>F)** is the p-value. Here it is 0.12, above 0.05, so we **cannot reject** the null hypothesis that samples from the two plates share the same centroid. That is what we hoped for. (It does not *prove* the plates are the same.)

## Dispersion: are the groups equally spread out?
PERMANOVA can also give a significant result when groups differ in **spread** (dispersion) rather than in center. `betadisper()` measures how spread out each group is, and `permutest()` tests whether the spreads differ:
```
permutest(betadisper(phyloseq::distance(mydata, method = "bray"),
                     data.frame(sample_data(mydata))$Plate))
```

Example output:
```
          Df   Sum Sq  Mean Sq     F N.Perm Pr(>F)
Groups     1 0.030151 0.030151 2.336    999  0.146
Residuals 18 0.232322 0.012907
```

The p-value (0.146) is above 0.05, so the two groups are about equally spread out. This test matters most when the PERMANOVA is **significant**. If both are significant, part of the PERMANOVA result may come from a difference in spread, and you should say so when you report it.

## Test your own hypothesis
Write a simple, testable hypothesis for your site, for example: *"The microbial community at this site differs between households with and without pets."* Then choose the metadata column that matches it (`sample_variables(mydata)` lists them all) and run the same three steps: ordination, PERMANOVA and dispersion.

Two things to check before you trust a result:

* **Blank answers.** Some households left survey questions blank. A sample with a blank answer can't be placed in a group, so the notebook leaves it out of the tests (`ps_drop_incomplete()`) and prints how many samples were tested.
* **Group sizes.** The notebook prints how many samples are in each group. A test can't tell you much when a group has only one or two samples. This year, for example, only one household has cats. `Household` (the number of people) splits a site into several small groups, so to compare household sizes use **`HouseholdSize2`** (`1–2` vs. `3+` people) instead.

# Step 4: Complex subsets
To look at several sites **and** several houses at once, make two lists and keep only the samples that match both. Use IDs from the lists in the notebook. The notebook stops and names any site or house that isn't in the data.
```
mysitelist  <- c("SiteX", "SiteY", "SiteZ")       # <-- sites to include
myhouselist <- c("AB8A", "C3D1", "E5F2", "G7H3")  # <-- houses to include

mycomplexdata <- subset_samples(mbQC, Site %in% mysitelist & House %in% myhouselist)
```

`%in%` means "is one of": `Site %in% mysitelist` is TRUE for every sample whose site is in your list.

## Ordination with two variables
The PCoA uses **color** for site and **shape** for house, so you can see both groupings in one plot. Do samples cluster by site, by house, or both?

## PERMANOVA with two variables
A model can include more than one variable, written `~ A + B`:
```
adonis2(mydata_dist ~ House + Site, data = sampledf, by = "terms")
adonis2(mydata_dist ~ Site + House, data = sampledf, by = "terms")
```

With `by = "terms"`, the variables are tested **in order**: each one after accounting for the ones before it. That is why the notebook runs the model both ways.

Example output (three sites and four houses from 2026, 12 samples):
```
         Df SumOfSqs      R2      F Pr(>F)
House     3   1.2570 0.26920 1.4160  0.076 .
Site      2   1.6371 0.35058 2.7661  0.001 ***
Residual  6   1.7755 0.38023
Total    11   4.6696 1.00000
```

Here **Site** is significant (p = 0.001) and explains about 35% of the variation. House explains about 27% but is not significant (p = 0.076). In this example both orders give the same numbers, because every house has one sample from every site. If some samples are missing, the two orders can disagree. When they do, be more cautious about the variable whose result changes.

Then check the dispersion for the variable that came out significant, as in Step 3. If the dispersion test is not significant, you can be more confident that the PERMANOVA reflects a real difference in community composition.

## Alpha diversity and bar plots with two variables
* **Alpha diversity:** `plot_richness()` puts site on the x-axis and colors the points by house. How would you change it to compare houses instead? See [Alpha Diversity](alpha_diversity_plot) for the options.
* **Bar plots with facets:** `facet_grid(rows ~ columns)` splits a bar plot into panels. The notebook shows the same data three ways: one row per house (`House ~ .`), one row per site (`Site ~ .`), and a grid of every house and site (`House ~ Site`). The bar plot code is the same as last week; see [Bar Plots of diversity at different scales](community_composition_plot). To browse color palettes, run `hcl_palettes(plot = TRUE)`.

# Step 5: Advanced QC, redux
What does quality control do to your results? `mbQC` and `mb_dirty` contain the same samples. The only difference is that `mbQC` has the rare, noisy ASVs removed. The notebook makes the same house subset from each and compares the number of ASVs and the richness of each sample. Try the same comparison with a PERMANOVA or an ordination: does quality control change your conclusions?

# Planning your research project
Your project has two analyses: **(1)** your home across its six sites, and **(2)** one site across all homes. At the end of the notebook, list which plots and tests from Data Analysis I and II you will use for each.

To review the basics, see [Phyloseq and R for analysis and visualization](phyloseq_analysis_visualization).
