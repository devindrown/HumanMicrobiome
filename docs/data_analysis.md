# Data Analysis Exploring the Class Microbiome Dataset
In this part of the lab, you will apply the analysis skills you learned last week to a new, larger dataset. The goal is to explore the microbial diversity within different houses and sites from our class project.

# Step 1: Setting Up Your RStudio Environment
First, let's get your R session ready to work with the class data.

1. Create a New RStudio Project: `File>New Project`. If RStudio asks you to “Save Current Workspace”, you should select “Don’t Save”. Select Existing Directory. On the next window, set the Project working directory to: `~/BIOL491_2024`
2. Load the Class Data: In the RStudio Console (the bottom-left panel), type the following command and press Enter. This script will load several pre-processed datasets into your environment.
```
source('treasurechest/LoadData.R')
```
3. Create a New R Script: Go to `File > New File > R Script`. This will open a blank script in the editor panel. Immediately save this file (`File > Save`). Name it something informative, like `house_analysis.R`. Remember to include the .R extension! All the code you write for today's lab should go into this script.

## Step 2: Understanding the Loaded Datasets
The `LoadData.R` script provided you with three phyloseq objects. Each is a slightly different version of the class dataset:

* `mb`: The complete, raw dataset containing all 187 samples, including negative controls.
* `mbQC`: A quality-controlled version. Negative controls have been removed, and low-abundance ASVs (Amplicon Sequence Variants) have been filtered out, leaving 120 samples. You will use this for most of your analysis.
* `mb_dirty`: This version has negative controls removed but still includes some known contaminating ASVs.

# Step 3: Create Your Personal Datasets

The full dataset is too large to analyze all at once. Your first task is to create smaller, manageable subsets to work with

1. Subset by House:
* First, see which houses are available by running this command in your console: 
```levels(mb@sam_data$House)```
* Choose one house ID (e.g., ab8a) and use the subset_samples() function to create a new phyloseq object containing only the samples from that house.
```
# This creates a new object 'myhouse' containing only samples where the 'House' column is "ab8a".
# Replace "ab8a" with the house ID you chose.
myhouse <- subset_samples(mbQC, House == "ab8a")
```

2. Subset by Site:
* Next, choose a specific sample site to investigate (e.g., `SiteA`, `SiteB`). You can see all available site labels by inspecting the metadata: ```view(mbQC@sam_data)```.
* Use the same subsetting technique to create a dataset for your chosen site.
```
# This creates 'mysite' containing only samples where 'Site' is "SiteZ".
# Replace "SiteZ" with the site you chose.
mysite <- subset_samples(mbQC, Site == "SiteZ")
```

**Pro Tip**: The analysis script from last week expects your data to be in an object named `mydata`. To easily reuse that code, you can copy your new subset into mydata like this:
```mydata <- myhouse``` or ```mydata <- mysite```.

# Step 4: Analyze Your Subsets

Now you are ready to analyze the diversity of your `myhouse` and `mysite` datasets. Refer back to the script and instructions from the previous lab for detailed guidance on how to perform the following analyses.

1. **Calculate Read Counts**: Determine the sequencing depth for each sample in your new subset. You may refer back to the previous lab, [Phyloseq and R for analysis and visualization](phyloseq_analysis_visualization)
2. **Alpha Diversity**: Calculate and plot alpha diversity metrics (e.g., Richness, Inverse Simpson). Instructions here: [Alpha Diversity](alpha_diversity_plot)
3. **Community Composition**: Create bar plots showing the taxonomic composition at the Phylum, Class, and Family levels. Instructions here: [Bar Plots of diversity at different scales](community_composition_plot)
4. Beta Diversity (Ordination): Create an NMDS ordination plot to visualize how sample communities relate to each other. Note: This will be most informative for your `mysite` dataset, which compares multiple houses at the same site. Hint: Check out the `treasurechest` for code.

**Checkpoint: Please show your instructor the set of figures you've generated for your house dataset before you proceed with analyzing the site data.**


# What's next

You now have a complete workflow for exploring microbiome data and a script full of your own analysis code. Next week, we will build on this foundation to perform statistical testing and formally test hypotheses.
