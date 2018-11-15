# Code Snacks

Below are some bits of code to help you work on your Microbiome Research Projects. These are not complete instructions, but rather some pointers on where you might go to help you explore the metadata.


I assume that you have a datasets called myhouse and mysite. In case you forgot, you can `subset` your data by typing
```
mysite <- subset_samples(mb, Site=="SiteC")
```

For simplicity, let's put the data you'd like to visualize into a common container, `mydata`
```
mydata <-mysite
```

# Filtering samples

We have 192 samples in our dataset. Not all of these are appropriate for all types of analysis. Sometimes we are also missing metadata.

## How to remove NA from your sample sets

Missing data in has been coded as `NA`, so it is easy to filter out.
```
mydataQC <- subset_samples(mydata, Undergraduate != "NA")
```

Typing in `mydata` and `mydataQC` shows you that 1 sample was removed
```
> mydata
sample_data() Sample Data:       [ 20 samples by 17 sample variables ]
> mydataQC
sample_data() Sample Data:       [ 19 samples by 17 sample variables ]
```

You can see the results of filtering out this unidentified sample. Compare the following plots of divisity
```
plot_richness(mydata, x = "Undergraduate", measures="Chao1")
plot_richness(mydataQC, x = "Undergraduate", measures="Chao1")
```

## How to remove negative controls from your sample sets

Let's look at another example where you might want to just focus on the samples from swabs rather than the controls. I made a column in the metadata called *Type* and for each negative control, that sampled was labeled as *QC*
```
myhouseQC <- subset_samples(myhouse, Type != "QC")
```
Here we are keeping all samples that **do not** have a *Type* of *QC*. 

You can see the results of filtering by comparing the following plots of divisity
```
plot_richness(myhouse, x = "Site", measures="Chao1")
plot_richness(myhouseQC, x = "Site", measures="Chao1")
```
You could use thise method on any column in the metadata.

# Visualizing Two variables

Let's go back to your Site data that you have stored in `mydataQC`

To generate a PCoA plot using Undergraduate status to color the points you would first generate teh distance matrix
```
# Calculate the distance matrix first
mydataQC_pcoa_bray <- ordinate(
  physeq = mydataQC, 
  method = "PCoA",
  #  weight=TRUE,
  distance = "bray"
)
```
The you can plot the data. Here you have to specific *Undergraduate* twice in the plot code
```
plot_ordination(
  physeq = mydataQC,
  ordination = mydataQC_pcoa_bray,
  title = "PCoA of mydata (bray)",
  color = "Undergraduate"
) + 
  geom_point(aes(color = Undergraduate), alpha = 0.7, size = 4)
```

Perhaps you want to explore both Undergraduate status and Dry Cabin living at the same time. Before you do this though, you shoudl `subset` your data again to make sure you don't have rows without data in both categories

```
mydataQC <- subset_samples(mydata, Undergraduate != "NA" & DryCabin != "NA")
```

Then you can generate the distance matrix again as above `mydataQC_pcoa_bray <- ordinate...`
Then run the plot command. You add a single line to the plot specificying the *symbol*. 

```
plot_ordination(
  physeq = mydataQC,
  ordination = mydataQC_pcoa_bray,
  title = "PCoA of mydata (bray)",
  shape = "DryCabin",  color = "Undergraduate"
) + 
  geom_point(aes(color = Undergraduate), alpha = 0.7, size = 4)
```


