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

# Filter out unwanted samples

We have 192 samples in our dataset. Not all of these are appropriate for all types of analysis. Sometimes we are also missing metadata.

## How to remove NA from your sample sets

Missing data in has been coded as `NA`, so it is easy to filter out.
```
mydataQC <- subset_samples(mydata, Undergraduate != "NA")
```

Compare the following plots of divisity
```
plot_richness(mydata, x = "Undergraduate", measures="Chao1")
plot_richness(mydataQC, x = "Undergraduate", measures="Chao1")
```

## How to remove negative controls from your sample sets
Let's look at another example where you 
