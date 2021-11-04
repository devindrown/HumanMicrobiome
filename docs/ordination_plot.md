# Ordination Plots, visualize beta diversity, community composition differences

In this example, we'll start with a data set containing a mix of houses and sites. We can copy the site into a generic container, **`mydata`**.
```
mydata <- mycomplexdata
```

Visualizating ordinations takes two separate steps. First, we have to calculate the distances among all of our samples. Second, we plot those distances.

### Calculate distances
As before we calculate distances  using the `ordinate` function. In this example, we'll generate a `PCoA` graph using bray-curtis distances, `bray`.
```
mydata_ord <- ordinate(
  physeq = mydata, 
  method = "PCoA",
  distance = "bray"
)
```
### Plot distances
Next, we using the `plot_ordination` function to visualizae those distances.
```
plot_ordination(
  physeq = mydata,
  ordination = mydata_ord
) 
```
In this example, my plot looks like the one below. All of the samples are colored the same.
![Raw Ord Plot](myord.demo.1.png)

In most cases, we are using this kind of plot to visualize the differences among at least two different groups. We can color each of our points by membership in that group using the `color = ` parameter. In this examples, we'll color by sequencing plate (a grouping with a null expectation of differences).
```
plot_ordination(
  physeq = mydata,
  ordination = mydata_ord,
  color = "Plate"
) 
```
![Raw Ord Plot](myord.demo.2.png)

Everything is still quite small. If we had a lot of points, it would be quite hard to see overlapping samples. We can fix both of these by increaseing the size (`size = `) and making the points transparent (`alpha = `) within the `geom_point` function of `ggplot`. We will modify the `phyloseq` object with the additional code.
```
plot_ordination(
  physeq = mydata,
  ordination = mydata_ord,
  color = "Plate"
) +
  geom_point(aes(color = Plate),size=5, alpha = 0.7) +
  labs(color='Plate') 
```
![Raw Ord Plot](myord.demo.3.png)

There may be times when you want to label individual points. This can be done with the `geom_text` function of `ggplot`. In the below example, we'll label each point with the `House` and move the labels horizontally (`hjust`) and vertically (`vjust`) so they don't overlap with the point.
```
plot_ordination(
  physeq = mydata,
  ordination = mydata_ord,
  color = "Plate"
) +
  geom_point(aes(color = Plate),size=5, alpha = 0.7) +
  labs(color='Plate')  +
  geom_text(aes(label = House), size=3, hjust=1.2, vjust=1.2)
```
![Raw Ord Plot](myord.demo.4.png)

### Adding an elipse
To help visualize the clustering of the different groups in our plot, we might want to provide an elipse around each group. We can do this pretty simply withte `stat_ellipse` function. We will assume a multivariate distribution with the parameter `type = norm`. In the below example, we will add the elipse first, and then plot the point on top. 
```
plot_ordination(
  physeq = mydata,
  ordination = mydata_ord,
  color = "Plate"
) +
  stat_ellipse(geom="polygon",type="norm", alpha=0.1, show.legend=F, aes(fill = Plate)) +
  geom_point(aes(color = Plate),size=5, alpha = 0.7) +
  labs(color='Plate')
```
![Raw Ord Plot](myord.demo.5.png)

The elipse is pretty wide. Since this is just a visualize aid and we'll do some stats in a moment, let's use the `level =` parameter to shrink that elipse down.
```
plot_ordination(
  physeq = mydata,
  ordination = mydata_ord,
  color = "Plate"
) +
  stat_ellipse(geom="polygon",type="norm", level=0.75, alpha=0.1, show.legend=F, aes(fill = Plate)) +
  geom_point(aes(color = Plate),size=5, alpha = 0.7) +
  labs(color='Plate')
```
![Raw Ord Plot](myord.demo.6.png)

### Final clean up
We can modify this plot with a few more tweaks. If we fix the coordinate system, then the x and y axes will be on the same scale. You can change the color palette to using `ggsci` package. You’ll need to install the ggsci package and load the library before this will work. Here we use the two functions scale_color_lancet() and scale_fill_lancet() to pick colors like the journal Lancet. We're going to get rid of the gray background with `theme_minimal`. Finally, let's standardize the text size with some additional `theme` adjustments.
```
plot_ordination(
  physeq = mydata,
  ordination = mydata_ord,
  color = "Plate"
) +
  coord_fixed() +
  stat_ellipse(geom="polygon",type="norm", level=0.75, alpha=0.1, show.legend=F, aes(fill = Plate)) +
  geom_point(aes(color = Plate),size=5, alpha = 0.7) +
  labs(color='Plate') +
  scale_color_lancet() + scale_fill_lancet() +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 12, face = "bold"),
    axis.title = element_text(size = 14,face = "bold"),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14),
  )
```
![Raw Ord Plot](myord.demo.7.png)





We will use a series of commands to combine OTU/ASVs by the taxonomic rank. 
First let's create a variable with the taxonomic rank at which you want to combine OTU/ASVs. We'll call this variable **`myTaxLevel`**. In this example, we'll be combining at the Phylum level.
```
myTaxLevel <- "Phylum"
```

You should also consider filtering out low abudance groups. We'll call this variable **`myfilter`**. In this example, we'll filter out any group that is present at less than 5% (`< 5%`).
```
myFilter <- 5
```

Now you can create a string variable with all of the information. You'll use this to label, **`myYaxis`** , for your vertical axis.
```
myYaxis <- paste("Relative Abundance (", myTaxLevel, " > ", myFilter, "%) \n")
```
Now you are ready to start reorganizing your data. The code below should work for any taxonomic level and any degree of filtering. You should be able to change the values of `myTaxLevel` and `myFilter` and then run the commands below to generate any number of different plots.


## Covert to relative abundance and combine

Start by converting the raw count data into relative abundance data as a percent (0 - 100 range). This creates a new data set called **`relmydata`**.
```
relmydata = transform_sample_counts(mydata,function(x) 100 * x / sum(x))
```
You've already set up your variables and filters so you can combine your OTU/ASVs. Your combined data will now be in the container called **`relmydata_grouped`**.
```
relmydata_grouped <- relmydata %>%
  tax_glom(taxrank = myTaxLevel) %>%        # group at your Taxonomic level
  psmelt() %>%                              # Melt to long format
  filter(Abundance > myFilter) %>%          # Filter out low abundance taxa
  arrange(myTaxLevel)                       # Sort data frame alphabetically by your Taxonomic level
```

We can clean up this data and remove unclassified OTU/ASVs labeled as `Bacteria_unclassified`. Your cleaned data will now be called **`relmydata_grouped_clean`**.
```
relmydata_grouped_clean <- subset(relmydata_grouped, relmydata_grouped[[myTaxLevel]] != "Bacteria_unclassified")
```

Next, create a color palette based on the number of groups within your taxonomic level. In this example, we'll use the Viridis palette to pull colors. The colors will be stored in a variable called **`mycolors`**.
```
mycolors <- sequential_hcl(length(unique(relmydata_grouped_clean[[myTaxLevel]])), palette = "viridis")
```

Finally, we can put this all together and use `ggplot` to graph the data.
```
ggplot(relmydata_grouped_clean, aes_string( x = "House", y = "Abundance", fill = myTaxLevel)) + 
  geom_bar(stat = "identity") +
  scale_fill_manual(values = mycolors) +
  ylab(myYaxis) +
  ggtitle("Community Composition")
```
![Raw Bar Plot](mysite.demo.1.png)


You can improve this graph. By adding classic theme (`theme_classic()`), you can remove the gray background on the graph and add some lines for the axes.
```
ggplot(relmydata_grouped_clean, aes_string( x = "House", y = "Abundance", fill = myTaxLevel)) + 
  geom_bar(stat = "identity") +
  scale_fill_manual(values = mycolors) +
  ylab(myYaxis) +
  ggtitle("Community Composition") +
  theme_classic()
```
![Raw Bar Plot](mysite.demo.2.png)


The labels on the horizontal axis are overlapping. You can rotate the text using the `axis.text.x` call in `theme()`.
```
ggplot(relmydata_grouped_clean, aes_string( x = "House", y = "Abundance", fill = myTaxLevel)) + 
  geom_bar(stat = "identity") +
  scale_fill_manual(values = mycolors) +
  ylab(myYaxis) +
  ggtitle("Community Composition") +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
  )
```
![Raw Bar Plot](mysite.demo.3.png)


The font size is still quite small. We can increase and standardize the size in the `theme`.
```
ggplot(relmydata_grouped_clean, aes_string( x = "House", y = "Abundance", fill = myTaxLevel)) + 
  geom_bar(stat = "identity") +
  scale_fill_manual(values = mycolors) +
  ylab(myYaxis) +
  ggtitle("Community Composition") +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
    axis.text = element_text(size = 12, face = "bold"),
    axis.title = element_text(size = 14,face = "bold"),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14),
  )
```
![Raw Bar Plot](mysite.demo.4.png)


If you want, you can add dashed grid lines to help see differences better with the `panel.grid` call in `theme()`.
```
ggplot(relmydata_grouped_clean, aes_string( x = "House", y = "Abundance", fill = myTaxLevel)) + 
  geom_bar(stat = "identity") +
  scale_fill_manual(values = mycolors) +
  ylab(myYaxis) +
  ggtitle("Community Composition") +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
    axis.text = element_text(size = 12, face = "bold"),
    axis.title = element_text(size = 14,face = "bold"),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14),
    axis.ticks = element_blank(),
    panel.grid = element_line(color = "#b4aea9"),
    panel.grid.major.y = element_line(linetype = "dashed"),
  )
```
![Raw Bar Plot](mysite.demo.5.png)

You change the color palette to using `ggsci` package. You'll need to install the `ggsci` package and load the library before this will work. Here we use the two functions `scale_color_lancet()` and `scale_fill_lancet()` to pick colors like the journal *Lancet*.
```
ggplot(relmydata_grouped_clean, aes_string( x = "House", y = "Abundance", fill = myTaxLevel)) + 
  geom_bar(stat = "identity") +
  # scale_fill_manual(values = mycolors) +
  scale_color_lancet() + scale_fill_lancet() +
  ylab(myYaxis) +
  ggtitle("Community Composition") +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
    axis.text = element_text(size = 12, face = "bold"),
    axis.title = element_text(size = 14,face = "bold"),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14),
    axis.ticks = element_blank(),
    panel.grid = element_line(color = "#b4aea9"),
    panel.grid.major.y = element_line(linetype = "dashed"),
  )
```
![Raw Bar Plot](mysite.demo.6.png)

If you'd like to be a little silly, you can change the color palette to Simpsons theme.
```
ggplot(relmydata_grouped_clean, aes_string( x = "House", y = "Abundance", fill = myTaxLevel)) + 
  geom_bar(stat = "identity") +
  # scale_fill_manual(values = mycolors) +
  # scale_color_lancet() + scale_fill_lancet() +
  scale_color_simpsons() + scale_fill_simpsons() +
  ylab(myYaxis) +
  ggtitle("Community Composition") +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
    axis.text = element_text(size = 12, face = "bold"),
    axis.title = element_text(size = 14,face = "bold"),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14),
    axis.ticks = element_blank(),
    panel.grid = element_line(color = "#b4aea9"),
    panel.grid.major.y = element_line(linetype = "dashed"),
  )
```
![Raw Bar Plot](mysite.demo.7.png)


