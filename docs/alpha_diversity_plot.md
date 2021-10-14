# Alpha Diversity Plot

In this example, we'll start with a single house data set. We can copy the house into a generic container
```
mydata <- myhouse
```

Start with the basic plot from Phyloseq
```
plot_richness(mydata, measures = "InvSimpson")
```
![Raw Bar Plot](myhouse.demo.1.png)

Make the point larger with `geom_point` and setting `size = 10`
```
plot_richness(mydata, measures = "InvSimpson") +
  geom_point(size = 10)
```
![Raw Bar Plot](myhouse.demo.2.png)

Next, add some color and separate some of the sites in the house with `color = "Source"`
```
plot_richness(mydata, measures = "InvSimpson",
              color = "Source") +
  geom_point(size = 10)
```
![Raw Bar Plot](myhouse.demo.3.png)

Now, fix the X axis labels with `x = "Site"`
```
plot_richness(mydata, measures = "InvSimpson",
              color = "Source",
              x = "Site") +
  geom_point(size = 10)
```
![Raw Bar Plot](myhouse.demo.4.png)

Let's continue to fix the labels with the `labs()` function
```
plot_richness(mydata, measures = "InvSimpson",
              color = "Source",
              x = "Site") +
  geom_point(size = 10) +
  labs( x = "", y = "Inv Simpson")
```
![Raw Bar Plot](myhouse.demo.5.png)

We have some extra labels at the top, strip text. We can remove that redundancy.
```
plot_richness(mydata, measures = "InvSimpson",
              color = "Source",
              x = "Site") +
  geom_point(size = 10) +
  labs( x = "", y = "Inv Simpson") +
  theme(
    strip.text.x = element_blank(),
    strip.background = element_blank(),
  )
```
![Raw Bar Plot](myhouse.demo.6.png)

The font size is still quite small. We can increase and standardize the size in the `theme`.
```
plot_richness(mydata, measures = "InvSimpson",
              color = "Source",
              x = "Site") +
  geom_point(size = 10) +
  labs( x = "", y = "Inv Simpson") +
  theme(
    strip.text.x = element_blank(),
    strip.background = element_blank(),
    axis.text = element_text(size = 12, face = "bold"),
    axis.title = element_text(size = 14,face = "bold"),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14),
  )
```
![Raw Bar Plot](myhouse.demo.7.png)

Finally, the X axis labels are facing the wrong direction. We cna flip the angle in the `theme`.
```
plot_richness(mydata, measures = "InvSimpson",
              color = "Source",
              x = "Site") +
  geom_point(size = 10) +
  labs( x = "", y = "Inv Simpson") +
  theme(
    strip.text.x = element_blank(),
    strip.background = element_blank(),
    axis.text = element_text(size = 12, face = "bold"),
    axis.title = element_text(size = 14,face = "bold"),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)
  )
```
![Raw Bar Plot](myhouse.demo.8.png)
