# Running R

Install the latest version of R. You can find the install [HERE](https://cran.rstudio.com/) You want to install the **base version**

One of the best tools for using R is called R-Studio. You can download the Free RStudio Desktop verison [HERE](https://www.rstudio.com/products/rstudio/download/#download)



# Gettting your R install ready

Install R packages

```
install.packages("ggplot2")
install.packages("vegan")
install.packages("dplyr")
install.packages("scales")
install.packages("reshape2")
install.packages("colorspace")
install.packages("ape")
```

Install Phyloseq

```
source('http://bioconductor.org/biocLite.R')
biocLite('phyloseq')
```

You'll have to load all of these packages into your workspace before you can use them.

I've put all of the commands we used today in a single R script here (demodataR.R)
