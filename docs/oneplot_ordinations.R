# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

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

pb <- plot_ordination(
  physeq = mydata,
  ordination = mydata_pcoa_bray,
  title = "PCoA bray",
  color = "SampleID"
) + 
  geom_point(aes(color = SampleID), alpha = 0.7, size = 4) +
  theme(legend.position="none")

pwu <- plot_ordination(mydata, mydata_pcoa_Wunifrac,
                       title = "PCoA Wt Unifrac)",
                       color = "SampleID") + 
  geom_point(aes(color = SampleID), alpha = 0.7, size = 4) +
  theme(legend.position="none")

puu <- plot_ordination(mydata, mydata_pcoa_unifrac,
                       title = "PCoA UnWt Unifrac",
                       color = "SampleID") + 
  geom_point(aes(color = SampleID), alpha = 0.7, size = 4) +
  theme(legend.position="none")

nb <- plot_ordination(
  physeq = mydata,
  ordination = mydata_NMDS_bray,              
  title = "NMDS bray",
  color = "SampleID"
) + 
  geom_point(aes(color = SampleID), alpha = 0.7, size = 4) +
  theme(legend.position="none")

nwu <- plot_ordination(mydata, mydata_NMDS_Wunifrac,
                       title = "NMDS Wt Unifrac",
                       color = "SampleID") + 
  geom_point(aes(color = SampleID), alpha = 0.7, size = 4) +
  theme(legend.position="none")

nuu <- plot_ordination(mydata, mydata_NMDS_unifrac,
                       title = "NMDS UnWt Unifrac",
                       color = "SampleID") + 
  geom_point(aes(color = SampleID), alpha = 0.7, size = 4) +
  theme(legend.position="none")

multiplot(pb, nb, pwu, nwu, puu, nuu, cols=3)
