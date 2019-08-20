#######################################
##           javier cordoba          ##
##             script for            ##
##  plotting GC content distribution ##
#######################################

###################
## load packages ##
###################
library(ggplot2)
 
#################################
## load Multiple plot function ##
#################################

# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.

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
#########################
## path specifications ##
#########################
dir="C:/Users/LANCIOS-PC/Documents/Oficina/Labo/Paper_Transcriptome/Data/new_assembly3/"

exp = c("EML","DRK","KLY","UNK","YSD","KLY.cl")

myList = list()

dev.off()
for (i in 1:length(exp)) {
  myList[[i]] <- read.csv(file=paste0(dir,exp[i],".txt"), sep="\t")
  data <- myList[[i]]
  data <- data[data$Length > 200,]
  h = mean(data$GC)
  data <- as.data.frame(myList[[i]])
  
  p1 <- ggplot(data=data, aes(myList[[i]]$GC, fill= ..count..)) + 
    geom_histogram(colour ="black", bins = 40) + 
    geom_vline(xintercept = h, col="red") + 
    labs(y="Frequency", x= '') + 
    theme_classic() +
    theme(legend.position = "none", text = element_text(size=15), axis.text.y = element_blank()) +
    scale_x_continuous(limits=c(0, 1), breaks= seq(from=0.2, to=0.8, by=0.2))
  
  p2 <- ggplot(data=data, aes(myList[[i]]$GC, myList[[i]]$Length)) + 
    geom_hex(colour="black", bins = 40) + 
    labs(y='', x= '') + 
    theme_classic() + 
    theme(legend.position = "none", text = element_text(size=15), axis.text.y = element_blank()) +
    scale_x_continuous(limits=c(0, 1), breaks= seq(from=0.2, to=0.8, by=0.2))
  
  # export to file
  png(file= paste0("C:/Users/LANCIOS-PC/Documents/Oficina/Labo/Paper_Transcriptome/Grph/new_assembly3/pval/",
                   "GC_content_", exp[i], ".png"),
      width=15, height=7, units="cm", res=300, pointsize=1)
  
  multiplot(p1, p2, cols=2)
  grid.text("GC content", x = 0.45, y=1, just = "top", vp = viewport(layout.pos.row = 1, layout.pos.col = 1:2), gp = gpar(fontsize = 10, fontface = "bold", fontfamily="Times"))
  grid.text(paste0(exp[i]," = ",round(h, digits = 2)), x= 0.45, y= 0.9, just = "bottom", vp = viewport(layout.pos.row = 1, layout.pos.col = 1:2), gp = gpar(fontsize = 8, fontface = "bold", fontfamily="Times"))
  grid.text("A", y=0.95, x= 0.08, just = "top", vp = viewport(layout.pos.row = 1, layout.pos.col = 1:2), gp = gpar(fontsize = 8, fontface = "bold", fontfamily="Times"))
  grid.text("B", y=0.95, x= 0.58, just = "top", vp = viewport(layout.pos.row = 1, layout.pos.col = 1:2), gp = gpar(fontsize = 8, fontface = "bold", fontfamily="Times"))
 
  dev.off()
}
