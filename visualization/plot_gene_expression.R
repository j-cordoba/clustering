gc()
#####################################
## plot "gene" expression patterns ##
#####################################
sek=seq(from=5, to=75, by=1) # do plot for clusters with number: from 5 to 75 (seq iterator)

dir="/path/to/my/directory"

rm(i)
# LOADING LOOP:
# load every kmeans objet created
for (i in sek){
  # read the correlation kmeans from computer
  km <- readRDS(file=paste0(dir,"Data/new_assembly3/pval/","km", i))
  
  cl = max(km$clustering) # define number of clusters
  
  # define graph limits, default -5 < 0 < 5
  # ylim.up = round(max(dataSCL))
  # ylim.lw = round(min(dataSCL))
  
  # background grey rectangles for sample division
  col_set2 = c("#FFFFFF", "#ABABAB")
  col_set2 <- sapply(col_set2, lighten)
  
  # colour blind friendly palette
  col_set = c("#0072B2", "#999999", "#E69F00", "#D55E00", "#F0E442", "#56B4E9", "#009E73", "#CC79A7")
  palette(col_set)
  
  rm(j)
  # LIST LOOP:
  for (j in 1:cl){
    # get gene identifiers from each cluster and pass to a list
    list <- rownames(dataSCL[km$clustering == j, ])
    write(list, ncolumns= 1, file= paste0(dir,"Grph/new_assembly3/cor/", "km", i, "list", j, ".ls"))
    }
  
  rm(j)
  # BOXPLOT LOOP: 
  # Export pdf of boxplot - expression patterns
  pdf(file= paste0(dir,"Grph/new_assembly3/cor/", "boxplot", i, ".pdf"))
  par(mfrow=c(3,3), mar=c(2,0.5,0,0.5))
  par(oma=c(2,3,3,1))
    
  for (j in 1:cl){
    # substraction of each cluster data group into a matrix
    M <- as.matrix(dataSCL[km$clustering == j,], ncol = ncol(dataSCL))
        
    # using boxplots
    boxplot(x=20,y=0, col= "white", ylim=c(-5,5), xaxt= "n", yaxt= "n", outpch = NA )
    
    # back ground preparation
    rect(xleft = c(0,3.5,6.5,7.5,8.5,9.5,10.5,13.5,16.5,19.5),
         xright =c(3.5,6.5,7.5,8.5,9.5,10.5,13.5,16.5,19.5,20.5),
         ybottom = par("usr")[2], ytop=par("usr")[3], col = (col_set2), border = "transparent") # background rectangles-sample separators
    
    # plot values
    boxplot(M, use.cols = TRUE, col= exp, ylim=c(-5,5), xaxt= "n", yaxt= "n", outpch = NA, add = TRUE) # overlap the second boxplot
    mtext(j, side= 3, line= -1.5, adj= 0.05, outer= FALSE, font= 2, cex=1) # assign number of cluster to plot
    axis(side= 2, at= seq(from= -5, to= 5, by= 2), cex.axis = 0.5, outer= T) # assign y axe scale
    axis(side= 1, at= c(1:23), labels= FALSE) # assign x axe ticks
    text(x= c(1:23), y= -6.5, tabla$Code ,xpd= TRUE, srt= 60, cex= 0.5) # assign names to the x ticks
    title(main= paste("Expression patterns with",cl,"clusters"), outer= T) # title
  }
  dev.off()
  
  rm(j)
  # MATPLOT LOOP:
  # Export pdf of matplot expression patterns
  pdf(file= paste0(dir,"Grph/new_assembly3/cor/", "matplot", i, ".pdf"))
  par(mfrow=c(3,3), mar=c(2,0.5,0,0.5))
  par(oma=c(2,3,3,1))
    
  for (j in 1:cl){
    # extract cluster number j from data
    M <- as.matrix(dataSCL[km$clustering == j,], ncol = ncol(dataSCL))
    
    # extract centroid from cluster number j
    centroid = rownames(dataSCL)[km$medoids[j]]
    centroid <- as.vector(dataSCL[centroid,])
    
    # extract Confidence Intervale for correlated/anticorrelated genes in cluster
    ci_vectors <- CI_vecs(M)
    
    # get submatrix of genes correlated/anticorrelated
    positive <- ci_vectors[[5]]
    negative <- ci_vectors[[6]]
    
    # get index of random subsamples from correlated/anticorrelated genes
    # avoid plots overcharged of lines
    sub.pos <- sample(dim(positive)[1], round(dim(positive)[1]/10), replace=T)
    sub.neg <- sample(dim(negative)[1], round(dim(negative)[1]/10), replace=T)
    
    # using matplot (matrix plots of lines)
    matplot(central, type = "l", col = "#56B4E9", lwd = 1, ylim=c(-5,5), xaxt= "n", yaxt= "n")
    
    # back ground preparation
    rect(xleft = c(0,3.5,6.5,7.5,8.5,9.5,10.5,13.5,16.5,19.5),
         xright =c(3.5,6.5,7.5,8.5,9.5,10.5,13.5,16.5,19.5,20.5),
         ybottom = par("usr")[2], ytop=par("usr")[3], col = (col_set2), border = "transparent") # background rectangles-sample separators
    boxplot(x=20,y=0,ylim=c(-5,5), col="white", xaxt= "n", yaxt= "n", add = TRUE) # problem with top border, adding plot to solve
    
    # plot values
    matplot(t(positive[sub.pos,]), type = "l", col = "#56B4E9", lwd = 1, ylim=c(-5,5), xaxt= "n", yaxt= "n", add = TRUE)
    matplot(t(negative[sub.neg,]), type = "l", col = "#F0E442", lwd = 1, ylim=c(-5,5), xaxt= "n", yaxt= "n", add = TRUE)
    # plot lines of confidence intervales for correlated and anticorrelated genes
    lines(x =c(1:23), y= ci_vectors[[1]] , type= "l", lty= "dotdash",lwd=2, col="#0072B2")
    lines(x =c(1:23), y= ci_vectors[[2]] , type= "l", lty= "dotdash",lwd=2, col="#0072B2")
    lines(x =c(1:23), y= ci_vectors[[3]] , type= "l", lty= "dotdash",lwd=2, col="#D55E00")
    lines(x =c(1:23), y= ci_vectors[[4]] , type= "l", lty= "dotdash",lwd=2, col="#D55E00")
    mtext(j, side= 3, line= -1.5, adj= 0.05, outer= FALSE, font= 2, cex=1) # assign number of cluster to plot
    axis(side= 2, at= seq(from= -5, to= 5, by= 2), cex.axis = 0.5, outer= T) # assign y axe scale
    axis(side= 1, at= c(1:23), labels= FALSE) # assign x axe ticks
    text(x= c(1:23), y= -6.5, tabla$Code ,xpd= TRUE, srt= 60, cex= 0.5) # assign names to the x ticks
    title(main= paste("Expression patterns within",cl,"clusters"), outer= T) # title
  }
  dev.off()
}
