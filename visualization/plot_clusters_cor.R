gc()
######################
## loading packages ##
######################
#source("https://bioconductor.org/biocLite.R")
#biocLite("sva")
library("Biobase")
library("sva")
library(viridisLite) # color scale easy to read by those with colorblindness
library(viridis) # color scale easy to read by those with colorblindness
library(gplots)
library(ggplot2)
library(methods) # Load basic functions
library(cluster) # Loads the cluster library for kmedoids method


gc()
####################
## load functions ##
####################
# Calculate confidence intervale
confidence_interval <- function(vector, interval) {
  # Standard deviation of sample
  vec_sd <- sd(vector)
  # Sample size
  n <- length(vector)
  # Mean of sample
  vec_mean <- mean(vector)
  # Error according to t distribution
  error <- qt((interval + 1)/2, df = n - 1) * vec_sd / sqrt(n)
  # Confidence interval as a vector
  result <- c("lower" = vec_mean - error, "upper" = vec_mean + error)
  return(result)
}

# function to get maximun and minimun values from data
top2 <- function(M){
  n <-  length(M)
  
  max.1 <- -Inf
  max.2 <- -Inf
  
  for(i in 1:n){
    if(M[i] > max.1){
      max.2 <- max.1
      max.1 <- M[i]
    }
    else if(M[i] > max.2){
      max.2 <- M[i]
    }
  }
  
  min.1 <- +Inf
  min.2 <- +Inf
  
  for(i in 1:n){
    if(M[i] < min.1){
      min.2 <- min.1
      min.1 <- M[i]
    }
    else if(M[i] < min.2){
      min.2 <- M[i]
    }
  }
  return(c("max.1" = max.1, "max.2" = max.2, "min.1" = min.1, "min.2" = min.2))
}


gc()
################################
##  Create colours for Plots  ##
################################

# Stablishg colour of graph; colour blind friendly palette

# funcion para crear lighter colors, para combinar con pch 21:25 y bg
lighten <- function(color, factor=1.4){
  col <- col2rgb(color)
  col <- col*factor
  col <- rgb(t(as.matrix(apply(col, 1, function(x) if (x > 255) 255 else x))), maxColorValue=255)
  col
}

# color para background gris, divisor entre muestras del mismo experimento
col_set2 = c("#FFFFFF", "#ABABAB")
col_set2 <- sapply(col_set2, lighten)
col <- colorRampPalette(col_set2)

# Blid friendly colour pallete; inner boxplot 
col_set = c("#0072B2", "#999999", "#E69F00", "#D55E00", "#F0E442", "#56B4E9", "#009E73", "#CC79A7")
palette(col_set)


gc()
#########################
## path specifications ##
#########################
# stablish your own directory
dir = "C:/Users/LANCIOS-PC/Documents/Oficina/Labo/Paper_Transcriptome/"

# Create a folder to export the results
dir.create(path=paste0(dir,"Data/new_assembly3/pval/"), showWarnings = FALSE)
dir.create(path=paste0(dir,"Grph/new_assembly3/pval/"), showWarnings = FALSE)


gc()
#################################################
##  compute intra - inter cluster correlation  ##
#################################################
# PLOT INTRA CLUSTER VALUES
# INCRIEASING NUMBER OF CLUSTERS VALUES TEND TO BY BIGGER AND BIGGER
mycl <- c(5:75)
# create a matrix to write INTRA cluster values
results.intra <- matrix(nrow = 9, ncol = length(mycl))
colnames(results.intra) <- mycl
rownames(results.intra) <- c("max inrar cor","2max intra cor",
                             "min intra cor", "2min intra cor",
                             "dif intra cor", "dif intra cor", 
                             "avg intra cor",
                             "90% intra cor","10% intra cor")
# create a matrix to write INTER cluster values
results.inter <- matrix(nrow = 9, ncol = length(mycl))
colnames(results.inter) <- mycl
rownames(results.inter) <- c("max inter cor","2max inter cor",
                             "min inter cor", "2min inter cor",
                             "dif inter cor", "dif inter cor", 
                             "avg inter cor",
                             "90% inter cor","10% inter cor")

# create a counter to pass iteration argument to marix
count = 0
# loop for loading files results.pxxx.xx
# and write average values to one single table called results ???
for (i in mycl){
  # load results.pcor output file
  y <- readRDS(file = paste0(dir, "Data/new_assembly3/pval/", "results.pcor.", i))
  
  count = count + 1
  results.intra[1, count] = top2(diag(y))[1] # max intra cor
  results.intra[2, count] = top2(diag(y))[2] # 2max intra cor
  results.intra[3, count] = top2(diag(y))[3] # min intra cor
  results.intra[4, count] = top2(diag(y))[4] # 2min intra cor
  results.intra[5, count] = mean(c(top2(diag(y))[1], top2(diag(y))[2])) # average max
  results.intra[6, count] = (top2(diag(y))[2]) - (top2(diag(y))[4]) # 2max - min inter corr
  results.intra[7, count] = mean(diag(y)) # avg inter corr
  results.intra[8, count] = confidence_interval(diag(y), 0.9)["upper"] # 90% inter corr
  results.intra[9, count] = confidence_interval(diag(y), 0.9)["lower"] # 10% inter corr
}
saveRDS(results.intra, file = paste0(dir, "Data/new_assembly3/pval/", "results.intra"))
# create a counter to pass iteration argument to marix
count = 0
# loop for loading files results.pxxx.xx
# and write average values to one single table called results
for (i in mycl){
  # load results from p valor correlation test inter - intra cluster
  y <- readRDS(file = paste0(dir,"Data/new_assembly3/pval/", "results.pcor.", i))
  
  count = count + 1
  results.inter[1, count] = top2(y[upper.tri(y)])[1] # max inter cor
  results.inter[2, count] = top2(y[upper.tri(y)])[2] # 2max inter cor
  results.inter[3, count] = top2(y[upper.tri(y)])["min.1"] # min inter cor
  results.inter[4, count] = top2(y[upper.tri(y)])["min.2"] # 2min inter cor
  results.inter[5, count] = mean(c(top2(y[upper.tri(y)])["min.1"]),(top2(y[upper.tri(y)])["min.2"])) # max - min inter cor
  results.inter[6, count] = (top2(y[upper.tri(y)])[2]) - (top2(y[upper.tri(y)])["min.2"]) # 2max - min inter cor
  results.inter[7, count] = mean(y[upper.tri(y)]) # avg inter cor
  results.inter[8, count] = confidence_interval((y[upper.tri(y)]),0.9)["upper"] # 90% inter cor
  results.inter[9, count] = confidence_interval((y[upper.tri(y)]),0.9)["lower"] # 10% inter cor
}
saveRDS(results.inter, file = paste0(dir, "Data/new_assembly3/pval/", "results.inter"))


gc()
######################################################################
##  check min intra cluster- max inter cluster euclidean distances  ##
######################################################################

ratio <- results.inter[5,]/results.intra[5,]
ratio[which.min(ratio)]

ls <- loess(ratio~mycl)
pr.loess <- predict(ls)
k = 1:length(mycl)

dev.off()
################################
##  plot intra - inter ratio  ##
################################
png(file= paste0(dir,"Grph/new_assembly3/pval/","Cluster_segregation.png"),width=15, height=7,
    units="cm", res=300, pointsize=10)
par(mar=c(1,1,1,1), oma=c(2,4,2,2), cex.lab=1.5, cex.axis=1.5, cex.main=1)

plot(x=c(1:length(mycl)), y = ratio, ylim =c(0,0.1),
     xaxt= "n", yaxt= "n", xlab= NA, ylab=NA, 
     axes=T, pch = 21, las=1, cex=0.5,
     bg=col_set2[2], col= "black")

segments(x0= c(1:length(mycl)), y0=results.intra[8,]/10,
         x1= c(1:length(mycl)), y1=results.intra[9,]/10,
         lend=1, col = col_set2[2]) # plot IC 90% and 10% inter pval
segments(x0= c(1:length(mycl)), y0=results.inter[8,]/10,
         x1= c(1:length(mycl)), y1=results.inter[9,]/10,
         lend=1, col = col_set2[2]) # plot IC 90% and 10% inter pval

points(x= c(1:length(mycl)),y=results.intra[7,]/10,
       pch=21, xlab='', ylab='', las=1, cex=0.5, xaxt='n',
       bg=col_set2[1], col= col_set[1]) # plot points of avg intra corr### between 25(5^2) and 169 (13^2)
points(x= c(1:length(mycl)),y=results.inter[7,]/10,
       pch=21, xlab='', ylab='', las=1, cex=0.5, xaxt='n',
       bg=col_set2[3], col= col_set[3]) # plot points of avg intra corr### between 25(5^2) and 169 (13^2)
lines(pr.loess~k, col=col_set[2], lwd=1)
points(x= which.min(ratio), y=ratio[which.min(ratio)],
       pch=1, las=1, cex=2, col= "red")
points(x=c(1:length(mycl)), y = ratio, 
       pch = 21, las=1, cex=0.5,
       bg=col_set2[2], col= "black")
axis(side = 2, at= seq(from= 0, to= 0.1, by= 0.01), 
     labels = seq(from= 0, to= 1, by= 0.1), cex.axis=0.8)
axis(side=1,at=(1:length(mycl)+0.1),labels=colnames(results.intra), cex.axis=0.6)
mtext("correlation" # "ratio"
      , side=2, line=3, cex.lab=1,las=3, col="black")
mtext("Optimal cluster segregation", outer = TRUE, cex = 1.5)
legend("topright", inset= 0.05,
       col = c("black", col_set[1], col_set[3]),
       legend = c("between - within ratio","max avg within-cluster", "min avg between-cluster"),
       pch=21, cex = 0.5, y.intersp = 0.8, xjust = 1,
       bty = 'n')

dev.off()
