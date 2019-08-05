#########################
## path specifications ##
#########################
# stablish your own directory
dir = "/media/vol2/home/jcordoba/Transcriptomics/ALL/new_assembly3/okay/NO_FINAL_set/publicset/abundance/"

# Create a folder to export the results
dir.create(path=paste0(dir,"cor/"), showWarnings = FALSE)
dir.create(path=paste0(dir,"cor/","pval/"), showWarnings = FALSE)
dir.create(path=paste0(dir,"cor/","pval/","grph/"), showWarnings = FALSE)


gc()
##########################
##  compute correlation ##
##########################
# Calculate correlation distance square for each pair of genes, across samples.
data.cor <- as.matrix(1-(cor(t(data), method="pearson"))^2)
saveRDS(data.cor, file=paste0(dir, "cor/", "data.cor"))


gc()
########################
##  compute clusters  ##
########################
# Build clusters based on correlation distance square for each pair of genes, across samples.
k.values <- [% nofk %] # number of clusters, explore different number

# calculate kmeans for k = nofk variable
km <- pam(data.cor, k = k.values, diss = TRUE)
# export kmeans objet
saveRDS(km, file= paste0(dir,"cor/", "pval/" ,"km[% nofk %]"))
print("kmeans saved")

# define method, (max) number of clusters for writing output
method = "kmedoids" # "kmedoids" the method used to clusterize our data
cl = max(km$clustering)
mycl = km$clustering


gc()
######################################################
##  compute gene - genes pairs cluster correlation  ##
######################################################
# create a table for results base on number of clusters
results.pcor.[% nofk %] = matrix(nrow = cl, ncol = cl)

# Run a loop to extract a matrix made of each defined cluster
for (m1 in 1:cl){
  M1 <- matrix(data[km$clustering == m1,], ncol= ncol(data))
  # function to compute euclidean distance between cluster matrix and central cluster vector (cluster matrix's mean), only for kmeans
  # mydist <- function(mat, vec){return(apply((mat - vec)^2, 1, sum))}
  # vec1 <- colMeans(M1)
  # cent.M1 <- M1[which.min(mydist(M1, vec1)),]
  for (m2 in m1:cl){
    M2 <- matrix(data[km$clustering == m2,], ncol= ncol(data))
    # vec2 <- colMeans(M2) # mean vector for centroid distance calculation
    # cent.M2 <- M2[which.min(mydist(M2, vec2)),] # calculate centroid point of cluster represented in M2

    # Create empty matrix for correlations
    p.cor = matrix(nrow = nrow(M1), ncol = nrow(M2))

    # Run a loop to calculate correlation of each gene - gene pairs
    for (i in 1:nrow(M1)){
      for (j in 1:nrow(M2)){
        p.cor[i,j] <- abs(cor(M1[i,], M2[j,]))
      }
    }

    if (m1 == m2) { # if the cluster matrix that we are comparing are the same m1 = m2
      if (nrow(M1) < 2) { # in case of 1 row matrix (matrix coherced to vector, upper.tri functions doesn't work, so unique value = 1
        results.pcor.[% nofk %][m1,m2] <- 1
      } else { # in case of square matrix (because m1 = m2), upper.tri functions gives mean
        results.pcor.[% nofk %][m1,m2] <- mean(p.cor[upper.tri(p.cor)])
      }
    } else {

      # in case of two different matrix m1 != m2 , total mean is calculated
      results.pcor.[% nofk %][m1,m2] <- abs(cor((data[km$medoids[m1],]), (data[km$medoids[m2],])))
    }
  }
} print("loop done")

# write results.pcor output file
colnames(results.pcor.[% nofk %]) <- seq(from=1, to=cl)
rownames(results.pcor.[% nofk %]) <- seq(from=1, to=cl)
saveRDS(results.pcor.[% nofk %], file= paste0(dir,"cor/", "pval/", "results.pcor.",[% nofk %]))
print("results saved")


gc()
#################################################
##  compute intra - inter cluster correlation  ##
#################################################
# create a matrix to write INTRA cluster values
results.intra <- matrix(nrow = 9, ncol = length(mycl))
colnames(results.intra) <- mycl
rownames(results.intra) <- c("max inter cor","2max inter cor",
                             "min inter cor", "2min inter cor",
                             "dif inter cor", "dif inter cor", 
                             "avg inter cor",
                             "90% inter cor","10% inter cor")
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
  y <- readRDS(file = paste0(dir, "cor/", "pval/", "results.pcor.", i))
  
  count = count + 1
  results.intra[1, count] = top2(diag(y))[1] # max inter cor
  results.intra[2, count] = top2(diag(y))[2] # 2max inter cor
  results.intra[3, count] = top2(diag(y))[3] # min inter cor
  results.intra[4, count] = top2(diag(y))[4] # 2min inter cor
  results.intra[5, count] = (top2(diag(y))[1]) - (top2(diag(y))[3]) # max - min inter corr
  results.intra[6, count] = (top2(diag(y))[2]) - (top2(diag(y))[4]) # 2max - min inter corr
  results.intra[7, count] = mean(diag(y)) # avg inter corr
  results.intra[8, count] = confidence_interval(diag(y), 0.9)["upper"] # 90% inter corr
  results.intra[9, count] = confidence_interval(diag(y), 0.9)["lower"] # 10% inter corr
}

# create a counter to pass iteration argument to marix
count = 0
# loop for loading files results.pxxx.xx
# and write average values to one single table called results
for (i in mycl){
  # load results from p valor correlation test inter - intra cluster
  y <- readRDS(file = paste0(dir,"Data/scaled/", "results.pcor.", i))
  
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


gc()
######################################################################
##  check min intra cluster- max inter cluster euclidean distances  ##
######################################################################
A <- mean(c((1-results.intra[1,]),(1-results.intra[2,])))/sqrt(1-results.inter[5,])
# B <- sqrt(1- results.inter[5,])/mean(c((1-results.intra[1,]),(1-results.intra[2,])))


