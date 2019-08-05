##########################
##  compute correlation ##
##########################

# Calculate correlation distance square for each pair of genes, across samples.
data.cor <- as.matrix(1-(cor(t(data), method="pearson"))^2)
saveRDS(data.cor, file=paste0(dir, "Data/", "data.cor"))

k.values <- [% nofk %]

# calculate kmeans for k = nofk variable
km <- pam(data.per, k = k.values, diss = TRUE)


# define method, (max) number of clusters for writing output
# the name of the graphs and other features of graph colours
method = "kmedoids" # "kmedoids" the method used to clusterize our data
cl = max(km$clustering)
mycl = km$clustering

# export kmeans objet
saveRDS(km, file= paste0(dir,"cor/", "pval/" ,"km[% nofk %]"))

gc()
print("kmeans saved")

#####################################################################
##  compute p-values and correlation of genes intra/inter clusters ##
#####################################################################
# function to compute euclidean distance between cluster matrix and central cluster vector (cluster matrix's mean)
#mydist <- function(mat, vec){return(apply((mat - vec)^2, 1, sum))}

# create a table for results base on number of clusters
#results.pval.[% nofk %] = matrix(nrow = cl, ncol = cl)
results.pcor.[% nofk %] = matrix(nrow = cl, ncol = cl)

# Run a loop to extract a matrix made of each defined cluster
for (m1 in 1:cl){
  M1 <- matrix(data[km$clustering == m1,], ncol= ncol(data))
 # vec1 <- colMeans(M1)
 # cent.M1 <- M1[which.min(mydist(M1, vec1)),]

  for (m2 in m1:cl){
    M2 <- matrix(data[km$clustering == m2,], ncol= ncol(data))
   # vec2 <- colMeans(M2) # mean vector for centroid distance calculation
   # cent.M2 <- M2[which.min(mydist(M2, vec2)),] # calculate centroid point of cluster represented in M2

    # Create empty matrix for correlations and p.values
    p.cor = matrix(nrow = nrow(M1), ncol = nrow(M2))
    #p.val = matrix(nrow = nrow(M1), ncol = nrow(M2))

    # Run a loop to calculate p.values of each gene - gene
    for (i in 1:nrow(M1)){
      for (j in 1:nrow(M2)){
        p.cor[i,j] <- abs(cor(M1[i,], M2[j,]))
        #p.val[i,j]<- t.test(M1[i,], M2[j,])$p.value
      }
    }


    if (m1 == m2) { # if the cluster matrix that we are comparing are the same m1 = m2

      if (nrow(M1) < 2) { # in case of 1 row matrix, upper.tri functions doesn't work, so unique value = 1
        results.pcor.[% nofk %][m1,m2] <- 1
     #   results.pval.[% nofk %][m1,m2] <- 1
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
saveRDS(results.pcor.[% nofk %], file= paste0(dir,"cor/","pval/","results.pcor.",[% nofk %]))
print("results saved")


