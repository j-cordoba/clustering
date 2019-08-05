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
saveRDS(data.cor, file=paste0(dir, "Data/", "data.cor"))


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
saveRDS(results.pcor.[% nofk %], file= paste0(dir,"cor/","pval/","results.pcor.",[% nofk %]))
print("results saved")


