##########################
##  compute correlation ##
##########################

# Calculate correlation distance square for each pair of genes, across samples.
data.cor <- as.matrix(1-(cor(t(data), method="pearson"))^2)
saveRDS(data.cor, file=paste0(dir, "Data/", "data.cor"))



