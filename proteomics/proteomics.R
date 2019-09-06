source("https://bioconductor.org/biocLite.R")
biocLite("limma")
biocLite("Biobase")

library(Biobase)
library(RColorBrewer)
library(gplots)
library(limma)
library(edgeR)
library(ggplot2)
library(reshape2)
library(data.table)

####################################
######  READ PROTEOMIC DATA   ######
####################################

# Import data and convert to matrix

dir = "/media/javier/Datas/Proteomics/Euglena/FINAL/new_assembly2/"
#dir = "/media/javier/STORE N GO/"
#dir = "D:/Rstudio/"

##############################################
### CONDENSE DUPLICATE ID CALCULATING MEAN ###
##############################################

## how to condense rows from a table that has the same identifier. We want to calculate the mean of the duplicate ID
## we invoke the library (from data.camp)

## we read our table
dat <- read.table(file=paste0(dir,'ratio'), sep='\t', header=TRUE)
## we set our key term, such as ID (in this case)
keys <- colnames(dat[1])
## we transform our data into a table
X <- as.data.table(dat)
## we calculate the mean of duplicate rows (!!!how it works??)
ratio <- X[,lapply(.SD,mean),keys]
## we transform the output into dataframe
proteins_cond.rat <- as.data.frame(ratio)
## we export our new table
write.table(proteins_cond.rat, file=paste0(dir,'proteins_cond.rat'), quote=FALSE, row.names=FALSE,  sep='\t')

# Para excluir lineas (por que corresponden a contaminantes), necesitamos una lista de ID a excluir y la tabla general con IDs buenos y malos.
# awk -F "\t" 'FNR==NR{a[$1]=$1;next}!($1 in a){print $0}' conta.idl proteins.rat > proteins_clean.rat

###########################
### LOG TRANSFORMATION  ###
###########################

# data preparation for statistics and graphical representation
# import data
pro_rat <- read.table(file=paste0(dir,'proteins_cond_clean.rat'), sep='\t', header=TRUE)
class(pro_rat)
head(pro_rat)
rownames(pro_rat) <- pro_rat$Main.Accession # add the names of proteins with row name

pro_mtx <- as.matrix(pro_rat[,2:5]) # extract numerical values as a matrix
pro_mtx[pro_mtx == 0] <- NA # transfrom value 0 to NA, otherwise log transformation will give infinite value of log 0
pro_log <- log2(pro_mtx) # log transfromation, now the ratio data looks symetrical
pro_log[is.na(pro_log)] <- 0 # replace NA values by 0
colnames(pro_log[1]) <- c("ID")
write.table(pro_log, file='/media/javier/Datas/Proteomics/Euglena/ALLclean_MTCwt/new_assembly2/proteins_cond_clean_log.rat', quote=FALSE, row.names=TRUE,  sep='\t')

heatmap(pro_log)


######################
###  CORRELATION   ###
######################

cormat <- round(cor(pro_log),2)
library(reshape2)
melted_cormat <- melt(cormat)
head(melted_cormat)
library(ggplot2)
ggplot(data = melted_cormat, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile()
# Get lower triangle of the correlation matrix
get_lower_tri<-function(cormat){
  cormat[upper.tri(cormat)] <- NA
  return(cormat)
}
# Get upper triangle of the correlation matrix
get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}
upper_tri <- get_upper_tri(cormat)
upper_tri
# Melt the correlation matrix
library(reshape2)
melted_cormat <- melt(upper_tri, na.rm = TRUE)
# Heatmap
library(ggplot2)
ggplot(data = melted_cormat, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Pearson\nCorrelation") +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1))+
  coord_fixed()
reorder_cormat <- function(cormat){
  # Use correlation between variables as distance
  dd <- as.dist((1-cormat)/2)
  hc <- hclust(dd)
  cormat <-cormat[hc$order, hc$order]
}
# Reorder the correlation matrix
cormat <- reorder_cormat(cormat)
upper_tri <- get_upper_tri(cormat)
# Melt the correlation matrix
melted_cormat <- melt(upper_tri, na.rm = TRUE)
# Create a ggheatmap
ggheatmap <- ggplot(melted_cormat, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Pearson\nCorrelation") +
  theme_minimal()+ # minimal theme
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1))+
  coord_fixed()
# Print the heatmap
print(ggheatmap)
ggheatmap + 
  geom_text(aes(Var2, Var1, label = value), color = "black", size = 4) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    legend.justification = c(1, 0),
    legend.position = c(0.6, 0.8),
    legend.direction = "horizontal")+
  guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
                               title.position = "top", title.hjust = 0.5))


######################
### CLUSTERIZATION ###
######################

## PER FACTOR
d <- dist(t(pro_log)) # calculate distance between samples
h <- hclust(d) # clusterize ditances
plot(h)
cutree(h, 1)

mds <- cmdscale(d)
cols <- as.factor(colnames(pro_log))
plot(mds, col=as.numeric(cols))
legend("top", levels(cols), col=seq(along=levels(cols)), pch=seq(along=levels(cols)))

## PER GENE
d <- dist(pro_log) # calculate distance between samples
h <- hclust(d) # clusterize ditances
plot(h)
cutree(h, 6)

mds <- cmdscale(d)
cols <- as.factor(colnames(pro_log))
plot(mds, col=as.numeric(cols))
legend("top", levels(cols), col=seq(along=levels(cols)), pch=levels(cols))

###############################
### HEAT MAP REPRESENTATION ###
###############################
pro_rat <- read.table(file='/media/javier/Datas/Proteomics/Euglena/ALLclean_MTCwt/new_assembly2/proteins_cond_clean.rat', sep='\t', header=TRUE)
class(pro_rat)
rownames(pro_rat) <- pro_rat$Main.Accession

dat <- as.matrix(t(pro_rat[,2:5])) # numerical rows (trasposed)

df_molten_dat <- melt(as.matrix(dat)) # reshape into dataframe
names(df_molten_dat)[c(1:2)] <- c("Exp", "Protein")


g <- ggplot(data = df_molten_dat,
            aes(x = Protein, y = Exp, fill = value)) + 
  geom_raster() +
  scale_fill_distiller(palette = "RdYlBu", trans = "log10") +
  theme(plot.title = element_text(size= 20, hjust= 0.5, face= "bold", colour= "black"),
        axis.title.x = element_text(size= 14,face= "bold",vjust = 0.2, colour = "black"),
        axis.title.y = element_text(size= 14,face= "bold",vjust = 0.2, colour = "black"),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size= 12, hjust = 1),
        #legend.title = element_text(colour="black", size=8, face="bold"),
        legend.title = element_blank(),
        legend.text = element_text(colour="black", size=12, face="bold"),
        #legend.position= c(0.2, 0.5),
        legend.position=("right"),
        #scale_fill_discrete(labels=c("dirty", "clean")),
        
        axis.line.x = element_line(colour = "black",size = 0.25,linetype = "solid"),
        axis.line.y = element_line(colour = "black",size = 0.25,linetype = "solid"),
        
        #axis.ticks = element_line(colour = "black",size = 0.25,linetype = "solid"),
        #axis.ticks.length = unit(0.25,"cm"),
        panel.background = element_rect(fill="gray95", colour = "black"),
        #panel.grid.major = element_line(colour = "grey", size= 0.25, linetype = "dashed")
        element_blank(), panel.grid.minor = element_blank()) + 
  ggtitle("Protein Expression Ratios")

ppi <- 300

png(filename="/media/javier/Datas/Transcriptomics/resumen/Graphs/pro_heatmap.png",width=15*ppi, height=10*ppi, res=ppi)
plot(g)
dev.off()

###############################
### STATISTICS AND CONTRAST ###
###############################

pro_log <- read.table(file=paste0(dir,'proteins_cond_clean_log.rat'), sep='\t', header=TRUE)
mat <- as.matrix(pro_log[,2:5])
rownames(mat) <- pro_log$ID
groups <- colnames(mat)

eset <- ExpressionSet(assayData = mat)

#design <- model.matrix(~ 0 + groups)
#design <- model.matrix(~ 0 + c(2,2,2,2,1,1,1,1,2,2,2,2,2,2,2,2))

#colnames(design) <- colnames(mat)
a <- c("AC", "AC", "TMP", "TMP")
b <- c("LL", "HL", "LL", "LL")
luz <- factor(b)
medio <- factor(a)
#design <- model.matrix()
design.m <- model.matrix(~0+medio)
design.l <- model.matrix(~0+luz)
#cont <- makeContrasts(CN-(AC+CO+MN),levels=design)


lmfit.l <- lmFit(mat, design.l)
lmfit.eb.l <- eBayes(lmfit.l)
t <- lmfit.eb.l$p.value
t <- topTable(lmfit.eb.l, confint=TRUE)

lmfit.m <- lmFit(mat, design.m)
lmfit.eb.m <- eBayes(lmfit.m)
lmfit.eb.m$p.value
topTable(lmfit.eb.m)

listLUZ <- lmfit.eb.l[["p.value"]] < 0.05

HL <- as.data.frame(which(listLUZ[,1]))
HL$HL <- rownames(HL)
HL <- as.data.frame(HL$HL)
colnames(HL) <- c("ID")

LL <- as.data.frame(which(listLUZ[,2]))
LL$LL <- rownames(LL)
LL <- as.data.frame(LL$LL)
colnames(LL) <- c("ID")

listMEDIO <- lmfit.eb.m[["p.value"]] < 0.05

AC <- as.data.frame(which(listMEDIO[,1]))
AC$AC <- rownames(AC)
AC <- as.data.frame(AC$AC)
colnames(AC) <- c("ID")

TMP <- as.data.frame(which(listMEDIO[,2]))
TMP$TMP <- rownames(TMP)
TMP <- as.data.frame(TMP$TMP)
colnames(TMP) <- c("ID")

uniprot <- read.table(file='D:\\Rstudio\\proteins.ann', sep='\t', header=TRUE)
write.table(uniprot$uniprot, file = 'D:\\Rstudio\\sig\\uniprot.txt',
            sep = '\t', col.names = F , row.names = F, quote = F)
HL.uni <- merge(HL, uniprot)
LL.uni <- merge(LL, uniprot)
TMP.uni <- merge(TMP, uniprot)
AC.uni <- merge(AC, uniprot)


write.table(HL.uni, file = 'D:\\Rstudio\\sig\\HL.uni.txt',
            sep = '\t', row.names = F, quote = F)
write.table(LL.uni, file = 'D:\\Rstudio\\sig\\LL.uni.txt',
            sep = '\t', row.names = F, quote = F)
write.table(TMP.uni, file = 'D:\\Rstudio\\sig\\TMP.uni.txt',
            sep = '\t', row.names = F, quote = F)
write.table(AC.uni, file = 'D:\\Rstudio\\sig\\AC.uni.txt',
            sep = '\t', row.names = F, quote = F)


kegg <- read.csv(file = 'D:\\Rstudio\\kegg.csv', sep = ';', header = F)

HL.keg <- merge(HL, kegg)
LL.keg <- merge(LL, kegg)
TMP.keg <- merge(TMP, kegg)
AC.keg <- merge(AC, kegg)


write.table(HL.keg, file = 'D:\\Rstudio\\sig\\HL.keg.txt',
            sep = '\t', row.names = F, quote = F)
write.table(LL.keg, file = 'D:\\Rstudio\\sig\\LL.keg.txt',
            sep = '\t', row.names = F, quote = F)
write.table(TMP.keg, file = 'D:\\Rstudio\\sig\\TMP.keg.txt',
            sep = '\t', row.names = F, quote = F)
write.table(AC.keg, file = 'D:\\Rstudio\\sig\\AC.keg.txt',
            sep = '\t', row.names = F, quote = F)

entrez <- read.csv(file = 'D:\\Rstudio\\EntrezID.csv', sep = ';', header = F)

HL.enz <- merge(HL, entrez)
LL.enz <- merge(LL, entrez)
TMP.enz <- merge(TMP, entrez)
AC.enz <- merge(AC, entrez)


write.table(HL.enz, file = 'D:\\Rstudio\\sig\\HL.enz.txt',
            sep = '\t', row.names = F, quote = F)
write.table(LL.enz, file = 'D:\\Rstudio\\sig\\LL.enz.txt',
            sep = '\t', row.names = F, quote = F)
write.table(TMP.enz, file = 'D:\\Rstudio\\sig\\TMP.enz.txt',
            sep = '\t', row.names = F, quote = F)
write.table(AC.enz, file = 'D:\\Rstudio\\sig\\AC.enz.txt',
            sep = '\t', row.names = F, quote = F)

my_palette <- colorRampPalette(c("red", "white", "green"))(n = 299)


#heatmap(leg, labRow = FALSE, col =brewer.pal(8,"Spectral") )


acetate <- merge(AC, pro_log)
#a <- 
heatmap(as.matrix(acetate[,2:5]), labRow = FALSE, col =brewer.pal(8,"Spectral") )
#ppi <- 300
#png(filename='D:\\Rstudio\\sig\\AC_heatmap.png',width=15*ppi, height=10*ppi, res=ppi)
#plot(a)
#dev.off()


minimo <- merge(TMP, pro_log)
# b <-
heatmap(as.matrix(minimo[,2:5]), labRow = FALSE, col =brewer.pal(8,"Spectral") )
#ppi <- 300
#png(filename='D:\\Rstudio\\sig\\TMP_heatmap.png',width=15*ppi, height=10*ppi, res=ppi)
#plot(b)
#dev.off()



luz <- merge(HL, pro_log)
#c <- 
heatmap(as.matrix(luz[,2:5]), labRow = FALSE, col =brewer.pal(8,"Spectral") )
#ppi <- 300
#png(filename='D:\\Rstudio\\sig\\HL_heatmap.png',width=15*ppi, height=10*ppi, res=ppi)
#plot(c)
#dev.off()

#d <- 
heatmap.2(mat, labRow = FALSE, col =brewer.pal(8,"Spectral"))
#ppi <- 300
#png(filename='D:\\Rstudio\\sig\\MAIN_heatmap.png',width=15*ppi, height=10*ppi, res=ppi)
#plot(d)
#dev.off()

anotation <- read.csv(file = 'D:\\Rstudio\\anotation_original.tab', sep = '\t', header = T)
final <- merge(uniprot, pro_log)
FINAL <- merge(final, anotation)
FINAL <- B[,c(2,1,3:16)]

sig.Acetate <- FINAL[AC$ID,]
sig.Minimo <- FINAL[TMP$ID,]
sig.Light <- FINAL[HL$ID,]

write.table(sig.Acetate, file = 'D:\\Rstudio\\sig\\sig.Acetate.txt',
            sep = '\t', row.names = F, quote = F)
write.table(sig.Minimo, file = 'D:\\Rstudio\\sig\\sig.Minimo.txt',
            sep = '\t', row.names = F, quote = F)
write.table(sig.Light, file = 'D:\\Rstudio\\sig\\sig.Light.txt',
            sep = '\t', row.names = F, quote = F)

