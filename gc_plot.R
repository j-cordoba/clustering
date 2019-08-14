library(hexbin)
dir="C:/Users/LANCIOS-PC/Documents/Oficina/Labo/Paper_Transcriptome/Data/new_assembly3/"

exp = c("EML","DRK","KLY","UNK","YSD","KLY.cl")

myList = list()

dev.off()

for (i in 1:length(exp)) {
  myList[[i]] <- read.csv(file=paste0(dir,exp[i],".txt"), sep="\t")
  hist(myList[[i]]$GC, panel.first=grid(), xlim=c(0,1), yaxt= "n", main = exp[i], col="cornflowerblue")
  abline(v=mean(myList[[i]]$GC),col="red")
  
  plot(hexbin(myList[[i]]$Length~myList[[i]]$GC, xbins=500, shape = 0.8),
                         yaxt="n", style= "color",
                         colramp = function(n) gray.colors(15),
                         legend=0, xlab="", ylab="")
  
}

grid.newpage()
pushViewport(viewport(layout=grid.layout(2, 3)))

for (k in 1:6) {
  
  pushViewport(viewport(layout.pos.col = as.integer((k-1) / 2) + 1,
                        layout.pos.row = (k-1) %% 2 + 1))
  plot(hexbin(RG$R[,k],RG$G[,k],xbin=50), ..., newpage = FALSE)
  popViewport(1)
  
}
