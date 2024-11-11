library("NbClust")

otu <- read.csv("out_decontam.csv", header = TRUE, sep = ",", quote = "", row.names = 1)
otu <- otu[,-c(28:ncol(otu))]
otu <- otu[,-c(2,10)]
group <- read.csv("group.csv", header = TRUE, sep = ",", quote = "", row.names = 1)
otu <- t(otu)

rownames(otu) <- gsub("CCeC", "CON", rownames(otu))
rownames(otu) <- gsub("EJCeC", "MET", rownames(otu))
rownames(otu) <- gsub("MCeC", "T2DM", rownames(otu))
rownames(group) <- gsub("CCeC", "CON", rownames(group))
rownames(group) <- gsub("EJCeC", "MET", rownames(group))
rownames(group) <- gsub("MCeC", "T2DM", rownames(group))

otu.scaled <- scale(otu)##标准化数据
dist1<-dist(otu.scaled)###计算欧氏距离
dist2<-dist(otu)###计算欧氏距离
fit1<-hclust(dist2,method = "average")

plot(fit1,hang = -1,cex=.8,main = "title")

devAskNewPage(ask = T)
nc<-NbClust(otu.scaled,distance = "euclidean",min.nc = 2,max.nc = 15,
            method = "average")
table(nc$Best.n[1,])
barplot(table(nc$Best.n[1,]),xlab = "No. of cluster")



library(factoextra)
library(igraph)
fviz_dend(fit1)

fviz_dend(fit1,k=3,rect =T,rect_fill = T,
          rect_border = c("#2E9FDF", "#00AFBB", "#E7B800", "#FC4E07"))##K为聚类个数，rect_border为区域颜色填充

fviz_dend(fit1,k=3,rect =T,rect_fill = T,horiz = TRUE,
          rect_border = c("#2E9FDF", "#E7B800", "#FC4E07"))



