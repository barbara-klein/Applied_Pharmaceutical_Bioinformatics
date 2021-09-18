# load("~/Documents/education/UU_pahrma_bioinformatics/Applied Pharmaceutical Bioinformatics/week3_unsupervised_machine_learning/unsupervised.Rdata")
getwd()
list.files()
# [1] "CYP.csv"                "descriptors.Rdata"      "ESLII_print12_toc.pdf" 
# [4] "PB_FirstStepsWithR.Rmd" "PCA_with_R"             "unsupervised_with_R"   
# [7] "unsupervised.pdf" 

list.files('unsupervised_with_R')
# relative path
load("unsupervised_with_R/unsupervised.Rdata")

# hierarchical clustering
hc.fp <- hclust(as.dist(compounds.fp.dist))
plot(hc.fp)
rect.hclust(hc.fp ,k=4)

# dendrogram plotting
dend <- as.dendrogram(hc.fp)
str(dend)
plot(dend)
plot(dend, type="triangle")
dend.h0.8 <- cut(dend, h=0.8)
dend.h0.8$upper
dend.h0.8$lower
plot(dend.h0.8$upper , main = "Upper tree of cut at h=0.8")
plot(dend.h0.8$lower [[1]] , main = "First branch of lower tree with cut at h=0.8")
str( dend.h0.8$lower[[1]])
length( dend.h0.8$lower )
str( dend.h0.8$lower[[7]] )

#more cluster visualisation
ddim <- cmdscale(as.dist(compounds.fp.dist), k = 2)
clusters.fp <- cutree(hc.fp , k=5)
clusters.fp
install.packages ("cluster")
library (cluster)
clusplot(ddim ,clusters.fp , color=T, labels =4, lines =0)

#Hierarchical clustering on molecular descriptors
hc.descr <- hclust(dist(compounds.descrs))
plot(hc.descr)
rect.hclust(hc.descr ,k=4)
plot(compounds.descrs$MW ,compounds.descrs$TopoPSA)
clusters.descrs <-cutree(hc.descr , k=5)
plot(compounds.descrs$MW ,compounds.descrs$TopoPSA , col=clusters.descr)

#k-means clustering
RNGkind (sample.kind = "Rounding")
set.seed (11) 
# set.seed() kills the randomness and make the algorythm deterministic in this case
kmeans.fp <- kmeans(compounds.fp.dist , centers = 5)
kmeans.fp
kmeans.fp$cluster
table(kmeans.fp$cluster)
kmeans.fp$betweenss /kmeans.fp$totss
kmeans.fp <- kmeans(compounds.fp.dist ,centers = 5, nstart = 15)
table(kmeans.fp$cluster)

#Visualisation of clusters
ddim.fp <- cmdscale(as.dist(compounds.fp.dist), k = 2)
library (cluster)
clusplot(ddim.fp ,kmeans.fp$cluster , color=T, labels =4, lines =0)
RNGkind (sample.kind = "Rounding")
set.seed (11)
kmeans.descr <- kmeans(compounds.descrs ,centers = 5, nstart =15)
ddim.descr <- cmdscale(dist(compounds.descrs), k = 2)
clusplot(ddim.descr ,kmeans.descr$cluster , color=T, labels =4, lines =0)
plot(compounds.descrs[c("MW","TopoPSA")], col=kmeans.descr$cluster)
points(kmeans.descr$centers[,c("MW","TopoPSA")], col =1:5, pch=8, cex =2)

# Choice of k
kmeans.descr$withinss
wss <- NULL
for (i in 1:15) wss[i] <- sum(kmeans(compounds.descrs ,centers=i)$withinss)
plot (1:15 , wss , type="b", xlab="Number of Clusters",ylab="Within groups sum of squares")

# Outliers detection
boxplot (scale(comp100.descrs))
boxplot(comp100.descrs)
par(mfrow = c(1, 3))
boxplot (comp100.descrs$TopoPSA)
boxplot (comp100.descrs$MW)
boxplot (comp100.descrs$AMR)
boxplot.stats(comp100.descrs$TopoPSA)
boxplot.stats(comp100.descrs$TopoPSA)$out

# finding index of the (TPSA outliers) observations 
tpsaOut <- boxplot.stats(comp100.descrs$TopoPSA)$out
tpsaOut <- which(comp100.descrs$TopoPSA %in% tpsaOut)
tpsaOut

# finding index of the (MW outliers) observations
MWOut <- which(comp100.descrs$MW %in% boxplot.stats(comp100.descrs$MW)$out)
MWOut

# visualizing outliers
par(mfrow = c(1,1))
plot(comp100.descrs$MW , comp100.descrs$TopoPSA)
points(comp100.descrs[tpsaOut ,c("MW","TopoPSA")], col="red", pch="+", cex =2)

plot(comp100.descrs$MW , comp100.descrs$TopoPSA)
points(comp100.descrs[union(tpsaOut ,MWOut),c("MW","TopoPSA")], col="red", pch="+", cex =2)

# Outlier detection by clustering (detect outliers is by performing k-means clustering)
RNGkind (sample.kind = "Rounding")
set.seed (11)
km.d <- kmeans(comp100.descrs ,centers = 5)
centers.d <- km.d$centers[km.d$cluster ,]
dist.d <- sqrt(rowSums (( comp100.descrs -centers.d)^2))
km.d$centers
km.d$cluster
centers.d

outliers.d <- order(dist.d,decreasing = T)[1:5]
outliers.d
comp100.descrs[outliers.d,]

plot(comp100.descrs$MW , comp100.descrs$TopoPSA , col=km.d$cluster , cex =0.3)
points(km.d$centers[,c("MW","TopoPSA")], col = 1:5, cex =1.5, pch =8)
points(comp100.descrs[outliers.d,c("MW","TopoPSA")], col=6, cex=1, pch =16)
