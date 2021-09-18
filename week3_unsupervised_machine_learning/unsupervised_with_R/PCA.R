# PCA
#set working directory
url <- "https://courses.pharmb.io/apb/solu.descr.Rdata"
download.file(url, destfile = "solu.descr.Rdata")

load("solu.descr.Rdata", verbose = TRUE)

# PCA
head(solu.descr)
names(solu.descr)
nrow(solu.descr)
pairs(solu.descr ,gap = 0.3, cex = 0.5)
round(cor(solu.descr) ,2)

# Scaling the data
pca <- prcomp(solu.descr , scale. = T)
names(pca)
pca$x
summary(pca)
plot(pca)
plot(pca$x[,1:2])

# variables working
pca$rotation
pca$rotation[,1]^2
sum(pca$rotation[,1]^2)
pca$rotation[,1]
biplot(pca)
apply(solu.descr,2,var)
pca.nsc <- prcomp(solu.descr, scale. = F)
biplot(pca.nsc)

# Explanatory data analysis / Solubility data
url.s <- "https://courses.pharmb.io/apb/solubility.csv"
download.file(url.s, destfile = "solubility.csv")
solu <- read.csv("solubility.csv", stringsAsFactors = F)
str(solu)
summary(solu)
names(solu)
names(solu)[2:3] <- c("measured","predicted")
names(solu)
cor(solu$predicted ,solu$measured)
plot(solu$predicted ~solu$measured)
hist(solu$measured)
summary(solu$measured)

# Add features
RNGkind(sample.kind = "Rounding")
set.seed(123)
solu.300 <- solu[sample(nrow(solu),size = 300),]
hist(solu.300$measured)
solu.300.mols <- parse.smiles(solu.300$SMILES)
solu.300.descr <- eval.desc(solu.300.mols ,
                              c("org.openscience.cdk.qsar.descriptors.molecular.TPSADescriptor",
                                "org.openscience.cdk.qsar.descriptors.molecular.HBondDonorCountDescriptor",
                                "org.openscience.cdk.qsar.descriptors.molecular.HBondAcceptorCountDescriptor",
                                "org.openscience.cdk.qsar.descriptors.molecular.CarbonTypesDescriptor",
                                "org.openscience.cdk.qsar.descriptors.molecular.XLogPDescriptor",
                                "org.openscience.cdk.qsar.descriptors.molecular.WeightDescriptor",
                                "org.openscience.cdk.qsar.descriptors.molecular.AromaticBondsCountDescriptor"))

str(solu.300.descr)
pairs(solu.300.descr,gap = 0.3, cex=0.3)
wss <- NULL
for (i in 1:15) wss[i] <- sum(kmeans(solu.300.descr, centers=i, nstart =15)$withinss)
plot(1:15, wss, type="b", xlab="Number of Clusters",ylab="Sum of squares")

RNGkind (sample.kind = "Rounding")
set.seed (55)
km.descr <- kmeans(solu.300.descr, centers =4, nstart =15)
pairs(solu.300.descr, col=km.descr$cluster, gap = 0.3, cex = 0.3)

hc.descr <- hclust(dist(solu.300.descr))
plot(hc.descr)
rect.hclust(hc.descr, k=4)

hc.descr.clusters <- cutree(hc.descr, k=4)
pairs(solu.300.descr, col = hc.descr.clusters, gap = 0.3, cex = 0.3)

table(km.descr$cluster)
table(hc.descr.clusters)
round(apply(solu.300.descr, 2, var) ,2)
table(solu.300.descr$C1SP1)
table(solu.300.descr$C2SP1)

pca <- prcomp(scale(solu.300.descr))
biplot(pca, cex =0.6, scale =0, ylim = c(-4,6))
pca$rotation
summary (pca)
screeplot(pca) # plot(pca) produces the same plot
