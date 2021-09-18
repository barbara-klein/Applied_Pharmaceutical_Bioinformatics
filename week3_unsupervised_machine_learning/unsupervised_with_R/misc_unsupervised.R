# This is a my first script containing instructions prior dataset preparation 
#for unsupervised learning methods 

url.CYP <- "https://courses.pharmb.io/apb/CYP.csv"
download.file(url.CYP, destfile = "CYP.csv")
cyp <- read.csv("CYP.csv")
dim(cyp)
head(cyp)
names(cyp)

# Preparing data with 30 items: 
RNGkind(sample.kind = "Rounding")
set.seed(66)
cyp.30 <- cyp[sample(nrow(cyp) ,30) ,] # sample of 30 items randomly picked with the sample()
dim(cyp.30)
summary(cyp.30)
cyp.30$Substance.SID <- as.factor(cyp.30$Substance.SID) # the summary(cyp.30) shows the difference
str(cyp.30) # one line summary for each variable

# calculating molecular descriptors with the rcdk package:
install.packages("rcdk")
library(rcdk)
library(fingerprint)
cyp <- read.csv("CYP.csv", stringsAsFactors = F)
RNGkind(sample.kind = "Rounding")
set.seed(66)
cyp.30 <- cyp[sample(nrow(cyp) ,30) ,]
str(cyp.30)
compounds.mols <- parse.smiles(cyp.30[,2])

# convert implicit hydrogen to explicit hydrogen:
lapply(compounds.mols , convert.implicit.to.explicit)

# calculating fingerprints with the get.fingerprint() function for the whole set:
compounds.mols.fps <- lapply(compounds.mols, get.fingerprint, type='extended') # it is a list of binary fingerprints for each of the 30 items

# Tanimoto similarity is a way to measure how similar the molecules are in the space defined by the binary fingerprints, 
# Calculate Tanimoto similarity matrix for the cluster analysis:
compounds.mols.fp.sim <- fp.sim.matrix(compounds.mols.fps, method = 'tanimoto')

# Calculate the distance matrix (e.g distance between each two data points)
compounds.fp.dist <- 1 - compounds.mols.fp.sim

# calculating molecular descriptors available in CDK divided in 5 categories: 
get.desc.categories() 
# TPSADescriptor, HBondDonorCountDescriptor, HBondAcceptorCountDescriptor, WeightDescriptor, ALOGDescriptor

compounds.descrs <- eval.desc(compounds.mols, 
          c("org.openscience.cdk.qsar.descriptors.molecular.TPSADescriptor",
            "org.openscience.cdk.qsar.descriptors.molecular.HBondDonorCountDescriptor", 
            "org.openscience.cdk.qsar.descriptors.molecular.HBondAcceptorCountDescriptor", 
            "org.openscience.cdk.qsar.descriptors.molecular.WeightDescriptor", 
            "org.openscience.cdk.qsar.descriptors.molecular.ALOGPDescriptor"))
summary(compounds.descrs)
row.names(compounds.descrs) <- NULL


# Preparing data with 100 items:
cyp.100 <- cyp[sample(nrow(cyp) ,100) ,] # sample of 100 items randomly picked with the sample()
dim(cyp.100)
summary(cyp.100)
cyp.100$Substance.SID <- as.factor(cyp.100$Substance.SID) # the summary(cyp.100) shows the difference
str(cyp.100) # one line summary for each variable

# calculating molecular descriptors with the rcdk package:
RNGkind(sample.kind = "Rounding")
set.seed(66)
cyp.100 <- cyp[sample(nrow(cyp) ,100) ,]
str(cyp.100)
comp100.mols <- parse.smiles(cyp.100[,2])

# convert implicit hydrogen to explicit hydrogen:
lapply(comp100.mols , convert.implicit.to.explicit)

# calculating fingerprints with the get.fingerprint() function for the whole set:
comp100.mols.fps <- lapply(comp100.mols, get.fingerprint, type='extended') # it is a list of binary fingerprints for each of the 30 items

# Tanimoto similarity is a way to measure how similar the molecules are in the space defined by the binary fingerprints, 
# Calculate Tanimoto similarity matrix for the cluster analysis:
comp100.mols.fp.sim <- fp.sim.matrix(comp100.mols.fps, method = 'tanimoto')

# Calculate the distance matrix (e.g distance between each two data points)
comp100.mols.fp.dist <- 1 - comp100.mols.fp.sim

# calculating molecular descriptors available in CDK divided in 5 categories: 
# TPSADescriptor, HBondDonorCountDescriptor, HBondAcceptorCountDescriptor, WeightDescriptor, ALOGDescriptor

comp100.descrs <- eval.desc(comp100.mols, 
                              c("org.openscience.cdk.qsar.descriptors.molecular.TPSADescriptor",
                                "org.openscience.cdk.qsar.descriptors.molecular.WeightDescriptor", 
                                "org.openscience.cdk.qsar.descriptors.molecular.ALOGPDescriptor"))
summary(comp100.descrs)
row.names(comp100.descrs) <- NULL

# saving unsupervised Rdata
save(file = "unsupervised.Rdata", compounds.fp.dist, compounds.descrs, comp100.descrs)

