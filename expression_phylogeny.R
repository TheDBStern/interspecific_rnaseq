library(matrixStats)
library(lsa)
library(ape)
library(parallel)
library(picante)

#############
# Estimate NJ and bootstrapped NJ trees
#############

# read in data
dat <- read.table('data/orthogroups.TMM.EXPR.matrix',sep='\t', header=T, row.names=1)
# square-root transform
dat <-as.matrix(sqrt(dat))

# take the mean expression value for each species. Unfortunately this code is very specific to this dataset
ccryp <- dat[,1]
cdubi <- rowMeans2(dat, cols=c(2:4))
cgray <- rowMeans2(dat, cols=c(5:7))
chamu <- rowMeans2(dat, cols=c(8:10))
cnert <- rowMeans2(dat, cols=c(11:13))
crust <- rowMeans2(dat, cols=c(14:15))
cseto <- rowMeans2(dat, cols=c(16:18))
ctene <- rowMeans2(dat, cols=c(19:21))
oaust <- rowMeans2(dat, cols=c(22:24))
oinco <- rowMeans2(dat, cols=c(25:26))
pfall <- rowMeans2(dat, cols=c(27:29))
phors <- rowMeans2(dat, cols=c(30:31))
pluci <- rowMeans2(dat, cols=c(32:33))
ppall <- dat[,34]

#create a new data matrix from the species means
dat <- data.frame(CCRYP=ccryp,CDUBI=cdubi, CGRAY=cgray,CHAMU=chamu,CNERT=cnert,CRUST=crust,CSETO=cseto,CTENE=ctene,OAUST=oaust,OINCO=oinco,PFALL=pfall,PHORS=phors,PLUCI=pluci,PPALL=ppall)
dat <- as.matrix(dat)
# create a distance matrix using pairwise peason distances
mat <- as.dist(1-cor(dat,method='pearson'))
# generate the NJ tree
mynj <- nj(mat)

# generate the 10000 boostrapp NJ trees
dat <- t(dat)
bstrees <- boot.phylo(phy = mynj, x = dat, FUN = function(xx) nj(1-cor(t(xx),method='pearson')), B = 10000, trees=T, mc.cores=4)$trees

write.tree(bstrees,'boot.pearson.species.tre')