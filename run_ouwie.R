##############
library(OUwie)
library(parallel)

####
# Theta shifts (LRT)
####
outfile <- "OUWie_results.theta.csv"
craytree <- read.tree('data/crayfish.nodelabels.tre')
tlist <- list.files(path = 'data/ouwie_dat',pattern='.sqrt.csv')

run.ouwie <- function(genedat){
	gene <- strsplit(genedat,'\\.')[[1]][1]
	print(gene)
	dat <- read.csv(paste('data/ouwie_dat',genedat,sep='/'), header=T)
	mod <- OUwie(craytree, dat, model='OU1', mserr="known")
	modm <- OUwie(craytree, dat, model='OUM', mserr="known")
	lrt <- 2 * (modm$loglik - mod$loglik)
	df <- data.frame(gene, lrt, pchisq(lrt,1, lower.tail = FALSE))
    write.table(df, file=outfile, row.names=FALSE, col.names=FALSE, sep=",", append=TRUE)
}

subset <- 1:3560
junk <- mclapply(tlist[subset], run.ouwie, mc.cores=16)

####
# Alpha estimates
####

### bootstrap runs
run.ouwie <- function(genedat){
	gene <- strsplit(genedat,'\\.')[[1]][1]
	print(gene)
	dat <- read.csv(paste('data/ouwie_dat',genedat,sep='/'), header=T)
	modma <- OUwie(craytree, dat, model='OUMA', mserr="known")
	boot.reps <- OUwie.boot(craytree,dat,model="OUMA", nboot=1000, mserr="known", alpha=modma$solution[1,], sigma.sq=modma$solution[2,],theta=modma$theta[,1], theta0=modma$theta[2,1])
	wiltest <- wilcox.test(boot.reps[,1],boot.reps[,2], alternative='less')
	df <- data.frame(gene, median(boot.reps[,1]), median(boot.reps[,2]), wiltest$p.value)
    write.table(df, file=outfile, row.names=FALSE, col.names=FALSE, sep=",", append=TRUE)
}

subset <- 1:3560
junk <- mclapply(tlist[subset], run.ouwie, mc.cores=16)

####
# Sigma^2 estimates
####

### bootstrap runs
run.ouwie <- function(genedat){
	gene <- strsplit(genedat,'\\.')[[1]][1]
	print(gene)
	dat <- read.csv(paste('data/ouwie_dat',genedat,sep='/'), header=T)
	modma <- OUwie(craytree, dat, model='OUMV', mserr="known")
	boot.reps <- OUwie.boot(craytree,dat,model="OUMV", nboot=1000, mserr="known", alpha=modma$solution[1,], sigma.sq=modma$solution[2,],theta=modma$theta[,1], theta0=modma$theta[2,1])
	wiltest <- wilcox.test(boot.reps[,3],boot.reps[,4], alternative='less')
	df <- data.frame(gene, median(boot.reps[,3]), median(boot.reps[,4]), wiltest$p.value)
    write.table(df, file=outfile, row.names=FALSE, col.names=FALSE, sep=",", append=TRUE)
}

subset <- 1:3560
junk <- mclapply(tlist[subset], run.ouwie, mc.cores=16)