##############
library(OUwie)
library(parallel)

####
# Theta shifts
####
outfile <- "OUWie_results.theta.csv"
craytree <- read.tree('data/crayfish.nodelabels.tre')
tlist <- list.files(path = '../ouwie_dat',pattern='.sqrt.csv')

run.ouwie <- function(genedat){
	gene <- strsplit(genedat,'\\.')[[1]][1]
	print(gene)
	dat <- read.csv(paste('../ouwie_dat',genedat,sep='/'), header=T)
	mod <- OUwie(craytree, dat, model='OU1', mserr="known")
	modm <- OUwie(craytree, dat, model='OUM', mserr="known")
	lrt <- 2 * (modm$loglik - mod$loglik)
#	boot.reps <- OUwie.boot(craytree,dat,model="OUM", nboot=1000, mserr="known", alpha=modm$solution[1,], sigma.sq=modm$solution[2,],theta=modm$theta[,1], theta0=modm$theta[2,1])
#	wiltest <- wilcox.test(boot.reps[,5],boot.reps[,6])
#	df <- data.frame(gene, lrt, pchisq(lrt,1, lower.tail = FALSE), median(boot.reps[,5]), median(boot.reps[,6]), wiltest$p.value)
	df <- data.frame(gene, lrt, pchisq(lrt,1, lower.tail = FALSE))
    write.table(df, file=outfile, row.names=FALSE, col.names=FALSE, sep=",", append=TRUE)
}

subset <- 1:566
#junk <- mclapply(tlist[subset], run.ouwie, mc.cores=16)
junk <- lapply(tlist, run.ouwie)

outfile <- "OUWie_results.csv"
craytree <- read.tree('data/crayfish.nodelabels.tre')

tlist <- list.files(path = '../ouwie_dat',pattern='.sqrt.csv')


### bootstrap runs
run.ouwie <- function(genedat){
	gene <- strsplit(genedat,'\\.')[[1]][1]
	print(gene)
	dat <- read.csv(genedat, header=T)
	modma <- OUwie(craytree, dat, model='OUMA', mserr="known")
	boot.reps <- OUwie.boot(craytree,dat,model="OUMA", nboot=1000, mserr="known", alpha=modma$solution[1,], sigma.sq=modma$solution[2,],theta=modma$theta[,1], theta0=modma$theta[2,1])
	wiltest <- wilcox.test(boot.reps[,1],boot.reps[,2], alternative='less')
	df <- data.frame(gene, median(boot.reps[,1]), median(boot.reps[,2]), wiltest$p.value)
    write.table(df, file=outfile, row.names=FALSE, col.names=FALSE, sep=",", append=TRUE)
}

subset <- 1:566
junk <- mclapply(tlist[subset], run.ouwie, mc.cores=16)

mod <- OUwie(craytree, dat, model='OU1', mserr="known")
modm <- OUwie(craytree, dat, model='OUM', mserr="known")
lrt <- 2 * (modm$loglik - mod$loglik)
boot.reps <- OUwie.boot(craytree,dat,model="OUM", nboot=100, mserr="known", alpha=modm$solution[1,], sigma.sq=modm$solution[2,],theta=modm$theta[,1], theta0=modm$theta[2,1])
wiltest <- wilcox.test(boot.reps[,5],boot.reps[,6])


## lrt runs
run.ouwie <- function(genedat){
	gene <- strsplit(genedat,'\\.')[[1]][1]
	print(gene)
	dat <- read.csv(paste('../ouwie_dat',genedat,sep='/'), header=T)
	modm <- OUwie(craytree, dat, model='OUM', mserr="known")
	modma <- OUwie(craytree, dat, model='OUMA', mserr="known")
	lrt <- 2 * (modma$loglik - modm$loglik)
	df <- data.frame(gene, modma$solution[1,][1],modma$solution[1,][2],lrt,pchisq(lrt,1, lower.tail = FALSE))
    write.table(df, file=outfile, row.names=FALSE, col.names=FALSE, sep=",", append=TRUE)
}

junk <- lapply(tlist, run.ouwie)

