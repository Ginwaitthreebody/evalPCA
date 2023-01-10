
## a function that wraps the snpStats plink file reader
source("http://popgen.dk/albrecht/open/online.R")
palette(RColorBrewer::brewer.pal("Set1", n=3))




f <- "genos/sims_k2branch03_Nsites5e5"

pl <- plink(f)

###########################################
################### PCA ###################
###########################################

## normalise genotypes with frequency
avg1<-colMeans(pl$geno)

# remove fixed sites
k <- avg1>0 & avg1<2
avg <- avg1[k]

pl$geno <- pl$geno[,k]
g<-t(pl$geno)-avg

## decompose G=UDV 
udv<-svd(g)
nInd <- ncol(g)

## plot
outpng <- "recadmixPCA12.png"
bitmap(outpng, h=4,w=4,res=300)
plot(udv$v[,1:2], xlab="PC1", ylab="PC2", pch=21, bg=2)
dev.off()


## use top k PCA to recontruct individual allele frequencies
k <- 1
A<-cbind(avg, t(t(udv$u[,1:k])*udv$d[1:k]))
S<-cbind(1,udv$v[,1:k])
piPCA<- t(S%*%t(A)/2)



###########################################
################### genetic PCA ###########
###########################################

## normalise genotypes with frequency
avg<-colMeans(pl$geno)
sd <- sqrt(avg*(1-avg/2))
g<- (t(pl$geno)-avg) / sd


## decompose G=UDV 
udv_gen <-svd(g)
nInd <- ncol(g)

## plot
outpng <- "recadmixPCA12_genetic.png"
bitmap(outpng, h=4,w=4,res=300)
plot(udv$v[,1:2], xlab="PC1", ylab="PC2", pch=21, bg=2)
dev.off()


## use top k PCA to recontruct individual allele frequencies
k <- 1
SIGMA <-diag(udv$d)
rawPI <- as.matrix(udv$u[,1:k])%*%SIGMA[1:k,1:k]%*%t(udv$v[,1:k])
piPCA_gen <- (rawPI*sd+avg)/2





###############################################
############### ADMIXTURE #####################
###############################################

## path to results
ffile <- "sims_k2branch03_Nsites5e5_bestlike.2.P"
qfile <- "sims_k2branch03_Nsites5e5_bestlike.2.Q"

## read data
f <- as.matrix(read.table(ffile))
q <- as.matrix(read.table(qfile))


## get individual allele freq from ADMIXTURE results
piADMIX <-t(t(q))%*%t(f)

outpng <- "recadmixAdmixK2.png"
bitmap(outpng, h=4,w=8,res=300)
barplot(t(q), border=NA, space=0, col=2:3, ylab="Admixture proportions", cex.axis=1.5, cex.lab=1.5)
dev.off()


save(udv, piPCA, udv_gen, piPCA_gen, f, q, piADMIX, file="recadmix_simulation.RData")

