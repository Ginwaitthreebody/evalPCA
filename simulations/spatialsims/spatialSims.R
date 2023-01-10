
## a function that wraps the snpStats plink file reader
source("http://popgen.dk/albrecht/open/online.R")
palette(RColorBrewer::brewer.pal("Set1", n=3))


f <- "genos/genosSampleContinuousVar"

pl <- plink(f)

###########################################
################### PCA ###################
###########################################

## normalise genotypes with frequency
avg<-colMeans(pl$geno)
g<-t(pl$geno)-avg

## decompose G=UDV 
udv<-svd(g)
nInd <- ncol(g)

## plot
outpng <- "spatialSimPCA12.png"
bitmap(outpng, h=4,w=4,res=300)
plot(udv$v[,1:2], xlab="PC1", ylab="PC2", pch=21, bg=2)
dev.off()


## use top k PCA to recontruct individual allele frequencies
k <- 2
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
outpng <- "spatialSimPCA12_genetic.png"
bitmap(outpng, h=4,w=4,res=300)
plot(udv$v[,1:2], xlab="PC1", ylab="PC2", pch=21, bg=2)
dev.off()


## use top k PCA to recontruct individual allele frequencies
k <- 2
SIGMA <-diag(udv$d)
rawPI <- udv$u[,1:k]%*%SIGMA[1:k,1:k]%*%t(udv$v[,1:k])
piPCA_gen <- (rawPI*sd+avg)/2





###############################################
############### ADMIXTURE #####################
###############################################

## path to results
ffile <- "genosSampleContinuousVar.2.P"
qfile <- "genosSampleContinuousVar.2.Q"

## read data
f <- as.matrix(read.table(ffile))
q <- as.matrix(read.table(qfile))


## get individual allele freq from ADMIXTURE results
piADMIX <-t(t(q))%*%t(f)

outpng <- "spatialSimAdmixK2.png"
bitmap(outpng, h=4,w=8,res=300)
barplot(t(q), border=NA, space=0, col=2:3, ylab="Admixture proportions", cex.axis=1.5, cex.lab=1.5)
dev.off()

################################
###### ADMIXUTRE WIHT K=3 ###### 
################################

## path to results
ffile <- "admixture/genosSampleContinuousVar.3.P"
qfile <- "admixture/genosSampleContinuousVar.3.Q"

## read data
fk3 <- as.matrix(read.table(ffile))
qk3 <- as.matrix(read.table(qfile))


## get individual allele freq from ADMIXTURE results
piADMIXk3 <-t(t(qk3))%*%t(fk3)

outpng <- "spatialSimAdmixK3.png"
bitmap(outpng, h=4,w=8,res=300)
barplot(t(qk3), border=NA, space=0, col=1:3, ylab="Admixture proportions", cex.axis=1.5, cex.lab=1.5)
dev.off()


save(udv, piPCA, udv_gen, piPCA_gen, f, q, piADMIX, fk3, qk3, piADMIXk3, file="spatial_simulation.RData")

