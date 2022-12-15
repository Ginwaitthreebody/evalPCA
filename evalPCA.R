
#build Q,F,PI,G for Scenario 1 and 2
##build Q1 for Scenario 1(specific for N_1,..,N_{K+1})
Vn<-c(0,20,40,60)
Vn<-c(0,10,30,60)

N<-Vn[length(Vn)]
M<-500000
K<-length(Vn)-1

Q<-vector()
dq<-vector()
for (p in 1:K) {
  for (i in (Vn[p]+1):Vn[p+1]) {
    for (k in 1:K) {
      Q[k]<-0;
      Q[p]<-1;
    }
    dq[i]<-data.frame(Q);
  }
}
l <- list()
for(i in 1:N){
  l[[i]] <- dq[[i]]
}
Q<-do.call(rbind,l)
Q<-t(Q)

##build Q2 for Scenario 2(specific for N_1,..,N_{K+1})
Vn<-c(0,20,40,60)
Vn<-c(0,10,30,60)

N<-Vn[length(Vn)]
M<-500000
K<-2
Q<-vector()
dq<-vector()
for (i in (Vn[1]+1):Vn[2]) {
  for (k in 1:K) {
    Q[k]<-0;
    Q[1]<-1;
  }
  dq[i]<-data.frame(Q);
}
for (i in (Vn[2]+1):Vn[3]) {
  for (k in 1:K) {
    Q[k]<-1/2;
  }
  dq[i]<-data.frame(Q);
}
for (i in (Vn[3]+1):Vn[4]) {
  for (k in 1:K) {
    Q[k]<-0;
    Q[2]<-1;
  }
  dq[i]<-data.frame(Q);
}
l <- list()
for(i in 1:N){
  l[[i]] <- dq[[i]]
}
Q<-do.call(rbind,l)
Q<-t(Q)

##build F
df<-vector()
for (s in 1:M) {
  F<-runif(K,0,1);
  df[s]<-data.frame(F);
}
l <- list()
for(s in 1:M){
  l[[s]] <- df[[s]]
}
F<-as.matrix(do.call(rbind,l))

##build PI
pi<-F%*%Q
M<-dim(pi)[1]
N<-dim(pi)[2]
#build G
g<-vector()
dg<-vector()
for (s in 1:M) {
  for (i in 1:N) {
    g[i]<-rbinom(1,2,pi[s,i]);
  }
  dg[s]<-data.frame(g);
}
l <- list()
for(s in 1:M){
  l[[s]] <- dg[[s]]
}
G<-do.call(rbind,l)
dim(G)
##save data
write.csv(G, file = "/Users/xx/Desktop/simulated data/G(Scenario x).csv")
##load data
G<-as.matrix(read.csv("/Users/xx/Desktop/simulated data/G(Scenario x).csv", stringsAsFactors = F)[,-1])

#methods

##Using PCA 1
###D_hat_M
delta<-0
delta_hat<-0
for(i in 1:N){
  for(s in 1:M){
    delta[s]<-G[s,i]*(2-G[s,i])
  }
  delta_hat[i]<-mean(delta)
}
###or
delta_hat<-colMeans(G*(2-G))

D_hat_M<-diag(delta_hat)
###R_hat_M & eigenvectors
M<-dim(G)[1]
N<-dim(G)[2]

R_hat_M<-t(G)%*%G/M-D_hat_M
r_hat<-eigen(R_hat_M)
r_hat$values

###Q_hat_M and P_hat
k_prime<-3
Q_hat_M<-t(r_hat$vec[,1:k_prime])
P_hat<-t(Q_hat_M)%*%Q_hat_M

##Using an estimate of Q_{k'}
###scenario_2.1b.2.Q is estimated from software ADMIXTURE
q <- as.matrix(read.table("/Users/xx/Desktop/scenario_2.1b.2.Q"))
Q_hat_K<-t(q)
library('pracma')
P_hat<-t(Q_hat_K)%*%inv(Q_hat_K%*%t(Q_hat_K))%*%Q_hat_K

##Using an estimate of Pi_{k'}
pi_hat<-1/2*G%*%t(Q_hat_M)%*%Q_hat_M
P_hat<-t(Q_hat_M)%*%Q_hat_M

##Using PCA 2
avg<-rowMeans(G)
g1<-G-avg
dim(g1)
### decompose G=UDV 
udv <- svd(g1)
###evaluate
k <- 2
S<-cbind(udv$v[,1:k-1],1)
S_K<-t(S)
P_hat<-t(S_K)%*%inv(S_K%*%t(S_K))%*%S_K

##Using PCA 3
dim(g1)
sd <- apply(g1,1,sd)
g2 <- g1/sd
dim(g2)
### remove sd=0
g2 <- g2[sd>0,] 
### decompose G=UDV 
udv_gen <-svd(g2)
### use top k-1 PCA 
k <- 3
S<-cbind(udv_gen$v[,1:k-1],1)
S_K<-t(S)
P_hat<-t(S_K)%*%inv(S_K%*%t(S_K))%*%S_K

##Calculate the residuals
N<-dim(P_hat)[1]
R_hat<-G%*%(diag(N)-P_hat)
dim(R_hat)
##hat(b)
cor_emp<-cor(R_hat)
##hat(c)
###cov
delta_hat<-colMeans(G*(2-G))
D_hat_M<-diag(delta_hat)

cov_est<-2*(diag(N)-P_hat)%*%D_hat_M%*%(diag(N)-P_hat)
###cor
corr<-vector()
dcorr<-vector()
for (i in 1:N) {
  for (j in 1:N) {
    corr[j]<-cov_est[i,j]/sqrt(cov_est[i,i]*cov_est[j,j])
  }
  dcorr[i]<-data.frame(corr);
}
l <- list()
for(s in 1:N){
  l[[s]] <- dcorr[[s]]
}
cor_est<-do.call(rbind,l)

##hat(b)-hat(c)
diff_0<-cor_emp-cor_est
diag(diff_0)<-NA
##The mean(sd) of hat(b) and hat(b)-hat(c)
diag(cor_emp)<-NA
mean_sd(cor_emp, pop = pop, is.ord=T)
mean_sd(diff_0, pop = pop, is.ord=T)
##Preparing for Visualization
###The upper triangle is hat(b) and the lower triangle is hat(b)-hat(c)
for(i2 in 1:(N-1)){
  for(i1 in (i2+1):N){
    diff_0[i1, i2] <- cor_emp[i1, i2]
    
  }
}
write.csv(diff_0, file = "/Users/xx/Desktop/corrs.csv")
###define pop label for scenario 1
pop<-c(rep("pop1",10),rep("pop2",20),rep("pop3",30))
pop<-c(rep("pop1",20),rep("pop2",20),rep("pop3",20))
###define pop label for scenario 2
pop <- c(rep("pop1",20), rep("admixed",20), rep("pop3", 20))
pop <- c(rep("pop1",10), rep("admixed",20), rep("pop3", 30))
##load data
corResAdmix <- as.matrix(read.csv("/Users/xx/Desktop/corrs.csv", stringsAsFactors = F)[,-1])
##Visualization
palette("default")
pdf("/Users/xx/Desktop/corrs.pdf",
    width =  8, height = 7.5, paper = "special")
plotCorRes(corResAdmix, pop=pop, is.ord=T, title="", plot_legend = T, pop_labels = c(T,T),
           min_z=-0.2, max_z=0.2, cex.main=3, cex.lab=2,cex.legend=2)
dev.off() 

#theoretical correlation
P_K<-t(Q)%*%inv(Q%*%t(Q))%*%Q
mu<-matrix(0.5,1,3)
sigma<-matrix(c(1/12,0,0,0,1/12,0,0,0,1/12),3,3)
mu_hat<-mu%*%Q-(mu%*%Q)^2
sigma_hat<-t(Q)%*%sigma%*%Q
delta<-mu_hat-diag(sigma_hat)
D<-diag(delta[1,],length(delta),length(delta))
cov<-2*(diag(length(delta))-P_K)%*%D%*%(diag(length(delta))-P_K)
##cor
corr<-vector()
dcorr<-vector()
for (i in 1:length(delta)) {
  for (j in 1:length(delta)) {
    corr[j]<-cov[i,j]/sqrt(cov[i,i]*cov[j,j])
  }
  dcorr[i]<-data.frame(corr);
}
l <- list()
for(s in 1:N){
  l[[s]] <- dcorr[[s]]
}
cor<-do.call(rbind,l)
diag(cor)<-NA
mean_sd(cor, pop = pop, is.ord=T)

#plink data
##simulated data
source("http://popgen.dk/albrecht/open/online.R")
pl1<-plink("/Users/xx/Desktop/plink data/ghost1")
###load data
G <- t(pl1$geno)
dim(G)
###read pop label
pop <- pl1$fam$V2

##real data
pl1<-plink("/Users/xx/Desktop/1000G")
G <- t(pl1$geno)
dim(G)
G<-na.omit(G)
###read pop label
pop <- pl1$fam$V2
pop_order <- c('YRI', 'ASW', 'CEU', 'MXL', 'CHB')
###order by pop
ord1 <- order(match(pop, pop_order))
pop <- pop[ord1]
### Q and corrs
### Q estimated from software ADMIXTURE by setting k=4
qk4 <- as.matrix(read.table("/Users/xx/Desktop/1000G/1000G5popsHumanOriginsSitesHg19Maf005.4.Q"))[ord1,]

corresk4 <- as.matrix(read.csv("/Users/xx/Desktop/1000G/corrreal5p(c&s,4+).csv", stringsAsFactors = F)[,-1])[ord1,ord1]
###order within pop by min to max main admixture 
for(p in pop_order){
  idx <- which(pop==p)
  main_k <- which.max(apply(qk4[idx,],2,mean))
  neword <- order(qk4[idx,main_k])
  pop[idx] <- pop[idx][neword]
  qk4[idx,] <- qk4[idx,][neword,]
  corresk4[idx,] <- corresk4[idx,][neword,]
  corresk4[,idx] <- corresk4[,idx][,neword]
}
### order ancestral populations
ordk <- function(q, k, refpops = c("YRI", "CEU", "CHB","MXL"), pop){
  
  refpops <- refpops[1:k]
  
  kord <- integer(0)
  for(p in refpops){
    
    kord <- c(kord, which.max(apply(q[pop==p,],2,mean)))
    
  }
  
  kord
}
qk4 <- qk4[, ordk(qk4,4,pop=pop)]
###plot qk4
pdf("/Users/xx/Desktop/adreal+.pdf",
    width =  8.5, height = 7.5, paper = "special")
par(mar=c(5, 4, 4, 2), mgp=c(2.5, 1, 0))
barplot(t(qk4), border=NA, space=0, ylab="admixture proportion", cex.axis=1.5, cex.lab=2, col=c(1,2,3,4),font.lab=2,font.axis=2)
mtext(at=c(54.5, 139.0, 219.0, 300.5, 384.0),side=1,
      text=unique(pop),xpd=T, padj=0.8, font = 2, cex=2)
dev.off()

