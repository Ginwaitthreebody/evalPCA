# evalPCA
Scripts used for the paper "Evaluation of model fit" 2023

If you already have the data, you can skip to the methods section.

##build Q,F,PI,G for Scenario 1 and 2

##build Q1 for Scenario 1(specific for N_1,..,N_{K+1})

##build Q2 for Scenario 2(specific for N_1,..,N_{K+1})

##build F

##build PI

methods

For PCA1, the main thing is to compute P_hat.
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

So the residual part of the calculation we started with the covariance.
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

We also have some methods for naming the population tag to which the sample individual belongs. For example,
###define pop label for scenario 1
pop<-c(rep("pop1",10),rep("pop2",20),rep("pop3",30))
pop<-c(rep("pop1",20),rep("pop2",20),rep("pop3",20))
###define pop label for scenario 2
pop <- c(rep("pop1",20), rep("admixed",20), rep("pop3", 20))
pop <- c(rep("pop1",10), rep("admixed",20), rep("pop3", 30))

The end result is our visual diagram, which can be combined and compared if needed.
##load data
corResAdmix <- as.matrix(read.csv("/Users/xx/Desktop/corrs.csv", stringsAsFactors = F)[,-1])
##Visualization
palette("default")
pdf("/Users/xx/Desktop/corrs.pdf",
    width =  8, height = 7.5, paper = "special")
plotCorRes(corResAdmix, pop=pop, is.ord=T, title="", plot_legend = T, pop_labels = c(T,T),
           min_z=-0.2, max_z=0.2, cex.main=3, cex.lab=2,cex.legend=2)
dev.off() 
