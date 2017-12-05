
# This script simulates data for the vignette.

ls <- ls()
set.seed(1)

n <- 100
q <- 5
p <- 5
mean <- rep(0,times=n)
mean[sample(51:100,size=10)] <- 5
y <- stats::rnorm(n=n,mean=mean)
Y <- matrix(stats::rnorm(n*q),nrow=n,ncol=q)
snp <- rep(c(0,1,2),times=c(50,40,10))
z <- rep(NA,times=n)
z[snp==0] <- 0
SNPs <- matrix(sample(c(0,1,2),replace=TRUE,size=n*p),nrow=n,ncol=p)
Z <- matrix(NA,nrow=n,ncol=p)
Z[SNPs==0] <- 0

toydata <- list(Y=Y,SNPs=SNPs,Z=Z,y=y,snp=snp,z=z)

rm(list=setdiff(ls(),c(ls,"toydata")))