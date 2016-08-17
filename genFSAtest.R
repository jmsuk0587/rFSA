dat <- read.delim(file.choose())
vec<-rep(NA,dim(dat)[2]-1)
for(i in 1:(dim(dat)[2]-1)){
  vec[i]<-length(levels(dat[,i]))
}

## Remove SNPs that only have one genotype
if (length(which(vec<=1))>0){ print('yes')
  keeps<-(1:(dim(dat)[2]-1))[-which(vec<=1)]
  dat<-dat[,c(keeps,dim(dat)[2])]
}

library(devtools)
install_github("jmsuk0587/rFSA")
library(rFSA)
fit<-genFSA(yname = "phen",data = dat,quad = F,m=2,cores = 2,numrs = 10,checknum = 100)
