bFSA=function(yname,data,numrs=1,cores=1){
  func<-function(x1,x2,y){
  X <- cbind(x1,x2)
  a <- X[which(y==0),]
  b <- X[which(y==1),]
  
  mu11 <- mean(a[,1])
  mu12 <- mean(a[,2])
  mu21 <- mean(b[,1])
  mu22 <- mean(b[,2])
  
  var1_11 <- var(a[,1])
  var1_22 <- var(a[,2])
  var1_12 <- cov(a[,1],a[,2])
  
  var2_11 <- var(b[,1])
  var2_22 <- var(b[,2])
  var2_12 <- cov(b[,1],b[,2])
  
  mu1 <- c(mu11,mu12)
  mu2 <- c(mu21,mu22)
  sig1 <- matrix(c(var1_11,var1_12,var1_12,var1_22),ncol=2)
  sig2 <- matrix(c(var2_11,var2_12,var2_12,var2_22),ncol=2)
  
  sig <- (sig1+sig2)/2
  a <- sig[1,1]
  b <- sig[1,2]
  c <- sig[2,1]
  d <- sig[2,2]
  temp <- matrix(c(d,-c,-b,a),ncol=2)
  sig.inv <- 1/(a*d-b*c)*temp
  distance <- 1/8*(t(mu1-mu2)%*%sig.inv%*%(mu1-mu2))+1/2*log(det(sig)/sqrt(det(sig1)*det(sig2)))
  return(distance) 
  }
  
  fsa<-function(i,history,...){
    cur<-history[i,1:2]
    last<-rep(NA,2)
    numswap<-0
    memswap<-NULL
    checks<-0
    last.criterion<-(-Inf)
    while(!identical(cur,last)){
      last<-cur
      if(numswap==0){moves<<-swaps(cur = cur,n = dim(newdata[,-1])[2],quad=F)}
      if(numswap>0){moves<<-nextswap(curpos = cur,n = dim(newdata[,-1])[2],quad=F,prevpos =memswap)$nswaps}

      tmp<-mclapply(X = 1:dim(moves)[2],FUN = function(k) func(newdata[,-1][,moves[1,k]],newdata[,-1][,moves[2,k]],newdata[,1]) ,mc.cores=cores)
      checks<-checks+dim(moves)[2]
      cur<-moves[,which.max.na(unlist(tmp))[1]]
      cur.criterion<-unlist(tmp[which.max.na(unlist(tmp))[1]])
      if(last.criterion>cur.criterion){cur<-last.pos
      cur.criterion<-last.criterion}
      
      
      
      
      numswap<-numswap+1
      last1<-last
      last.criterion<-cur.criterion
      last.pos<-cur
      memswap<-unique(c(memswap,last1))
    }
    history[i,(1+2):(2*2)]<-cur
    history[i,(dim(history)[2]-2)]<-cur.criterion
    history[i,(dim(history)[2]-1)]<-numswap-1
    history[i,(dim(history)[2])]<-checks
    return(history[i,])
  }
  
  ypos<-which(colnames(data)==yname)
  newdata<-data.frame(cbind(data[,ypos],data[,-ypos]))
  colnames(newdata)[1]<-yname
  
  history<-matrix(rep(NA,numrs*(2*2+3)),ncol=((2*2+3)))
  history[,1:2]<-rstart(m=2,nvars=(dim(newdata)[2]-1),numrs=numrs)
  
  
  solutions<-matrix(unlist(lapply(1:numrs,FUN =function(i) fsa(i,history))),ncol=dim(history)[2],byrow = T)
  solutions[,1:(2*2)]<-matrix(colnames(newdata)[c(solutions[,1:(2*2)]+1)],ncol=(2*2))
  solutions<-data.frame(solutions)
  
  colnames(solutions)[dim(solutions)[2]:(dim(solutions)[2]-2)]=c("checks","swaps","criterion")
  colnames(solutions)[1:2]=paste("start",1:2,sep=".")
  colnames(solutions)[(2+1):(2*2)]=paste("best",1:2,sep=".")
  solutions$criterion<-as.numeric(levels(solutions$criterion))[solutions$criterion]
  solutions$swaps<-as.numeric(levels(solutions$swaps))[solutions$swaps]
  solutions$checks<-as.numeric(levels(solutions$checks))[solutions$checks]
  solutions<<-solutions
  a<-solutions[,(2+1):(2+1+2)]
  b<-unique(t(apply(a,sort,MARGIN = 1)),MARGIN = 1)
  a<-t(apply(a,sort,MARGIN = 1))
  c<-cbind(b,0)
  for(i in 1:dim(b)[1]){
    for(j in 1:dim(a)[1]){
      c[i,(2+2)]<-sum(as.numeric(c[i,(2+2)])+as.numeric(identical(a[j,],b[i,])))
    }
  }
  tableres<-data.frame(cbind(c,NA),stringsAsFactors = F)
  colnames(tableres)[(dim(tableres)[2]-1)]<-"times"
  colnames(tableres)[2:(dim(tableres)[2]-2)]<-paste("Var",2:(dim(tableres)[2]-2),sep="")
  colnames(tableres)[1] <- "criterion"
  colnames(tableres)[dim(tableres)[2]]<-"warnings"
  tableres$warnings<-"NA"
  return(list(solutions=solutions,table=tableres,efficiency=paste("You did:",sum(solutions$checks)," model checks compared to ",choose(n = (dim(newdata)[2]-1),k = 2)," checks you would have done with exahstive search.")))
}