
newalg<-function(yname,dat,parmax,intmax,nrs,fitfunc,...){
  res2<-NULL
  res2<-list(1)
  for(u in 2:parmax){
    res2[[u]]<-fitfunc(yname=yname,data=dat,m=u,interactions=F,numrs=nrs,cores=4,...)$table
    res2[[u]]<-res2[[u]][order(res2[[u]][,1],decreasing = T),]
    }
  
  res3<-rep("--",parmax+2*intmax+3)
  for(pnum in 2:parmax){
      for(ints in 2:intmax){
        for(z in 1:min(c(2,dim(res2[[pnum]])[1]))){
          tmps<-NULL
          show(c(pnum,ints,z))
          tmps<-fitfunc(yname=yname,data=dat,fixvar=res2[[pnum]][z,2:(pnum+1)],m=ints,interactions=T,numrs=nrs,cores=4,...)$solutions
          s<-unlist(lapply(X = 1:dim(tmps)[2],FUN = function(ii){paste(as.character(tmps[order(tmps$criterion,decreasing = T),][1,ii]))}))
          s[(pnum+1):(pnum+ints)]<-"--"
          if((pnum+2*ints+3+1)<=(parmax+2*intmax+3)){s[(pnum+2*ints+3+1):(parmax+2*intmax+3)]<-"--"}
          res3<-rbind(res3,s)
      }
    }
  }
  
  return(list(results=res3))
}

set.seed(12345)
dat<-matrix(rnorm(100*100,sd=1/10),ncol=100)
dat<-data.frame(dat)

d<-newalg(yname = "X1",dat=dat,parmax=5,intmax=3,nrs=3,fitfunc=lmFSA,criterion=adj.r.squared,minmax="max")
d<-newalg(yname = "X1",dat=dat,parmax=5,intmax=3,nrs=3,fitfunc=lmFSA,criterion=BIC,minmax="min")
