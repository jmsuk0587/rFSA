#' FSA Linear Models
#' @description  A function using a Feasible Soultion Algorithm to find a set of feasible solutions for the linear model including mth-order interactions (Note that these solutions are optimal in the sense that no one swap to any of the variables will increase the criterion function.)
#' @param formula an object of class "formula" (or one that can be coerced to that class): a symbolic description of the model to be fitted. The details of model specification are given under "Details".
#' @param data a data frame, list or environment (or object coercible by as.data.frame to a data frame) containing the variables in the model. If not found in data, the variables are taken from environment(formula), typically the environment from which lm is called.
#' @param fixvar=NULL a variable to fix in the model. Usually a covariate that should always be included (Age, Sex, ...). Will still consider it with interactions.
#' @param quad to include quadratic terms or not.
#' @param m=2 order of interaction to look at. Defults to 2.
#' @param numrs number of random starts to perform.
#' @param save_solutions whether to save the solutions in the current working directory as 'FSAsolutions.csv'.
#' @param cores number of cores to use while running. Note: Windows can only use 1 core. See mclapply for details.
#' @param ... arguments to be passed to the lm function
#'
#' @details PLEASE NOTE: make sure categorical variables are factors or characters otherwise answers will not reflect the variable being treated as a continuous variable.
#' @return returns a list of solutions and table of unique solutions.
#' $solutions is a matrix of fixed terms, start position, feasible solution, criterion function value, and number of swaps to solution.
#' $table is a matrix of the unique feasible solutions and how many times they occured out of the number of random starts chosen.
#' @export
#'
#' @examples
#' #use mtcars package see help(mtcars)
#' data(mtcars)
#' colnames(mtcars)
#' lmFSA(mpg~cyl*disp,data=mtcars,fixvar="hp",quad=F,m=2,numrs=10,save_solutions=F,cores=1)
#' fit<-lm(mpg~hp*wt,data=mtcars) #this is the most common answer from lmFSA.
#' summary(fit) #review
lmFSA=function(formula,data,fixvar=NULL,quad=F,m=2,numrs=1,save_solutions=T,cores=1,...){
  originalnames<-colnames(data)
  data<-data.frame(data)
  lhsvar<-lhs(formula)
  startvar<-get.vars()
  fit=lm(formula,data=data,...)
  ypos<-which(colnames(data)==lhsvar)

  xdata<-data[,-ypos]
  ydata<-data[,ypos]
  newdata<-data.frame(cbind(ydata,xdata))
  fixpos<-which(colnames(xdata) %in% fixvar)
  if(length(fixpos)==0){fixpos=NULL}

  history<-matrix(rep(NA,numrs*(2*m+2)),ncol=(2*m+2))
  history[,1:m]<-rstart(m=m,nvars=(dim(newdata)[2]-1),numrs=numrs)
  curpos<-which(colnames(xdata) %in% startvar[-1])
  if(length(curpos)!=0){history<-rbind(c(curpos,rep(NA,length(curpos)+2)),history)}

  fsa<-function(i,history,...){
    cur<-history[i,1:m]
    last<-rep(NA,m)
    numswap<-0
    while(!identical(cur,last)){
      last<-cur
      moves<-swaps(cur = cur,n = dim(xdata)[2],quad=quad)
      form<-function(j) formula(paste0(colnames(newdata)[1],"~",paste0(fixvar,sep="+"),paste(colnames(xdata)[moves[,j]],collapse = "*")),sep="")
      tmp<-mclapply(X = 1:dim(moves)[2],FUN = function(k) summary(lm(form(k),data=newdata,...))$r.squared,mc.cores=cores)
      cur<-moves[,which.max(tmp)]
      r.sq<-unlist(tmp[which.max(tmp)])
      numswap<-numswap+1
    }
    history[i,(1+m):(2*m)]<-cur
    history[i,(dim(history)[2]-1)]<-r.sq
    history[i,(dim(history)[2])]<-numswap-1
    return(history[i,])
  }
  solutions<-matrix(unlist(lapply(1:numrs,FUN =function(i) fsa(i,history))),ncol=dim(history)[2],byrow = T)
  solutions[,1:(2*m)]<-matrix(colnames(newdata)[c(solutions[,1:(2*m)]+1)],ncol=(2*m))
  solutions<-data.frame(solutions)
  colnames(solutions)[dim(solutions)[2]:(dim(solutions)[2]-1)]=c("swaps","r.sq")
  colnames(solutions)[1:m]=paste("start",1:m,sep=".")
  colnames(solutions)[(m+1):(m*2)]=paste("best",1:m,sep=".")
  solutions$r.sq<-as.numeric(levels(solutions$r.sq))[solutions$r.sq]
  solutions$swaps<-as.numeric(levels(solutions$swaps))[solutions$swaps]
  if(length(fixvar)!=0){solutions<-data.frame(fixvar=matrix(rep(x=fixvar,dim(solutions)[1]),nrow=dim(solutions)[1],byrow=T),solutions)}
  if(save_solutions==T){write.csv(solutions,paste0(getwd(),"/FSAsolutions",".csv"))}
  solutions<<-solutions
  a<-solutions[,(length(fixvar)+m+1):(length(fixvar)+m+1+m-1)]
  b<-unique(t(apply(a,sort,MARGIN = 1)),MARGIN = 1)
  a<-t(apply(a,sort,MARGIN = 1))
  c<-cbind(b,0)
  for(i in 1:dim(b)[1]){
    for(j in 1:dim(a)[1]){
      c[i,(m+1)]<-sum(as.numeric(c[i,(m+1)])+as.numeric(identical(a[j,],b[i,])))
    }
  }
  tableres<-c
  colnames(tableres)[dim(c)[2]]<-"times"
  return(list(solutions=solutions,table=tableres))
}
