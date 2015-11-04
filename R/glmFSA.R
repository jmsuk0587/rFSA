#' FSA Generalized Linear Models
#' @description  A function using a Feasible Soultion Algorithm to find a set of feasible (optimal in that no one swap to any of the variables will increase the criterion function) mth order interaction Generalized Linear Model
#' @param formula an object of class "formula" (or one that can be coerced to that class): a symbolic description of the model to be fitted. The details of model specification are given under "Details".
#' @param data a data frame, list or environment (or object coercible by as.data.frame to a data frame) containing the variables in the model. If not found in data, the variables are taken from environment(formula), typically the environment from which lm is called.
#' @param fixvar=NULL a variable to fix in the model. Usually a covariate that should always be included (Age, Sex, ...). Will still consider it with interactions.
#' @param quad to include quadratic terms or not.
#' @param m=2 order of interaction to look at. Defults to 2.
#' @param numrs number of random starts to perform.
#' @param save_results whether to save the results in the current working directory as 'FSAresults.csv'.
#' @param cores number of cores to use while running. Note: Windows can only use 1 core. See mclapply for details.
#' @param fam family argument passed to glm. a description of the error distribution and link function to be used in the model. This can be a character string naming a family function, a family function or the result of a call to a family function.
#' @param ... arguments to be passed to the glm function
#'
#' @details PLEASE NOTE: make sure categorical variables are factors or characters otherwise answers will not reflect the variable being treated as a continuous variable.
#' @return a list of results and a table of unique results.
#' @export
#'
#' @examples
#' #use mtcars package see help(mtcars)
#' data(mtcars)
#' colnames(mtcars)
#' df<-mtcars
#' df[,c(2,8:11)]<-lapply(mtcars[,c(2,8:11)],factor)
#' glmFSA(am~cyl*disp,data=df,fixvar="hp",quad=F,m=2,numrs=10,save_solutions=F,fam="binomial",cores=1)
#' fit<-glm(am~hp*qsec,data=df,family="binomial") #this is the most common answer from lmFSA.
#' summary(fit) #review
glmFSA=function(formula,data,fixvar=NULL,quad=F,m=2,numrs=1,save_solutions=T,fam,cores=detectCores()/3,...){
  originalnames<-colnames(data)
  data<-data.frame(data)
  lhsvar<-lhs(formula)
  startvar<-get.vars()
  fit=glm(formula,data=data,family=fam,...)
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
      tmp<-mclapply(X = 1:dim(moves)[2],FUN = function(k){
        fit<-glm(form(k),data=newdata,family=fam,...)
        anv<-anova(fit,test="Chisq")
        c(summary(fit)$aic,anv[dim(anv)[1],dim(anv)[2]])
      },mc.cores=cores)

      tmp<-matrix(unlist(tmp),byrow=T,ncol=2)
      cur<-moves[,which.min(tmp[,2])]
      r.sq<-tmp[which.min(tmp[,2]),2]
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
  colnames(solutions)[dim(solutions)[2]:(dim(solutions)[2]-1)]=c("swaps","p-value")
  return(list(solutions=solutions,table=tableres))
}
