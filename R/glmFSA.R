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
#' @param interactions T or F for whether to include interactions in model. 
#' @param ... arguments to be passed to the glm function
#' @details PLEASE NOTE: make sure categorical variables are factors or characters otherwise answers will not reflect the variable being treated as a continuous variable.
#' @return returns a list of solutions and table of unique solutions.
#' $solutions is a matrix of fixed terms, start position, feasible solution, criterion function value (p-value of interaction), and number of swaps to solution.
#' $table is a matrix of the unique feasible solutions and how many times they occured out of the number of random starts chosen. It also returns any warning messages with these solutions in the last column.
#' @export
#' @examples
#' dat<-read.csv("https://archive.ics.uci.edu/ml/machine-learning-databases/hepatitis/hepatitis.data",header = F)
#' colnames(dat)<-c("Class","Age","Sex","Sterioid","Antivirals","Fatigue","Malaise","Anorexia","Liver Big","Liver Firm","Spleen Palpable","Spiders","Ascites","Varices","Bilirubin","Alk Phosphate","Sgot","Albumin","Protime","Histology")
#' dat<-as.matrix(dat)
#' dat[which(dat=="?")]=NA
#' dat<-data.frame(dat)
#' dat[,c(2,15,16,17,18,19)]<-lapply(X = dat[,c(2,15,16,17,18,19)],as.numeric)
#' colnames(dat)
#' glmFSA(yname="Class",data=dat,fixvar="Age",quad=F,m=2,numrs=10,save_solutions = F,fam="binomial",cores=1)

glmFSA=function(yname,data,fixvar=NULL,quad=F,m=2,numrs=1,save_solutions=F,cores=1,interactions=T,criterion=AIC,minmax="max",fam="binomial",...){
  originalnames<-colnames(data)
  data<-data.frame(data)
  lhsvar<-yname

  ypos<-which(colnames(data)==lhsvar)
  
  xdata<-data[,-ypos]
  ydata<-data[,ypos]
  newdata<-data.frame(cbind(ydata,xdata))
  fixpos<-which(colnames(xdata) %in% fixvar)
  if(length(fixpos)==0){fixpos=NULL}
  
  history<-matrix(rep(NA,numrs*(2*m+3)),ncol=((2*m+3)))
  history[,1:m]<-rstart(m=m,nvars=(dim(newdata)[2]-1),numrs=numrs)
  curpos<-which(colnames(xdata) %in% startvar[-1])
  if(length(curpos)!=0){history<-rbind(c(curpos,rep(NA,length(curpos)+2)),history)}
  
  fsa<-function(i,history,...){
    cur<-history[i,1:m]
    last<-rep(NA,m)
    numswap<-0
    memswap<-NULL
    if(minmax=="max"){last.criterion<-(-Inf)}
    if(minmax=="min"){last.criterion<-(Inf)}
    checks<-0
    while(!identical(cur,last)){
      last<-cur
      if(numswap==0){moves<<-swaps(cur = cur,n = dim(xdata)[2],quad=quad)}
      if(numswap>0){moves<<-nextswap(curpos = cur,n = dim(xdata)[2],quad=quad,prevpos =memswap)$nswaps
      }
      if(interactions==T){form<-function(j) formula(paste0(colnames(newdata)[1],"~",paste0(fixvar,sep="+"),paste(colnames(xdata)[moves[,j]],collapse = "*")),sep="")}
      if(interactions==F){form<-function(j) formula(paste0(colnames(newdata)[1],"~",paste0(fixvar,sep="+"),paste(colnames(xdata)[moves[,j]],collapse = "+")),sep="")}
      tmp<-mclapply(X = 1:dim(moves)[2],FUN = function(k) criterion(glm(form(k),data=newdata,family=fam,...)),mc.cores=cores)
      checks<-checks+dim(moves)[2]
      if(minmax=="max"){cur<-moves[,which.max.na(unlist(tmp))[1]]
      cur.criterion<-unlist(tmp[which.max.na(unlist(tmp))[1]])
      if(last.criterion>cur.criterion){cur<-last.pos
      cur.criterion<-last.criterion}
      }
      if(minmax=="min"){cur<-moves[,which.min.na(unlist(tmp))[1]]
      cur.criterion<-unlist(tmp[which.min.na(unlist(tmp))[1]])
      if(last.criterion<cur.criterion){cur<-last.pos
      cur.criterion<-last.criterion}
      }
      
      
      numswap<-numswap+1
      last1<-last
      last.criterion<-cur.criterion
      last.pos<-cur
      memswap<-unique(c(memswap,last1))
    }
    history[i,(1+m):(2*m)]<-cur
    history[i,(dim(history)[2]-2)]<-cur.criterion
    history[i,(dim(history)[2]-1)]<-numswap-1
    history[i,(dim(history)[2])]<-checks
    return(history[i,])
  }
  solutions<-matrix(unlist(lapply(1:numrs,FUN =function(i) fsa(i,history))),ncol=dim(history)[2],byrow = T)
  solutions[,1:(2*m)]<-matrix(colnames(newdata)[c(solutions[,1:(2*m)]+1)],ncol=(2*m))
  solutions<-data.frame(solutions)
  print
  colnames(solutions)[dim(solutions)[2]:(dim(solutions)[2]-2)]=c("checks","swaps","criterion")
  colnames(solutions)[1:m]=paste("start",1:m,sep=".")
  colnames(solutions)[(m+1):(m*2)]=paste("best",1:m,sep=".")
  solutions$criterion<-as.numeric(levels(solutions$criterion))[solutions$criterion]
  solutions$swaps<-as.numeric(levels(solutions$swaps))[solutions$swaps]
  solutions$checks<-as.numeric(levels(solutions$checks))[solutions$checks]
  if(length(fixvar)!=0){solutions<-data.frame(fixvar=matrix(rep(x=fixvar,dim(solutions)[1]),nrow=dim(solutions)[1],byrow=T),solutions)}
  if(save_solutions==T){write.csv(solutions,paste0(getwd(),"/FSAsolutions",".csv"))}
  solutions<<-solutions
  a<-solutions[,(length(fixvar)+m+1):(length(fixvar)+m+1+m)]
  b<-unique(t(apply(a,sort,MARGIN = 1)),MARGIN = 1)
  a<-t(apply(a,sort,MARGIN = 1))
  c<-cbind(b,0)
  for(i in 1:dim(b)[1]){
    for(j in 1:dim(a)[1]){
      c[i,(m+2)]<-sum(as.numeric(c[i,(m+2)])+as.numeric(identical(a[j,],b[i,])))
    }
  }
  tableres<-data.frame(cbind(c,NA),stringsAsFactors = F)
  colnames(tableres)[(dim(tableres)[2]-1)]<-"times"
  colnames(tableres)[2:(dim(tableres)[2]-2)]<-paste("Var",2:(dim(tableres)[2]-2),sep="")
  colnames(tableres)[1] <- "criterion"
  colnames(tableres)[dim(tableres)[2]]<-"warnings"
  withWarnings <- function(expr) {
    myWarnings <- NULL
    wHandler <- function(w) {
      myWarnings <<- c(myWarnings, list(w))
      invokeRestart("muffleWarning")
    }
    val <- withCallingHandlers(expr, warning = wHandler)
    list(value = val, warnings = myWarnings)
  }
  form<-function(j) formula(paste0(colnames(newdata)[1],"~",paste0(fixvar,sep="+"),paste(tableres[j,2:m+1],collapse = "*")),sep="")
  warns<-NULL
  for (i in 1:dim(tableres)[1]){
    ca<-as.character(withWarnings(glm(form(i),data=newdata,family=fam,...))$warnings[[1]])
    if(length(ca)==0){warns<-c(warns,NA)}
    else{warns<-c(warns,ca)}
  }
  tableres$warnings<-warns
  return(list(solutions=solutions,table=tableres,efficiency=paste("You did:",sum(solutions$checks)," model checks compared to ",choose(n = dim(xdata)[2],k = m)," checks you would have done with exahstive search.")))
}
