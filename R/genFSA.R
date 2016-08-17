#' Feasible Solution Algorithm (FSA) for Genetic Models
#' @description  A function using a Feasible Solution Algorithm to find a set of feasible solutions for a genetic model of a specific form that could include mth-order interactions (Note that these solutions are optimal in the sense that no one swap to any of the variables will increase the criterion function.)
#' @param yname the character form of the name of the y (response) variable. Please write the name of the y variable in quotes (example: "Disease").
#' @param data a data frame, list or environment (or object coercible by as.data.frame to a data frame) containing the variables in the model. 
#' @param fixvar=NULL a variable to fix in the model. Usually a covariate that should always be included (Example: Age, Sex). Will still consider it with interactions.
#' @param quad to include quadratic terms or not.
#' @param m=2 order of terms to potentially include. If interactions is set to TRUE then m is the order of interactions to be considered. Defaults to 2. For Subset selection (interaction=F), m is the size of the subset to examine.
#' @param numrs number of random starts to perform.
#' @param save_solutions whether to save the solutions in the current working directory as 'FSAsolutions.csv'.
#' @param cores number of cores to use while running. Note: Windows can only use 1 core. See mclapply for details.
#' @param fam family argument passed to glm. A description of the error distribution and link function to be used in the model. This can be a character string naming a family function, a family function or the result of a call to a family function.
#' @param interactions T or F for whether to include interactions in model. Defaults to FALSE. 
#' @param criterion which criterion function to either maximize or minimize. For linear models one can use: apress (Allens Press Statistic), int.p.val (Interaction p-value), AIC, BIC.
#' @param minmax whether to minimize or maximize the criterion function
#' @param checknum Number of models to fit at each swap. 
#' @param ... arguments to be passed to the glm function
#' @details PLEASE NOTE: make sure categorical variables are factors or characters otherwise answers will not reflect the variable being treated as a continuous variable.
#' @return returns a list of solutions and table of unique solutions.
#' $solutions is a matrix of fixed terms, start position, feasible solution, criterion function value (p-value of interaction), and number of swaps to solution.
#' $table is a matrix of the unique feasible solutions and how many times they occurred out of the number of random starts chosen. It also returns any warning messages with these solutions in the last column.
#' $efficiency is text comparing how many models you ran during your FSA search compared to how many you would have done with exhaustive search. Note: The FSA algorithm takes additional time to run on top of the model checks that were done during the algorithm. This additional time is approximately 15% more time than if you had just ran the model checks. 
#' @export

genFSA=function(yname,data,fixvar=NULL,quad=F,m=2,numrs=1,save_solutions=F,cores=1,interactions=F,criterion=AIC,minmax="min",fam="binomial",checknum=NA,...){

  originalnames<-colnames(data)
  data<-data.frame(data)
  lhsvar<-yname
  
  ypos<-which(colnames(data)==lhsvar)
  startvar<-NULL
  xdata<-data[,-ypos]
  ydata<-data[,ypos]

  newdata<<-data.frame(cbind(ydata,xdata))
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
    while(!identical(cur,last)&&!identical(c(cur[2],cur[1]),last)){
      last<-cur
      if(numswap==0){moves<<-swaps(cur = cur,n = dim(xdata)[2],quad=quad)}
      if(numswap>0){moves<<-nextswap(curpos = cur,n = dim(xdata)[2],quad=quad,prevpos =memswap)$nswaps
      }
      if(m == 2){
      vec<-rep(NA,dim(moves)[2])
      for(g in 1:(dim(moves)[2])){
        if(is.factor(xdata[,moves[1,g]]) && is.factor(xdata[,moves[2,g]])){
          int<-(paste(xdata[,moves[1,g]],xdata[,moves[2,g]]))
          int[is.na(xdata[,moves[1,g]])|is.na(xdata[,moves[2,g]])]<-NA
          int<-as.factor(int)
          h <- length(levels(int))
          int1 <- xdata[,moves[1,g]]
          int2 <- xdata[,moves[2,g]]
          l1 <- length(levels(int1))
          l2 <- length(levels(int2))
          if(l1*l2 != h){
            vec[g] <- 1
          }
        }
      }
      if(length(which(vec==1))>0){
        keeps<-(1:(dim(moves)[2]))[-which(vec==1)]
        moves <<- moves[,keeps]
      }
      }
      if(is.null(dim(moves))||dim(moves)[2]==0){cur=last;cur.criterion=last.criterion;numswap=numswap+1;break}
      
      if((!is.na(checknum)) && (checknum<dim(moves)[2])){moves<<-moves[,sample(1:dim(moves)[2],size=checknum,replace=F)]}
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
  colnames(solutions)[dim(solutions)[2]:(dim(solutions)[2]-2)]=c("checks","swapsn","criterion")
  colnames(solutions)[1:m]=paste("start",1:m,sep=".")
  colnames(solutions)[(m+1):(m*2)]=paste("best",1:m,sep=".")
  solutions$criterion<-as.numeric(levels(solutions$criterion))[solutions$criterion]
  solutions$swapsn<-as.numeric(levels(solutions$swapsn))[solutions$swapsn]
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
