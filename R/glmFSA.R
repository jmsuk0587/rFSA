#' rFSA: Feasible Solution Algorithm (FSA) for Generalized Linear Models
#' @description  A function using a Feasible Solution Algorithm to find a set of feasible solutions for a generalized linear model of a specific form that could include mth-order interactions (Note that these solutions are optimal in the sense that no one swap to any of the variables will increase the criterion function.)
#' @param formula an object of class "formula" (or one that can be coerced to that class): a symbolic description of the model to be fitted. See help(glm) for details.
#' @param data a data frame, list or environment (or object coercible by as.data.frame to a data frame) containing the variables in the model. 
#' @param fixvar a variable to fix in the model. Usually a covariate that should always be included (Example: Age, Sex). Will still consider it with interactions. Default is NULL.
#' @param quad to include quadratic terms or not.
#' @param m order of terms to potentially include. If interactions is set to TRUE then m is the order of interactions to be considered. Defaults to 2. For Subset selection (interaction=F), m is the size of the subset to examine. Default is 2. 
#' @param numrs number of random starts to perform.
#' @param cores number of cores to use while running. Note: Windows can only use 1 core. See mclapply for details.  If function detects a Windows user it will automatically set cores=1. 
#' @param family family argument passed to glm. A description of the error distribution and link function to be used in the model. This can be a character string naming a family function, a family function or the result of a call to a family function.
#' @param interactions T or F for whether to include interactions in model. Defaults to FALSE. 
#' @param criterion which criterion function to either maximize or minimize. For linear models one can use: apress (Allens Press Statistic), int.p.val (Interaction p-value), AIC, BIC.
#' @param minmax whether to minimize or maximize the criterion function
#' @param checkfeas A vector of a potential feasible solution
#' @param ... arguments to be passed to the glm function
#' @details PLEASE NOTE: make sure categorical variables are factors or characters otherwise answers will not reflect the variable being treated as a continuous variable.
#' @return returns a list of solutions and table of unique solutions.
#' $solutions is a matrix of fixed terms, start position, feasible solution, criterion function value (p-value of interaction), and number of swaps to solution.
#' $table is a matrix of the unique feasible solutions and how many times they occurred out of the number of random starts chosen. It also returns any warning messages with these solutions in the last column.
#' $efficiency is text comparing how many models you ran during your FSA search compared to how many you would have done with exhaustive search. Note: The FSA algorithm takes additional time to run on top of the model checks that were done during the algorithm. This additional time is approximately 15% more time than if you had just ran the model checks. 
#' @importFrom parallel mclapply
#' @importFrom graphics par plot
#' @importFrom methods show
#' @importFrom stats AIC anova as.formula cov fitted formula glm influence lm predict resid var
#' @importFrom utils capture.output tail
#' @export
#' @examples
#' dat<-read.csv("http://tinyurl.com/zq7l775",header = FALSE)
#' colnames(dat)<-c("Class","Age","Sex","Sterioid","Antivirals",
#'                  "Fatigue","Malaise","Anorexia","Liver Big",
#'                  "Liver Firm","Spleen Palpable","Spiders",
#'                  "Ascites","Varices","Bilirubin","Alk Phosphate",
#'                  "Sgot","Albumin","Protime","Histology")
#' dat<-as.matrix(dat)
#' dat[which(dat=="?")]=NA
#' dat<-data.frame(dat)
#' dat[,c(2,15,16,17,18,19)]<-lapply(X = dat[,c(2,15,16,17,18,19)],as.numeric)
#' colnames(dat)
#' fit<-glmFSA(formula="Class~Age+Sgot*Albumin",data=dat,fixvar="Age",quad=FALSE,m=2,
#'        numrs=10,family="binomial",cores=1)
glmFSA = function(formula,data,fixvar = NULL,quad = FALSE,m = 2,numrs = 1,cores = 1,
                  interactions = TRUE,criterion = AIC,minmax = "min",family =
                    "binomial",checkfeas=NULL,...) {
  if(identical(criterion,r.squared)|identical(criterion,r.squared)){return(show("Sorry the criterion function you listed cannot be used with glmFSA."))}
  formula <- as.formula(formula)
  fit <- glm(formula,data = data,family = family,...)
  yname <- all.vars(formula)
  if (!all(c(yname,fixvar) %in% colnames(data))) {
    return(
      show(
        "Sorry, one of the variables you specified in your formula or fixvar is not a name for a column in the data you specified. Please try again."
      )
    )
  }
  originalnames <- colnames(data)
  data <- data.frame(data)
  lhsvar <- yname[1]
  
  if (.Platform$OS.type == "unix") {
  } else {
    cores = 1
  }
  
  
  ypos <- which(colnames(data) == lhsvar)
  startvar <- NULL
  xdata <- data[,-ypos]
  ydata <- data[,ypos]
  newdata <- data.frame(cbind(ydata,xdata))
  fixpos <- which(colnames(xdata) %in% fixvar)
  if (length(fixpos) == 0) {
    fixpos = NULL
  }
  
  history <- matrix(rep(NA,numrs * (2 * m + 3)),ncol = ((2 * m + 3)))
  
  if(!is.null(checkfeas)){
    checkfeas<-which(colnames(xdata) %in% checkfeas)
    history[,1:m] <-rbind(rstart(m = m,nvars = (dim(newdata)[2] - 1),numrs = numrs-1),c(checkfeas[1:m]))
  } else  history[,1:m] <- rstart(m = m,nvars = (dim(newdata)[2] - 1),numrs = numrs)
  
  
  curpos <- which(colnames(xdata) %in% startvar[-1])
  if (length(curpos) != 0) {
    history <- rbind(c(curpos,rep(NA,length(curpos) + 2)),history)
  }
  
  fsa <- function(i,history,...) {
    cur <- history[i,1:m]
    last <- rep(NA,m)
    numswap <- 0
    memswap <- NULL
    if (minmax == "max") {
      last.criterion <- (-Inf)
    }
    if (minmax == "min") {
      last.criterion <- (Inf)
    }
    checks <- 0
    while (!identical(cur,last) && !identical(c(cur[2],cur[1]),last)) {
      last <- cur
      if (numswap == 0) {
        moves <- swaps(cur = cur,n = dim(xdata)[2],quad = quad)
      }
      if (numswap > 0) {
        moves <-
          nextswap(
            cur = cur,n = dim(xdata)[2],quad = quad
          )
      }
      if (dim(moves)[2] == 0) {
        moves <- t(t(last))
      }
      if (interactions == T) {
        form <-
          function(j)
            formula(paste0(
              colnames(newdata)[1],"~",paste(fixvar,collapse = "+"),"+",paste(colnames(xdata)[moves[,j]],collapse = "*")
            ),sep = "")
      }
      if (interactions == F) {
        form <-
          function(j)
            formula(paste0(
              colnames(newdata)[1],"~",paste(fixvar,collapse = "+"),"+",paste(colnames(xdata)[moves[,j]],collapse = "+")
            ),sep = "")
      }
      tmp <-
        parallel::mclapply(
          X = 1:dim(moves)[2],FUN = function(k){
            if((sum(complete.cases(cbind(ydata,xdata[,moves[,k]])))/length(ydata))>0.05){
              criterion(glm(
                form(k),data = newdata,family = family,...
              )) 
            } else {NA}
          },mc.cores = cores
        )
      checks <- checks + dim(moves)[2]
      if (minmax == "max") {
        cur <- moves[,which.max.na(unlist(tmp))[1]]
        cur.criterion <- unlist(tmp[which.max.na(unlist(tmp))[1]])
        if (last.criterion > cur.criterion) {
          cur <- last.pos
          cur.criterion <- last.criterion
        }
      }
      if (minmax == "min") {
        cur <- moves[,which.min.na(unlist(tmp))[1]]
        cur.criterion <- unlist(tmp[which.min.na(unlist(tmp))[1]])
        if (last.criterion < cur.criterion) {
          cur <- last.pos
          cur.criterion <- last.criterion
        }
      }
      numswap <- numswap + 1
      last1 <- last
      last.criterion <- cur.criterion
      last.pos <- cur
      memswap <- unique(c(memswap,last1))
    }
    history[i,(1 + m):(2 * m)] <- cur
    history[i,(dim(history)[2] - 2)] <- cur.criterion
    history[i,(dim(history)[2] - 1)] <- numswap - 1
    history[i,(dim(history)[2])] <- checks
    return(history[i,])
  }
  solutions <-
    matrix(unlist(lapply(
      1:numrs,FUN = function(i)
        fsa(i,history)
    )),ncol = dim(history)[2],byrow = T)
  solutions[,1:(2 * m)] <-
    matrix(colnames(newdata)[c(solutions[,1:(2 * m)] + 1)],ncol = (2 * m))
  solutions <- data.frame(solutions)
  colnames(solutions) <-
    c(
      paste("start",1:m,sep = "."),paste("best",1:m,sep = "."),"criterion","swaps","checks"
    )
  solutions$criterion <-
    as.numeric(levels(solutions$criterion))[solutions$criterion]
  solutions$swaps <-
    as.numeric(levels(solutions$swaps))[solutions$swaps]
  solutions$checks <-
    as.numeric(levels(solutions$checks))[solutions$checks]
  if (length(fixvar) != 0) {
    solutions <-
      data.frame(fixvar = matrix(
        rep(x = fixvar,dim(solutions)[1]),nrow = dim(solutions)[1],byrow = T
      ),solutions)
  }
  solutions <- solutions
  a <- solutions[,(length(fixvar) + m + 1):(length(fixvar) + m + 1 + m)]
  b <- unique(t(apply(a,sort,MARGIN = 1)),MARGIN = 1)
  a <- t(apply(a,sort,MARGIN = 1))
  c <- cbind(b,0)
  for (i in 1:dim(b)[1]) {
    for (j in 1:dim(a)[1]) {
      c[i,(m + 2)] <-
        sum(as.numeric(c[i,(m + 2)]) + as.numeric(identical(a[j,],b[i,])))
    }
  }
  tableres <- data.frame(cbind(c),stringsAsFactors = F)
  colnames(tableres)[(dim(tableres)[2])] <- "times"
  colnames(tableres)[2:(dim(tableres)[2] - 1)] <-
    paste("Var",1:m,sep = "")
  colnames(tableres)[1] <- "criterion"

  call <- mget(names(formals()),sys.frame(sys.nframe()))
  ls <-
    list(
      originalfit = fit,call = call,solutions = solutions,table = tableres,efficiency =
        paste(
          "You did:",sum(solutions$checks)," model checks compared to ",choose(n = dim(xdata)[2],k = m)," checks you would have done with exahstive search."
        )
    )
  class(ls) <- "FSA"
  invisible(print(ls))
  return(ls)
  
}
