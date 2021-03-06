% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/glmFSA.R
\name{glmFSA}
\alias{glmFSA}
\title{rFSA: Feasible Solution Algorithm (FSA) for Generalized Linear Models}
\usage{
glmFSA(formula, data, fixvar = NULL, quad = FALSE, m = 2, numrs = 1,
  cores = 1, interactions = TRUE, criterion = AIC, minmax = "min",
  family = "binomial", checkfeas = NULL, ...)
}
\arguments{
\item{formula}{an object of class "formula" (or one that can be coerced to that class): a symbolic description of the model to be fitted. See help(glm) for details.}

\item{data}{a data frame, list or environment (or object coercible by as.data.frame to a data frame) containing the variables in the model.}

\item{fixvar}{a variable to fix in the model. Usually a covariate that should always be included (Example: Age, Sex). Will still consider it with interactions. Default is NULL.}

\item{quad}{to include quadratic terms or not.}

\item{m}{order of terms to potentially include. If interactions is set to TRUE then m is the order of interactions to be considered. Defaults to 2. For Subset selection (interaction=F), m is the size of the subset to examine. Default is 2.}

\item{numrs}{number of random starts to perform.}

\item{cores}{number of cores to use while running. Note: Windows can only use 1 core. See mclapply for details.  If function detects a Windows user it will automatically set cores=1.}

\item{interactions}{T or F for whether to include interactions in model. Defaults to FALSE.}

\item{criterion}{which criterion function to either maximize or minimize. For linear models one can use: apress (Allens Press Statistic), int.p.val (Interaction p-value), AIC, BIC.}

\item{minmax}{whether to minimize or maximize the criterion function}

\item{family}{family argument passed to glm. A description of the error distribution and link function to be used in the model. This can be a character string naming a family function, a family function or the result of a call to a family function.}

\item{checkfeas}{A vector of a potential feasible solution}

\item{...}{arguments to be passed to the glm function}
}
\value{
returns a list of solutions and table of unique solutions.
$solutions is a matrix of fixed terms, start position, feasible solution, criterion function value (p-value of interaction), and number of swaps to solution.
$table is a matrix of the unique feasible solutions and how many times they occurred out of the number of random starts chosen. It also returns any warning messages with these solutions in the last column.
$efficiency is text comparing how many models you ran during your FSA search compared to how many you would have done with exhaustive search. Note: The FSA algorithm takes additional time to run on top of the model checks that were done during the algorithm. This additional time is approximately 15% more time than if you had just ran the model checks.
}
\description{
A function using a Feasible Solution Algorithm to find a set of feasible solutions for a generalized linear model of a specific form that could include mth-order interactions (Note that these solutions are optimal in the sense that no one swap to any of the variables will increase the criterion function.)
}
\details{
PLEASE NOTE: make sure categorical variables are factors or characters otherwise answers will not reflect the variable being treated as a continuous variable.
}
\examples{
dat<-read.csv("http://tinyurl.com/zq7l775",header = FALSE)
colnames(dat)<-c("Class","Age","Sex","Sterioid","Antivirals",
                 "Fatigue","Malaise","Anorexia","Liver Big",
                 "Liver Firm","Spleen Palpable","Spiders",
                 "Ascites","Varices","Bilirubin","Alk Phosphate",
                 "Sgot","Albumin","Protime","Histology")
dat<-as.matrix(dat)
dat[which(dat=="?")]=NA
dat<-data.frame(dat)
dat[,c(2,15,16,17,18,19)]<-lapply(X = dat[,c(2,15,16,17,18,19)],as.numeric)
colnames(dat)
fit<-glmFSA(formula="Class~Age+Sgot*Albumin",data=dat,fixvar="Age",quad=FALSE,m=2,
       numrs=10,family="binomial",cores=1)
}

