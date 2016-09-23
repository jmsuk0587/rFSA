#' @export
r.squared<-function(model){
  summary(model)$r.squared
}

#' @export
rmse<-function(model){
 sqrt(mean(model$residuals^2))
}

#' @export
adj.r.squared<-function(model){
  summary(model)$adj.r.squared
}

#' @export
cv5.lmFSA<-function(model){
  as.numeric(tail(capture.output(cv.lm(data=dat,form.lm = model,m = 5,printit =T,plotit = F)),n=1))
}

#' @export
cv10.lmFSA<-function(model){
  as.numeric(tail(capture.output(cv.lm(data=dat,form.lm = model,m = 10,printit =T,plotit = F)),n=1))
}

#' @export
apress<-function(model){
  press(model)
}

#' @export
which.max.na<-function(vec){
  maxval<-max(vec,na.rm=T)
  which(vec==maxval)
}

#' @export
which.min.na<-function(vec){
  minval<-min(vec,na.rm=T)
  which(vec==minval)
}

#' @export
int.p.val<-function(model){
show(tail(anova(model,test="LRT")$`Pr(>Chi)`,1))
  tail(anova(model,test="LRT")$`Pr(>Chi)`,1)
}

#' @export
bdist <- function(model){
tmp_dat<-eval(model$call$data) 
y<-tmp_dat[,all.vars(formula(model))[1]]
x1<-tmp_dat[,all.vars(formula(model))[2]]
x2<-tmp_dat[,all.vars(formula(model))[3]]

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

