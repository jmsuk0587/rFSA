#' @export
r.squared<-function(model){
  summary(model)$r.squared
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
  tail(fit$coefficients,n = 1)
}