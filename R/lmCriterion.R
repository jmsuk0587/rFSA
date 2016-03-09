r.squared<-function(model){
  summary(model)$r.squared
}

adj.r.squared<-function(model){
  summary(model)$adj.r.squared
}

cv5.lmFSA<-function(model){
  as.numeric(tail(capture.output(cv.lm(data=dat,form.lm = model,m = 5,printit =T,plotit = F)),n=1))
}
cv10.lmFSA<-function(model){
  as.numeric(tail(capture.output(cv.lm(data=dat,form.lm = model,m = 10,printit =T,plotit = F)),n=1))
}

apress<-function(model){
  press(model)
}


which.max.na<-function(vec){
  maxval<-max(vec,na.rm=T)
  which(vec==maxval)
}

which.min.na<-function(vec){
  minval<-min(vec,na.rm=T)
  which(vec==minval)
}