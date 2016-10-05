#' Variables to include in first steip of an mth order interaction model determined from the Feasible Soution Alorithm.
#'
#' Finds the swaps available given a current position.
#' @param cur A vector of length greater than 2 of what current explantory varialbes are being used in the model.
#' @param n The number of explanatory variables in available to swap.
#' @param quad Whether to include quadratic terms. ie (x1*x1) as potential swaps.
#' @return a matrix with the possible forms by column.
#' @export
swaps<- function(cur,n,quad=FALSE){
  m <- length(cur)
  l <- (n-m)*m
   possible<- do.call(rbind,parallel::mclapply(m:1,FUN=function(j){
      possible <- matrix(rep(cur,l/m),nrow = l/m,byrow = T)
      possible[,j]<-(1:n)[-which((1:n) %in% cur )]
      possible
    },mc.cores=m))
  if(quad & m==2){
   possible<-rbind(c(cur[1],cur[1]),possible,c(cur[2],cur[2]))
  }
  vec<-cbind(cur,t(possible))
  return(unique(vec,MARGIN = 2))
}

#' Variables to include in the >1st step of an mth order interaction model determined from the Feasible Soution Alorithm.
#'
#' Finds the swaps available given a current position given previous picks.
#' @param cur A vector of length greater than 2 of what current explantory varialbes are being used in the model.
#' @param n The number of explanatory variables in available to swap.
#' @param quad Whether to include quadratic terms. ie (x1*x1) as potential swaps.
#' @param prevpos A vector of previous best spots
#' @return a matrix with the possible forms by column.
#' @export
nextswap<-function(curpos,n,prevpos,quad=FALSE){
  swps<-swaps(curpos,n,quad)
  nextpos<-colSums(do.call(rbind,parallel::mclapply(1:dim(swps)[1],FUN = function(qq){
    (swps[qq,] %in% prevpos)
  },mc.cores=dim(swps)[1])))
  
  retSwps<-swps[,nextpos==(length(curpos)-2)]
  if(is.null(dim(retSwps))){retSwps<-t(t(retSwps))}
  return(list(nswaps=retSwps,prevpos=prevpos))
}
