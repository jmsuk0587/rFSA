#' Random starts with m variables
#'
#' @param m order of interaction to include.
#' @param nvars number of varialbes that could be included in the interaction.
#' @param quad whether to include quadratic terms or not.
#' @param numrs number of random starts to compute.
#'
#' @return A matrix of random starts by row
#' @export
#'
#' @examples
#' rstart(2,5,quad=F,numrs=3) #doesn't include quadratic terms
#' rstart(2,5,quad=T,numrs=3) #does include quadratic terms
rstart=function(m,nvars,quad=F,numrs=1){# m is number of variables to include, nvars is a matrix of variables
  t(replicate(numrs,sample(1:nvars,size=m,replace=quad)))
}
