#' Find the maximum-likelihood values of J and m for a community change dataset
#'
#' @description
#' Uses optim() to fit the J and m parameters in Hubbell's neutral theory to
#' a dataset, where J is community size and m is the rate at which vacancies in
#' the community are filled by migrants from a static metacommunity.
#'
#' @details
#' optim() takes the midpoint of the J and m search bounds as the starting value.
#'
#' @param occs matrix of the number of observations in each species at each time.
#' One column for each species, one row for each time slice. Time goes from oldest
#' at the bottom to youngest at the top.
#' @param ages vector containing the ages of each time slice, in years, from
#' oldest to youngest.
#' @param sampled boolean indicating whether occs represents a sampled
#' community (TRUE) or instead represents true species abundances (FALSE).
#' @param metacommunity vector of relative abundances of species in the metacommunity.
#' If this doesn't sum to 1, will be normalized to sum to 1.
#' @param generationtime time between generations, in years.
#'
#' @returns a list containing "loglik", "J", and "m"
#'
#' @examples
#' #simulate under neutral theory with migration
#' set.seed(10)
#' sim <- simDrift(c(1000,1000,1000,1000),ts=seq(0,2000,50),m=0.001,ss=1000)
#' ecoDrift::fitJm(occs=sim$simulation,ages=sim$times,metacommunity=rep(0.25,4))
fitJm <- function(occs,ages,metacommunity,sampled=TRUE,generationtime=1){
  metacommunity <- metacommunity/sum(metacommunity) #make it sum to 1
  if(dim(occs)[1] != length(ages)){stop("'ages' must have length equal to the number of rows of 'occs'")}
  if(dim(occs)[2] != length(metacommunity)){stop("length of 'metacommunity' must be equal to number of columns in 'occs'")}
  op <- stats::optim(par=c(5,-3),fn=xxprobm,method="Nelder-Mead",control=list(fnscale=-1),occs=occs,ages=ages,
                     sampled=sampled,metacommunity=metacommunity,generationtime=generationtime)
  out <- list("loglik"=op$value,"J"=10^op$par[1],"m"=10^op$par[2])
  return(out)}
