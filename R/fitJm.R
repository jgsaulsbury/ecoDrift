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
#' @param metacommunity vector of relative abundances of species in the metacommunity,
#' or list of several such vectors. If this doesn't sum to 1, will be normalized
#' to sum to 1. If metacommunity is left NA, xxprobm will use abundance in first timestep
#' as a guess.
#' @param m float (optional) to constrain m to
#' @param generationtime time between generations, in years.
#' @param condition.nonext boolean indicating whether likelihoods should be conditioned
#' on n2 not including any 0s or 1 (local extinction or monodominance). Passed to xprobm.
#' TRUE is default.
#'
#' @returns a list containing "loglik", "J", and "m"
#'
#' @examples
#' #simulate under neutral theory with migration
#' set.seed(10)
#' sim <- simDrift(c(1000,1000,1000,1000),ts=seq(0,2000,50),m=0.001,ss=1000)
#' ecoDrift:::fitJm(occs=sim$simulation,ages=sim$times,metacommunity=rep(0.25,4))
fitJm <- function(occs,ages,metacommunity=NA,m=NA,sampled=TRUE,generationtime=1,condition.nonext=TRUE){
  if(is.na(m)){
    op <- stats::optim(par=c(5,-7),fn=xxprobm,method="Nelder-Mead",control=list(fnscale=-1),occs=occs,ages=ages,
                       sampled=sampled,metacommunity=metacommunity,generationtime=generationtime)
    out <- list("loglik"=op$value,"J"=10^op$par[1],"m"=10^op$par[2])
  } else {
    op <- suppressWarnings(stats::optimize(function(x)
      xxprobm(log10Jm=c(x,log10(m)),occs=occs,ages=ages,sampled=sampled,generationtime=generationtime,condition.nonext=condition.nonext),
                interval=c(1,9),maximum=TRUE))
    out <- list("loglik"=op$objective,"J"=10^op$maximum,"m"=m)
  }

  return(out)}
