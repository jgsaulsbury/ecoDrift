#' Find the maximum-likelihood value of J for a community change dataset
#'
#' @description
#' Uses optimize() to fit the J parameter in Hubbell's neutral theory to a dataset.
#' Optionally, outputs a confidence interval.
#'
#' @details
#' optimize() sometimes tests J values which are so unlikely as to return -Inf.
#' This is fine, and suppressWarnings() just keeps this quiet.
#'
#'
#' @param occs matrix of the number of observations in each species at each time.
#' One column for each species, one row for each time slice. Time goes from oldest
#' at the bottom to youngest at the top.
#' @param ages vector containing the ages of each time slice, in years, from
#' oldest to youngest. Numbers can run from high to low (e.g., millions of years
#' ago) or from low to high (e.g., years from the start of a simulation).
#' @param sampled boolean indicating whether occs represents a sampled
#' community (TRUE) or instead represents true species abundances (FALSE).
#' @param generationtime time between generations, in years. Default is 1.
#' @param CI boolean indicating whether or not to calculate a confidence interval
#' on the maximum-likelihood estimate of J. Uses function CIfunc_J.
#' @param searchinterval consider values of log10J within this interval.
#' @param CI.df integer giving the number of degrees of freedom to use when
#' using the chi-square lookup table to calculate confidence intervals. Should be
#' 1 unless calculating a simultaneous confidence interval.
#'
#' @returns a list containing "loglik", "J", and (optionally) "CI"
#' @export
#'
#' @examples
#' mat <- matrix(data=c(520,1200,1600,1090,900,401,930,610,355),nrow=3)
#' fitJ(occs=mat,ages=c(200,100,0),CI=TRUE)
fitJ <- function(occs,ages,sampled=TRUE,generationtime=1,CI=FALSE,searchinterval=c(1,9),CI.df=1){
  if(dim(occs)[1] != length(ages)){stop("'ages' must have length equal to the number of rows of 'occs'")}
  op <- suppressWarnings(stats::optimize(xxprob,interval=c(searchinterval[1],searchinterval[2]),occs=occs,ages=ages,sampled=sampled,generationtime=generationtime,maximum=TRUE))
  out <- list("loglik"=op$objective,"J"=10^op$maximum)
  if(CI){
    left <- suppressWarnings(stats::optimize(CIfunc_J,interval=c(searchinterval[1],log10(out$J)),occs=occs,ages=ages,ML=out$loglik,
                            sampled=sampled,generationtime=generationtime,CI.df=CI.df,maximum=FALSE))$minimum
    right <- suppressWarnings(stats::optimize(CIfunc_J,interval=c(log10(out$J),searchinterval[2]),occs=occs,ages=ages,ML=out$loglik,
                             sampled=sampled,generationtime=generationtime,CI.df=CI.df,maximum=FALSE))$minimum
    out$CI <- 10^c(left,right)}
  return(out)}
