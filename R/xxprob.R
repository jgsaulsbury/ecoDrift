#' Probability of a set of transitions under neutral theory
#'
#' @description
#' Gives the likelihood of a community composition timeseries under Hubbell's
#' neutral theory for a given value of J. Unlike its helper function xprob,
#' xxprob takes data in the form of the number of individuals observed in each
#' species, rather than relative abundances. Can take either a single dataset or
#' a list of them.
#'
#' @details
#' Assumes non-overlapping generations. Because the user supplies the ages of each
#' time slice, the function needs to be given a generation time to know how many
#' generations separate each adjacent pair of time slices. Default is 1 year.
#'
#' @param log10J log J (community size). Easier for fitJ() to optimize in log space.
#' @param occs matrix of the number of observations in each species at each time.
#' One column for each species, one row for each time slice. Time goes from oldest
#' at the bottom to youngest at the top.
#' @param ages vector containing the ages of each time slice, in years, from
#' oldest to youngest.
#' @param sampled boolean indicating whether occs represents a sampled
#' community (TRUE) or instead represents true species abundances (FALSE).
#' @param generationtime time between generations, in years.
#'
#' @returns loglik value.
#' @export
#'
#' @examples
#' mat <- matrix(data=c(52,12,160,109,30,401,93,31,355),nrow=3)
#' xxprob(log10J=5,occs=mat,ages=c(200,100,0))
xxprob <- function(log10J,occs,ages,sampled=TRUE,generationtime=1){
  if(!is.list(ages)){ #if there's just one timeseries
    occs <- list(occs) #make it the only member of a list
    ages <- list(ages)} #and do the same to ages
  if(length(occs)!=length(ages)){stop("number of occurrence matrices and ages vectors do not match")}
  loglik <- 0
  for(i in length(occs)){ #for each member of the list of occurrence tables
    occ <- occs[[i]]
    age <- ages[[i]]
    if(dim(occ)[1] != length(age)){stop(paste("age vector",i,"does not match number of rows in occurrence matrix",i))}
    ss <- rowSums(occ)
    occs.prop <-  occ/ss #from occurrences to proportional abundance
    for(i in rev(seq(dim(occ)[1]-1))){ #for every transition (from oldest to youngest)
      t = abs(age[i+1]-age[i])/generationtime
      loglik <- loglik + ifelse(sampled,
                                xprob(n1=as.numeric(occs.prop[i+1,]),n2=as.numeric(occs.prop[i,]),Jt=(10^log10J)/t,ss=c(ss[i+1],ss[i])),
                                xprob(n1=as.numeric(occs.prop[i+1,]),n2=as.numeric(occs.prop[i,]),Jt=(10^log10J)/t,ss=NA))}}
  return(loglik)}
