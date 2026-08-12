#' Probability of a set of transitions under neutral theory with migration
#'
#' @description
#' Gives the likelihood of a community composition timeseries under Hubbell's
#' neutral theory for given values of J and m. Unlike its helper function xprobm,
#' xxprobm takes data in the form of the number of individuals observed in each
#' species, rather than relative abundances. Can take either a single dataset or
#' a list of them.
#'
#' @details
#' Assumes non-overlapping generations. Because the user supplies the ages of each
#' time slice, the function needs to be given a generation time to know how many
#' generations separate each adjacent pair of time slices. Default is 1 year.
#'
#' @param log10Jm a 2-length vector containing log10 J (community size) and log10 m
#' (migration rate).
#' @param occs matrix of the number of observations in each species at each time,
#' or else a list containing several such matrices.
#' In each matrix there is one column per species and one row per time slice.
#' Time goes from oldest at the bottom to youngest at the top.
#' @param ages vector containing the ages of each time slice, in years, from
#' oldest to youngest, or else a list of such vectors if there are multiple
#' occurrence matrices.
#' @param sampled boolean indicating whether occs represents a sampled
#' community (TRUE) or instead represents true species abundances (FALSE).
#' @param metacommunity relative abundances of each species in the metacommunity,
#' given in the same order as in occs, or else a list of such vectors. If metacommunity
#'  is left NA, xxprobm will use abundance in first timestep as a guess.
#' @param generationtime time between generations, in years.
#' @param condition.nonext boolean indicating whether likelihoods should be conditioned
#' on n2 not including any 0s or 1 (local extinction or monodominance). Passed to xprobm.
#' TRUE is default.
#' @returns loglik value.
#'
#' @examples
#' set.seed(10)
#' #simulate under neutral theory with migration
#' sim <- simDrift(c(1000,1000,1000,1000),ts=seq(0,2000,50),m=0.005,ss=10000)
#' par(mfrow=c(1,2))
#' plot_spindles(occs=sim$simulation,ages=sim$times,plot.ss=FALSE,revertsettings=FALSE)
#' #plot likelihood surface for migration rate
#' ms <- seq(-7,-1,0.1) #get likelihood for these values of m
#' liks <- c()
#' for(i in ms){
#'   liks <- c(liks,
#'     ecoDrift:::xxprobm(log10Jm = c(log10(4000),i),occs=sim$simulation,ages=sim$times,sampled=TRUE))}
#' plot(10**ms,liks,ylim=c(100,300),xlab="m",ylab="loglik",type='l',log='x')
#' lines(c(0.005,0.005),c(0,400),lty='dashed') #true m
xxprobm <- function(log10Jm,occs,ages,metacommunity,sampled=TRUE,generationtime=1,condition.nonext=TRUE){
  #catching errors
  if(!is.list(ages)){ #if there's just one timeseries
    occs <- list(occs) #make it the only member of a list
    ages <- list(ages) #and do the same to ages
    metacommunity <- list(metacommunity)} #and to metacommunity
  if(length(occs)!=length(ages)){stop("number of occurrence mats and ages vectors do not match")}
  if(length(metacommunity)==1){ #if only one metacommunity provided, duplicate for all timeseries
    metacommunity <- replicate(length(occs),metacommunity)}
  if(log10Jm[2]>0|10**log10Jm[1]==0){return(-Inf)} #to use optim() without constraining param vals
  #body
  loglik <- 0
  for(i in length(occs)){ #for each member of the list of occurrence tables
    occ <- occs[[i]]
    age <- ages[[i]]
    meta <- metacommunity[[i]]
    if(dim(occ)[1] != length(age)){stop(paste("age vector",i,
                                          "does not match number of rows in occurrence matrix",i))}
    ss <- rowSums(occ)
    occs.prop <-  occ/ss #from occurrences to proportional abundance
    if(any(is.na(meta))){ #if metacommunity composition not provided...
      meta <- occs.prop[length(ss),]} #use local abundance in the first timestep as a guess
    if(dim(occ)[2] != length(meta)){stop(paste("length of metacommunity",i,
                                               "does not match number of columns in occurrence timeseries"))}
    for(i in rev(seq(dim(occ)[1]-1))){ #for every transition (from oldest to youngest)
      t = abs(age[i+1]-age[i])/generationtime
      loglik <- loglik + ifelse(sampled,
                          xprobm(n1=as.numeric(occs.prop[i+1,]),n2=as.numeric(occs.prop[i,]),
                            nmeta=meta,J=10**log10Jm[1],m=10**log10Jm[2],t=t,ss=c(ss[i+1],ss[i]),condition.nonext=condition.nonext),
                          xprobm(n1=as.numeric(occs.prop[i+1,]),n2=as.numeric(occs.prop[i,]),
                            nmeta=meta,J=10**log10Jm[1],m=10**log10Jm[2],t=t,ss=NA),condition.nonext=condition.nonext)}}
  return(loglik)}
