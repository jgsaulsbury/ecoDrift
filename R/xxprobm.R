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
#' @param occs matrix of the number of observations in each species at each time.
#' One column for each species, one row for each time slice. Time goes from oldest
#' at the bottom to youngest at the top.
#' @param ages vector containing the ages of each time slice, in years, from
#' oldest to youngest.
#' @param sampled boolean indicating whether occs represents a sampled
#' community (TRUE) or instead represents true species abundances (FALSE).
#' @param metacommunity relative abundances of each species in the metacommunity,
#' given in the same order as in occs. If none given, takes the average species
#' abundances across the timeseries as a guess.
#' @param generationtime time between generations, in years.
#'
#' @returns loglik value.
#'
#' @examples
#' set.seed(10)
#' #simulate under neutral theory with migration
#' sim <- simDrift(c(1000,1000,1000,1000),ts=seq(0,2000,50),m=0.005,ss=10000)
#' par(mfrow=c(1,2))
#' plot_spindles(occs=sim$simulation,ages=sim$times,plot.ss=FALSE)
#' #plot likelihood surface for migration rate
#' ms <- seq(-7,-1,0.1) #get likelihood for these values of m
#' liks <- c()
#' for(i in ms){
#'   liks <- c(liks,
#'     ecoDrift::xxprobm(log10Jm = c(log10(4000),i),occs=sim$simulation,ages=sim$times,sampled=TRUE))}
#' plot(10**ms,liks,ylim=c(100,300),xlab="m",ylab="loglik",type='l',log='x')
#' lines(c(0.005,0.005),c(0,400),lty='dashed') #true m
xxprobm <- function(log10Jm,occs,ages,sampled=TRUE,metacommunity=NA,generationtime=1){
  if(!is.list(ages)){ #if there's just one timeseries
    occs <- list(occs) #make it the only member of a list
    ages <- list(ages)} #and do the same to ages
  if(length(occs)!=length(ages)){stop("number of occurrence mats and ages vectors do not match")}
  if(log10Jm[2]>0|10**log10Jm[1]==0){return(-Inf)} #to use optim() without constraining param vals
  loglik <- 0
  for(i in length(occs)){ #for each member of the list of occurrence tables
    occ <- occs[[i]]
    age <- ages[[i]]
    if(dim(occ)[1] != length(age)){stop(paste("age vector",i,
                                          "does not match number of rows in occurrence matrix",i))}
    ss <- rowSums(occ)
    occs.prop <-  occ/ss #from occurrences to proportional abundance
    if(any(!is.na(metacommunity))){ #if metacommunity composition given...
      nmeta <- metacommunity} else{
      nmeta <- occs.prop[length(ss),]} #local abundance in the first timestep as a guess at meta
    for(i in rev(seq(dim(occ)[1]-1))){ #for every transition (from oldest to youngest)
      t = abs(age[i+1]-age[i])/generationtime
      loglik <- loglik + ifelse(sampled,
                          xprobm(n1=as.numeric(occs.prop[i+1,]),n2=as.numeric(occs.prop[i,]),
                            nmeta=nmeta,J=10**log10Jm[1],m=10**log10Jm[2],t=t,ss=c(ss[i+1],ss[i])),
                          xprobm(n1=as.numeric(occs.prop[i+1,]),n2=as.numeric(occs.prop[i,]),
                            nmeta=nmeta,J=10**log10Jm[1],m=10**log10Jm[2],t=t,ss=NA))}}
  return(loglik)}
