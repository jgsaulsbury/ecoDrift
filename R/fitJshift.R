#' Fit a drift model with a shift in J
#'
#' @description
#' Takes a time-series and fits a 3-parameter model where J shifts from one value
#' to another at some point in the time-series. Also can output confidence intervals
#' on all three parameters.
#'
#' @details
#' optimize() sometimes tests J values which are so unlikely as to return -Inf.
#' This is fine, and suppressWarnings() just keeps this quiet. The confidence
#' interval on the shift age returns the oldest and youngest ages that yield likelihoods
#' within a specified distance of the maximum likelihood. Thus, there may be
#' shift ages within the confidence interval that are outside the 95% confidence
#' interval if the likelihood surface is not convex.
#'
#' @param occs matrix of the number of observations in each species at each time.
#' One column for each species, one row for each time slice. Time goes from oldest
#' at the bottom to youngest at the top.
#' @param ages vector containing the ages of each time slice, in years, from
#' oldest to youngest. Numbers can run from high to low (e.g., millions of years
#' ago) or from low to high (e.g., years from the start of a simulation)
#' @param liksurf boolean indicating whether to include the likelihood surface
#' for the shift age at each point along ages.
#' @param sampled boolean indicating whether occs represents a sampled
#' community (TRUE) or instead represents true species abundances (FALSE).
#' @param generationtime time between generations, in years. Default is 1.
#' @param CI boolean indicating whether or not to calculate confidence intervals
#' on the maximum-likelihood estimates of J1, J2, and the shift age. Uses function
#' CIfunc_J to find 95% confidence intervals for J. Provides simultaneous confidence
#' intervals: in theory, for a 95% confidence interval there is a 95% chance that
#' all parameter estimates are correct.
#' @param searchinterval_J consider values of log10J within this interval.
#'
#' @returns a list containing "loglik", "J1", "J2", "shiftage" and (optionally) three confidence
#' intervals: "CI_J1", "CI_J2", and "CI_shiftage"
#'
#' @export
#'
#' @examples
#' set.seed(1)
#' #first simulate a timeseries with an increase in J halfway through
#' J1 <- 1000 #first half of TS
#' J2 <- 10000 #second half of TS
#' tslength <- 500
#' every <- 10
#' nsp <- 8
#' ss <- 1000000
#' ages <- seq(0,tslength,every)
#' timeseries1 <- simDrift(startingabs=rep(J1/nsp,nsp),ts=ages[1:(length(ages)/2 + 1)],ss=ss)$simulation #first half
#' timeseries2 <- simDrift(startingabs=timeseries1[1,]*J2/ss,ts=ages[1:(length(ages)/2 + 1)],ss=ss)$simulation #second half
#' timeseries <- rbind(timeseries2,timeseries1[-1,])
#' par(mfrow=c(1,2))
#' plot_spindles(timeseries,ages,plot.ss=FALSE)
#' lines(c(-1,10),c(tslength/2,tslength/2),lty='dashed',lwd=2) #plot point where J switches
#' plot_Js(timeseries,ages) #fit rate through time
#' lines(c(1E-10,1E-1),c(tslength/2,tslength/2),lty='dashed',lwd=2)
#' #then use fitJshift
#' fit <- fitJshift(timeseries,ages,CI=TRUE)
#' #plot uncertainty envelope then ML estimate
#' #points go clockwise from bottom left
#' polygon(x=1/c(fit$CI_J1[2],fit$CI_J1[2],fit$CI_J2[2],fit$CI_J2[2],
#'   fit$CI_J2[1],fit$CI_J2[1],fit$CI_J1[1],fit$CI_J1[1]),
#'   y=c(ages[1],fit$CI_shiftage[1],fit$CI_shiftage[1],ages[length(ages)],
#'   ages[length(ages)],fit$CI_shiftage[2],fit$CI_shiftage[2],ages[1]),border=NA,col=rgb(0.7, 0.7, 0.7, 0.5))
#' lines(1/c(fit$CI_J1[2],fit$CI_J1[2],fit$CI_J2[2],fit$CI_J2[2]),
#'   c(ages[1],fit$CI_shiftage[1],fit$CI_shiftage[1],ages[length(ages)]),lwd=1,col="grey50")
#' lines(1/c(fit$CI_J2[1],fit$CI_J2[1],fit$CI_J1[1],fit$CI_J1[1]),
#'   c(ages[length(ages)],fit$CI_shiftage[2],fit$CI_shiftage[2],ages[1]),lwd=1,col="grey50")
#' lines(1/c(fit$J1,fit$J1,fit$J2,fit$J2),c(min(ages),fit$shiftage,fit$shiftage,max(ages)),lwd=2)
fitJshift <- function(occs,ages,liksurf=FALSE,sampled=TRUE,generationtime=1,CI=FALSE,searchinterval=c(1,9)){
  if(dim(occs)[1] != length(ages)){stop("'ages' must have length equal to the number of rows of 'occs'")}
  bestfits <- list() #list of best-fit models for every shift age (no CIs yet)
  #first age
  noshift <- fitJ(occs,ages,sampled=sampled,generationtime=generationtime,CI=FALSE,searchinterval=searchinterval)
  bestfits <- append(bestfits,list(list("loglik"=noshift$loglik,"J1"=NA,"J2"=noshift$J,"shiftage"=ages[1])))
  for(i in 2:(length(ages)-1)){ #consider all shift ages except first and last ages
    #lower half of occs
    occs1 <- occs[(length(ages)+1-i):length(ages),]
    ages1 <- ages[1:i]
    fit1 <- fitJ(occs1,ages1,sampled=sampled,generationtime=generationtime,CI=FALSE,searchinterval=searchinterval)
    #upper half of occs
    occs2 <- occs[1:(length(ages)+1-i),]
    ages2 <- ages[i:length(ages)]
    fit2 <- fitJ(occs2,ages2,sampled=sampled,generationtime=generationtime,CI=FALSE,searchinterval=searchinterval)
    #append to bestfits
    bestfits <- append(bestfits,list(list("loglik"=fit1$loglik+fit2$loglik,"J1"=fit1$J,"J2"=fit2$J,"shiftage"=ages[i])))}
  #last age
  bestfits <- append(bestfits,list(list("loglik"=noshift$loglik,"J1"=noshift$J,"J2"=NA,"shiftage"=ages[length(ages)])))
  #ML estimate
  liksurf <- unlist(lapply(bestfits, `[[`, 1))
  out <- bestfits[[which.max(liksurf)]]
  out$liksurf <- liksurf
  if(CI){
    #cut occs in half again
    occs1 <- occs[(length(ages)+1-which.max(liksurf)):length(ages),]
    ages1 <- ages[1:which.max(liksurf)]
    occs2 <- occs[1:(length(ages)+1-which.max(liksurf)),]
    ages2 <- ages[which.max(liksurf):length(ages)]
    #this is slightly inefficient
    fit1 <- fitJ(occs1,ages1,sampled=sampled,generationtime=generationtime,CI=TRUE,searchinterval=searchinterval,CI.df=3)
    fit2 <- fitJ(occs2,ages2,sampled=sampled,generationtime=generationtime,CI=TRUE,searchinterval=searchinterval,CI.df=3)
    out$CI_J1 <- fit1$CI
    out$CI_J2 <- fit2$CI
    threshold <- stats::qchisq(p=0.95,d=3)/2 #likelihood threshold for CI
    shiftage_indices <- which(max(liksurf)-liksurf < threshold)
    out$CI_shiftage <- c(ages[head(shiftage_indices,1)],ages[tail(shiftage_indices,1)])}
  return(out)}
