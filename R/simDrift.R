#' Simulate neutral community evolution
#'
#' @description
#' Simulate drift in a neutral community with non-overlapping generations. Includes
#' the option to stabilize the community via immigration from a static metacommunity.
#'
#' @param startingabs vector containing the starting abundance for each species.
#' @param ts vector of times at which to return the state of the community.
#' @param ss sample size for each time slice. Sampling is with replacement.
#' @param m rate of migration into the community from the metacommunity.
#' @param metacommunity relative abundances for species in the static metacommunity.
#' If none provided, uses the starting abundances for the local community.
#'
#' @returns a list with entries "simulation", a matrix containing the output of
#' the simulation with a row for each time slice (time goes from bottom to top)
#' and each column is a species; and "times", a vector that passes on the value of
#' ts.
#'
#' @export
#'
#' @examples
#' #First try with no metacommunity
#' J <- 1000
#' tslength <- 500
#' every <- 10
#' nsp <- 10
#' ages <- seq(0,tslength,every)
#' timeseries <- simDrift(startingabs=rep(J/nsp,nsp),ts=ages,ss=1000)
#' plot_spindles(occs=timeseries$simulation,ages=ages,linesevery=100)
#' #then try with migration
#' timeseries <- simDrift(startingabs=rep(J/nsp,nsp),ts=ages,ss=1000,m=0.01)
#' plot_spindles(occs=timeseries$simulation,ages=ages,linesevery=100)
simDrift <- function(startingabs,ts,ss=NA,m=0,metacommunity=NA){
  J <- sum(startingabs)
  ts <- ts-min(ts) #in case ts doesn't start at 0
  if(m>0 & any(is.na(metacommunity))){ #if no metacommunity provided...
    metacommunity <- startingabs/J} #use the starting abundances
  out <- list("simulation"=matrix(NA,ncol=length(startingabs),nrow=length(ts)),"times"=ts)
  colnames(out$simulation) <- paste("Species",seq(length(startingabs))) #give species names
  if(!is.na(ss)){
    out$simulation[length(ts),] <- stats::rmultinom(n=1,size=ss,prob=c(startingabs)/J)
  } else {
    out$simulation[length(ts),] <- startingabs}
  for(time in seq(max(ts))){
    migrants <- ifelse(m>0,stats::rbinom(1,J,m),0) #draw this many individuals from the metacommunity
    startingabs <- stats::rmultinom(n=1,size=J-migrants,prob=c(startingabs)/J)
    if(migrants>0){ #if there are migrants
      startingabs <- startingabs + stats::rmultinom(n=1,size=migrants,prob=metacommunity)}
    if(time %in% ts){ #add to the output if it's time
      if(!is.na(ss)){
        samp <- stats::rmultinom(n=1,size=ss,prob=c(startingabs)/J) #downsample
        out$simulation[length(ts)-which(ts==time)+1,] <- samp
      } else {
        out$simulation[length(ts)-which(ts==time)+1,] <- startingabs
      }}}
  return(out)}
