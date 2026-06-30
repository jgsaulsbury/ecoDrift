#' Probability of a single transition under neutral theory without migration
#'
#' @description
#' Gives the log-likelihood of a single observed transition in community composition
#' under Hubbell's neutral theory, for a given ratio of J (community size) to
#' t (# generations).
#'
#' @details
#' Used by xxprob. Sums of n1 and n2 must each be in the interval (0,1]. Allows
#' for incomplete sampling; if nothing is given for ss, assumes n1 and n2 represent
#' true relative abundances. By default ignores transitions that end in 0.
#'
#' @param n1 vector of relative abundances at start of transition.
#' @param n2 vector of relative abundances at end of transition.
#' @param Jt ratio of J (community size) to t (# generations).
#' @param ss single value or vector of length 2 giving number of samples at
#' times for n1 and n2. If a single value given, assumes it applies to both n1 and n2.
#' @param handle.ext string specifying how extinction (and monodominance) should be
#' handled. Three options: "ignore" simply passes over these transitions.
#' "condition" returns logliks conditional on non-extinction and non-monodominance.
#' "calculate" returns logliks that include species that hit abundance 0 or 1 in
#' the likelihood calculation. "condition" is the default.
#'
#' @importFrom cbinom dcbinom
#'
#' @returns Returns loglik value.
#' @export
#'
#' @examples
#' xprob(n1=c(0.1,0.1,0.8),n2=c(0.2,0.3,0.5),Jt=10,ss=c(1000,1000)) #-0.20906
#'
xprob <- function(n1,n2,Jt,ss=NA,handle.ext="condition"){
  #error handling
  if (is.null(handle.ext) || length(handle.ext) != 1L || !handle.ext %in% c("ignore","condition", "calculate")) {
    stop("'handle.ext' argument must be one of: 'ignore','condition','calculate'")}
  tol <- 1E-7
  if(sum(n1)<=0|sum(n2)<=0|sum(n1)>1+tol|sum(n2)>1+tol){
    stop("Sums of n1 and n2 must each be greater than 0 and no greater than 1")}
  if(length(n1)!=length(n2)){stop("n1 and n2 must have the same length")}
  if(!Jt>0){stop("Jt must be positive")}
  #body
  if(!any(is.na(ss))){ #if ss given...
    if(length(ss)==1){ss <- rep(ss,2)} #duplicate if a single value given
    sizes <- c(ss,Jt)
    size <- 1/sum(1/sizes)
  } else {
    size <- Jt}
  order <- rev(order(n1)) #sort by abundance of n1
  n1 <- n1[order]
  n2 <- n2[order]
  while(sum(n1) > 1-tol){ #remove last species if sum of n1 is 1, not needed for calculation
    n1 <- n1[-length(n1)]
    n2 <- n2[-length(n2)]}
  if(handle.ext != "calculate"){
    n1 <- n1[!(n2==0|n2==1)];n2 <- n2[!(n2==0|n2==1)]} #remove indices where n2 is zero or one (zeros in n1 removed in previous step)
  n2.adj <- n2*size + 0.5 #moving onto scale of cbinom (n2.adj)
  #for taxon 1:
  if(length(n1)==0){ #check to make sure there is at least one transition left
    out <- 0
  } else if(n2[1]==0){ #otherwise if species goes locally extinct...
    out <- cbinom::pcbinom(q=0.5,size=size,prob=n1[1],log=T)
  } else if(n2[1]==1){ #otherwise if species achieves monodominance...
    out <- log(1-cbinom::pcbinom(q=size+0.5,size=size,prob=n1[1]))
  } else { #otherwise...
    out <- cbinom::dcbinom(x=n2.adj[1],size=size,prob=n1[1],log=T)+log(size) #multiplied to make sense for prob densities on [0,1] interval
    #if conditioning on n2[1] being neither 0 nor 1...
    if(handle.ext=="condition"){
      denom <- cbinom::pcbinom(q=size+0.5,size=size,prob=n1[1])-
        cbinom::pcbinom(q=0.5,size=size,prob=n1[1])
      out <- out - log(denom)}
  }
  for(i in seq_len(length(n1))[-1]){ #for every other taxon i
    if(n2[i]==0){ #if species i goes locally extinct...
      out <- out + cbinom::pcbinom(q=0.5,size=size-sum(n2.adj[1:(i-1)]-0.5),
                          prob=n1[i]/(1-sum(n1[1:(i-1)])),log=T)
    } else if(n2[i]==1){ #else if species i achieves monodominance...
      out <- out + log(1-cbinom::pcbinom(q=size-sum(n2.adj[1:(i-1)]-0.5)+0.5,size=size-sum(n2.adj[1:(i-1)]-0.5),
                          prob=n1[i]/(1-sum(n1[1:(i-1)]))))
    } else{ #otherwise...
      out <- out + cbinom::dcbinom(x=n2.adj[i],size=size-sum(n2.adj[1:(i-1)]-0.5),
                          prob=n1[i]/(1-sum(n1[1:(i-1)])),log=T)+log(size)
      #if conditioning on n2[i] being neither 0 nor 1...
      if(handle.ext=="condition"){
        #get value of the CDF 0.5 units from the left and right boundaries
        cdf <- cbinom::pcbinom(q=c(0.5,size-sum(n2.adj[1:(i-1)]-0.5)+0.5),
                               size=size-sum(n2.adj[1:(i-1)]-0.5),
                               prob=n1[i]/(1-sum(n1[1:(i-1)])))
        denom <- diff(cdf)
        out <- out - log(denom)}
    }
  }
  return(out)}
