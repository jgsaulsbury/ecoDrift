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
#' @param ignore.ext boolean indicating whether transitions that end in 0 should be excluded
#' from likelihood calculation. TRUE by default.
#'
#' @importFrom cbinom dcbinom
#'
#' @returns Returns loglik value.
#' @export
#'
#' @examples
#' xprob(n1=c(0.1,0.1,0.8),n2=c(0.2,0.3,0.5),Jt=10,ss=c(1000,1000)) #-0.20906
#'
xprob <- function(n1,n2,Jt,ss=NA,ignore.ext=TRUE){
  #error handling
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
  if(ignore.ext){
    n1 <- n1[!n2==0];n2 <- n2[!n2==0]} #remove indices where n2 is zero (zeros in n1 removed in previous step)
  n2 <- n2*size + 0.5 #moving onto scale of cbinom
  out <- ifelse(length(n1)==0,0, #check to make sure there is at least one transition left
                cbinom::dcbinom(x=n2[1],size=size,prob=n1[1],log=T)+log(size)) #for taxon 1 (multiplied to make sense for prob densities)
  for(i in seq_len(length(n1))[-1]){ #for every other taxon i
      out <- out + dcbinom(x=n2[i],size=size-sum(n2[1:i-1]-0.5),prob=n1[i]/(1-sum(n1[1:i-1])),log=T)+log(size)}
  return(out)}
