#' Probability of a single transition under neutral theory with migration
#'
#' @description
#' Gives the log-likelihood of a single observed transition in community composition
#' under Hubbell's neutral theory, for a given ratio of J (community size) to
#' t (# generations) and a given migration rate m from a given metacommunity.
#'
#' @details
#' Used by xxprobm. Allows for incomplete sampling; if nothing is given for ss,
#' assumes n1 and n2 represent true relative abundances. m is the chance that a new
#' individual is drawn from the static metacommunity rather than the previous
#' generation of the local community.
#'
#' @param n1 vector of relative abundances at start of transition.
#' @param n2 vector of relative abundances at end of transition.
#' @param nmeta vector of relative abundances in the metacommunity. Should sum to 1;
#' if not, nmeta will be renormalized to sum to 1.
#' @param J community size
#' @param t time (# generations between n1 and n2)
#' @param m migration rate.
#' @param ss single value or vector of length 2 giving number of samples at
#' times for n1 and n2. If a single value given, assumes it applies to both n1 and n2.
#'
#' @importFrom cbinom dcbinom
#'
#' @returns loglik value.
#'
#' @examples
#' ecoDrift:::xprobm(n1=c(0.1,0.1,0.8),n2=c(0.3,0.3,0.4),
#'   nmeta=c(0.4,0.4,0.2),J=1000,t=20,m=0.05,ss=10000)
#' #4.494708
#' # shorter timescale:
#' ecoDrift:::xprobm(n1=c(0.1,0.1,0.8),n2=c(0.3,0.3,0.4),
#'   nmeta=c(0.4,0.4,0.2),J=1000,t=10,m=0.05,ss=10000)
#' #-2.949927
#' #smaller sample size:
#' ecoDrift:::xprobm(n1=c(0.1,0.1,0.8),n2=c(0.3,0.3,0.4),
#'   nmeta=c(0.4,0.4,0.2),J=1000,t=10,m=0.05,ss=100)
#' #1.501469
#'
xprobm <- function(n1,n2,nmeta,J,t,m,ss=NA){
  #error handling
  tol <- 1E-7
  if(sum(n1)<=0|sum(n2)<=0|sum(n1)>1+tol|sum(n2)>1+tol){
    stop("Sums of n1 and n2 must each be greater than 0 and no greater than 1")}
  if(length(n1)!=length(n2)){stop("n1 and n2 must have the same length")}
  if(!J/t>0){stop("Jt must be positive")}
  if(m<0 | m>1){stop("m must be between 0 and 1")}
  if(any(nmeta<0) | sum(nmeta)<=0){stop("nmeta must be all nonnegative and sum to >0")}
  #body
  if(sum(nmeta) != 1){nmeta <- nmeta/sum(nmeta)} #fix nmeta if needed
  order <- rev(order(n1)) #prepping n1 and n2
  n1 <- n1[order]
  n2 <- n2[order]
  nmeta <- nmeta[order]
  while(sum(n1) > 1-tol | sum(n2) > 1-tol){ #remove last species if sum is 1, not needed for calculation
    n1 <- n1[-length(n1)]
    n2 <- n2[-length(n2)]
    nmeta <- nmeta[-length(nmeta)]}
  n1 <- n1[!n2==0];n2 <- n2[!n2==0];nmeta <- nmeta[!nmeta==0] #remove indices where n2 is zero (zeros in n1 removed in previous step)
  #prep migration
  weight.local <- (1-m)**t
  prob <- n1*weight.local + nmeta*(1-weight.local)
  size <- J/t*weight.local + J*(1-(1-m)**2)
  #incorporate sample size
  if(!any(is.na(ss))){ #if ss given...
    if(length(ss)==1){ss <- rep(ss,2)} #duplicate if a single value given
    sizes <- c(ss,size)
    size <- min(sizes) / sum(min(sizes)/sizes)}
  n2 <- n2*size + 0.5 #moving onto scale of cbinom
  #start with species 1
  out <- ifelse(length(n1)==0,0,
    cbinom::dcbinom(x=n2[1],size=size,prob=prob[1],log=T)+log(size)) #for taxon 1
  for(i in seq_len(length(n1))[-1]){ #for every other taxon i
    out <- out + dcbinom(x=n2[i],size=size-sum(n2[1:i-1]-0.5),prob=prob[i]/(1-sum(prob[1:i-1])),log=T)+log(size)
    }
  return(out)}
