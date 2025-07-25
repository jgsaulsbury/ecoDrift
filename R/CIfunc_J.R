#Helper function for getting confidence intervals on fits of J.
#Used by fitJ. Similar to xxprob, but optimizes a function that has its optimum some distance
#below the maximum likelihood. That distance is 1.92 for a 95% confidence interval with 1 free parameter
#(based on a chi-square table). confidence indicates degree of confidence desired.
#Default is 95% (0.95). Returns the value of the helper function at the given logJ.
CIfunc_J <- function(log10J,occs,ages,ML,confidence=0.95,sampled=TRUE,generationtime=1){
  diff <- stats::qchisq(p=confidence,d=1)/2
  out <- xxprob(log10J,occs,ages,sampled,generationtime)
  return(abs(ML-diff-out))}
