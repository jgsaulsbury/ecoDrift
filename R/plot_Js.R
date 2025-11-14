#' Plot fitted J values for every transition in a single timeseries
#'
#' @description
#' Use to visualize how the rate of turnover changes across a community change
#' timeseries.
#'
#'
#' @param occs matrix of the number of observations in each species at each time.
#' One column for each species, one row for each time slice. Time goes from oldest
#' at the bottom to youngest at the top.
#' @param ages vector containing the ages of each time slice, in years, from
#' oldest to youngest.
#' @param xlim vector of length 2 giving the left and right bounds of the x axis.
#' @param linesevery value indicating the interval at which to draw horizontal gridlines.
#' Default is NA (don't draw them).
#' @param generationtime passed to fitJ().
#' @param searchinterval passed to fitJ()
#' @param ... additional plotting parameters passed to plot().
#'
#' @returns plots rate through time.
#' @export
#'
#' @examples
#' #simulate a timeseries with an increase in J halfway through
#' J1 <- 1000 #first half of TS
#' J2 <- 10000 #second half of TS
#' tslength <- 500
#' every <- 5
#' nsp <- 8
#' ss <- 1000000
#' ages <- seq(0,tslength,every)
#' timeseries1 <- simDrift(startingabs=rep(J1/nsp,nsp),ts=ages[1:(length(ages)/2 + 1)],ss=ss)$simulation #first half
#' timeseries2 <- simDrift(startingabs=timeseries1[1,]*J2/ss,ts=ages[1:(length(ages)/2 + 1)],ss=ss)$simulation #second half
#' timeseries <- rbind(timeseries2,timeseries1[-1,])
#' par(mfrow=c(1,2))
#' plot_spindles(timeseries,ages,plot.ss=FALSE)
#' lines(c(-1,10),c(tslength/2,tslength/2),lty='dashed',lwd=1) #plot point where J switches
#' plot_Js(timeseries,ages) #fit rate through time
#' lines(c(1E-10,1E-1),c(tslength/2,tslength/2),lty='dashed',lwd=1)
plot_Js <- function(occs,ages,xlim=NULL,linesevery=NA,sampled=TRUE,generationtime=1,searchinterval=c(1,9),ylab="Age, years",removeyaxis=FALSE,...){
  if(is.list(occs)){
    occs <- as.matrix(occs)}
  xages <- (head(ages,-1) + tail(ages,-1))/2 #midpoint of each transition
  Jhats <- c() #store ML values for each transition, oldest to youngest
  JLBs <- c() #lower bounds on J
  JUBs <- c() #upper bounds on J
  for(i in seq(length(ages)-1)){ #for every transition
    transition <- occs[(length(ages)-i):(length(ages)-i+1),]
    fit <- fitJ(occs=transition,ages=ages[i:(i+1)],sampled=sampled,generationtime=generationtime,CI=TRUE,searchinterval=searchinterval)
    Jhats <- c(Jhats,fit$J)
    JLBs <- c(JLBs,fit$CI[1])
    JUBs <- c(JUBs,fit$CI[2])}
  if(is.null(xlim)){
    xlim <- c(min(1/JUBs),max(1/JLBs))}
  if(xlim[1]>xlim[2]){xlim <- rev(xlim)}
  plot(1,type="n",xlim=xlim,ylim=c(ages[1],tail(ages,1)),log='x',xlab="1/J",ylab=ylab,xaxt='n',yaxt=ifelse(removeyaxis,"n","s"),...)
  if(!is.na(linesevery)){#horizontal lines depicting time
    for(t in seq(0,max(ages)+linesevery,linesevery)){
      graphics::lines(c(xlim[1]/10,xlim[2]*10),c(t,t),col="grey90")}}
  points(1/Jhats,xages)
  for(i in seq(length(xages))){ #for every transition
    graphics::lines(c(1/JLBs[i],1/JUBs[i]),c(xages[i],xages[i]))} #draw error bars
  axis(side=1,at=10^seq(floor(log10(xlim[1])),ceiling(log10(xlim[2])))) #tick marks
  for(i in seq(floor(log10(xlim[1])),floor(log10(xlim[2])))){
    axis(side=1,at=c(seq(10^i,10^(i+1),10^i)),labels = FALSE,tck=-0.01)
  }
}
