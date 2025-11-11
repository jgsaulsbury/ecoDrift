#' Plots a set of abundance trajectories as spindles
#'
#' @description
#' A convenient plotting function that takes community composition timeseries in
#' the same format as that used by fitJ and that produced by simNT. Plots a spindle
#' for each species, sorted by maximum abundance in the timeseries.
#'
#'
#' @param occs matrix of the number of observations in each species at each time.
#' One column for each species, one row for each time slice. Time goes from oldest
#' at the bottom to youngest at the top.
#' @param ages vector containing the ages of each time slice, in years, from
#' oldest to youngest.
#' @param plot.ss boolean indicating whether to plot the per-timeslice sample size
#' on the right side of the plot.
#' @param linesevery value indicating the interval at which to draw horizontal gridlines.
#' Default is NA (don't draw them).
#' @param ylab value indicating what to print on the y axis. Default is "Age, years"
#' @param removeyaxis boolean indicating whether to remove the y axis from the plot,
#' allowing the user to plot a custom y axis after the fact. Default is FALSE.
#' @returns plots spindle diagrams depicting the relative abundance of each species
#' through the timeseries.
#' @export
#'
#' @examples
#' J <- 10000
#' tslength <- 500
#' every <- 5
#' nsp <- 4
#' ages <- seq(0,tslength,every)
#' timeseries <- simDrift(startingabs=rep(J/nsp,nsp),ts=ages,ss=1000)
#' plot_spindles(occs=timeseries$simulation,ages=timeseries$times,linesevery=100)
plot_spindles <- function(occs,ages,plot.ss=TRUE,linesevery=NA,ylab="Age, years",removeyaxis=FALSE){
  buffer <- 0.05 #buffer between spindles
  ss <- rowSums(occs)
  occs.prop <-  occs/ss
  occs.prop <- occs.prop[,rev(order(apply(occs.prop,MARGIN=2,FUN=max)))] #resort occs.prop by peak relab
  peak.relabs <- apply(occs.prop, MARGIN=2, FUN=max) #stores peak relab of each sp
  plot(1,type='n',xlim=c(0,sum(peak.relabs)+buffer*(length(peak.relabs)-1)),
       ylim=c(max(ages),min(ages)),xaxt='n',ylab=ylab,xlab="",yaxt=ifelse(removeyaxis,"n","s"))
  if(plot.ss){
    graphics::mtext("Samples:",cex=0.7,line=0,adj=1.03,col='grey')
    graphics::par(las=2)
    graphics::axis(side=4,at=ages,labels=round(ss),cex.axis=0.7,tck=-0.01,hadj=0.5,col.axis='grey',col='grey')
    graphics::par(las=0)}
  graphics::box()
  if(!is.na(linesevery)){#horizontal lines depicting time
    for(t in seq(0,max(ages)+linesevery,linesevery)){
      graphics::lines(c(-1,dim(occs)[2]+1),c(t,t),col="grey90")}}
  for(i in seq(length(peak.relabs))){ #for every sp
    x <- ifelse(i==1,peak.relabs[i]/2,sum(utils::head(peak.relabs,i-1))+peak.relabs[i]/2 + buffer*(i-1))
    graphics::par(las=2)
    graphics::axis(side=1,at=x,labels="",tck=-0.01)
    taxnames <- colnames(occs.prop)[i]
    graphics::axis(side=1,at=x,labels=ifelse(!is.null(taxnames),taxnames,""),tick=FALSE,line=-0.5,cex.axis=0.5 + peak.relabs[i]/max(peak.relabs)*0.3)
    graphics::par(las=0)
    graphics::polygon(x=c(rev(x-occs.prop[,i]/2),x+occs.prop[,i]/2),y=c(rev(ages),ages),col='grey85')
    for(j in seq(length(ages))){ #for every age
      graphics::lines(c(x-occs.prop[j,i]/2,x+occs.prop[j,i]/2),c(ages[j],ages[j]),lwd=1)
    }
  }
}
