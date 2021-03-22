#' Plot DepthGram
#'
#' Function to plot the 3 DepthGram representations from the output of the depthGram function
#'
#' @param DG output of the depthGram function.
#' @param limits if TRUE the empirical limits for outlier detection are drawn. Default is FALSE.
#' @param ids labels for data points.
#' @param print if TRUE the graphical output is optimized for printed version. Default is FALSE.
#' @param plotly if TRUE the graphical output is displayed as an interactive plotly object. Default is FALSE.
#' @param plot.title character with the main title for plot.
#' @param shorten If labels must be shorten to 15 characters. Default is TRUE.
#' @param col colors for the plot. Default is NULL.
#' @param pch point shape.
#' @param sp point size.
#' @param st label size.
#' @param sa axis title sizes.
#' @param text.labels labels for the individuals. text.labels is overridden if limits=T, for which only outliers labels are shown.
#'
#' @return A list with the following items:
#' \itemize{
#' \item p - list with all the depthGram plots.
#' \item out - outliers detected.
#' \item colors - used colors for plotting.
#' }
#'
#' @references Aleman-Gomez, Y., Arribas-Gil, A., Desco, M. Elias-Fernandez, A., and Romo, J. (2021). "Depthgram: Visualizing Outliers in High Dimensional Functional Data with application to Task fMRI data exploration". <arXiv:2103.08874>
#'
#' @export
#'
#' @import gridExtra ggplot2
#' @importFrom plotly subplot layout
#' @importFrom grDevices hcl
#'
#' @examples
#'
#' N = 2e2
#' P = 1e3
#' grid = seq( 0, 1, length.out = P )
#' Cov = exp_cov_function( grid, alpha = 0.3, beta = 0.4 )
#'
#' Data <- list()
#' Data[[1]] <- generate_gauss_fdata( N, centerline = sin( 2 * pi * grid ), Cov = Cov )
#' Data[[2]] <- generate_gauss_fdata( N, centerline = sin( 2 * pi * grid ), Cov = Cov )
#' names <- paste0("id_", 1:nrow(Data[[1]]))
#' DG <- depthGram(Data, marg.out=TRUE, ids=names)
#' depthGramPlot(DG)
depthGramPlot <- function(DG, limits=F, ids=NULL, print=F, plotly=F, plot.title="", shorten=T, col = NULL, pch=19, sp=2, st=4, sa=10, text.labels=""){

  n=length(DG$mei.mbd.d)

  type=c("Dimensions DepthGram","Time DepthGram","Time/Correlation DepthGram")

  if(is.null(ids)){
    ids=as.character(1:n)
  }

  DG=data.frame(ID=rep(ids,3), mei.mbd=c(DG$mei.mbd.d,DG$mei.mbd.t,DG$mei.mbd.t2), mbd.mei=c(DG$mbd.mei.d, DG$mbd.mei.t, DG$mbd.mei.t2),
                type=rep(type,each=n))


  if(is.null(col)){
    hues = seq(15, 375, length=n+1)
    color <- grDevices::hcl(h=hues, l=65, c=100)[1:n]
  }else{
    color=col
  }

  out=NULL
  if(limits){
    P2<-function(x,n){a0=2/n;a2=-n/(2*(n-1));return(a0+x+a2*x^2)}
    meis=seq(0,1, ,n)
    DG$meis=rep(meis,3)
    DG$par=P2(DG$meis,n)
    distp=DG$mbd.mei-P2(1-DG$mei.mbd,n)
    q3=quantile(distp,0.75)
    q1=quantile(distp,0.25)
    DG$par2=DG$par+q3+1.5*(q3-q1)
    out<-unique(which(distp>q3+1.5*(q3-q1))%%n)
    pch=rep(1,n) #empty circle
    pch[out]=19 #solid circle
    text.labels=rep("",n)
    if (!plotly){
      if (shorten){ #shorten text labels
        text.labels[out]=sapply(ids[out],function(x) substr(x,1,min(15,nchar(x))))
      }else{
        text.labels[out]=ids[out]
      }
    }
    if(is.null(col)){
      hues = seq(15, 375, length=length(out)+1)
      color.out <- hcl(h=hues, l=65, c=100)[1:length(out)]
      color<-rep(8,n)
      color[out]<-color.out
    }
  }

  if(print){
    sp=3
    st=5
    sa=14
  }

  plots=list()

  for (i in 1:3){

    dat=DG[which(DG$type==type[i]),]

    plots[[i]]<- ggplot2::ggplot(dat, aes(x=1-.data$mei.mbd, y=.data$mbd.mei)) +
      ggplot2::geom_point(aes(x=1-.data$mei.mbd, y=.data$mbd.mei, group=.data$ID),color=color,size=sp,shape=pch) +
      ggplot2::geom_text(aes(x=1-.data$mei.mbd, y=.data$mbd.mei),label=text.labels,color=color,hjust=-0.15, vjust=-0.15,size=st)+
      ggplot2::xlim(c(0,1.005))+
      ggplot2::ylim(c(0,0.525)) +
      ggplot2::theme_minimal()+
      ggplot2::theme(axis.title= element_text(size=sa),axis.text= element_text(size=sa-2),title=element_text(size=sa))

    if(limits==TRUE){
      plots[[i]]<- plots[[i]] + ggplot2::geom_line(aes(x=.data$meis, y=.data$par),col=1,na.rm=T)+
        ggplot2::geom_line(aes(x=.data$meis, y=.data$par2),col=1,lty=2,na.rm=T)
    }

    if (i==1){pt<-plot.title}else{pt<-""}

    if (plotly){
      plots[[i]]<-plots[[i]]+ggplot2::facet_wrap(~.data$type)+ggplot2::ggtitle(pt)
    }else{
      if(i<=2){
        plots[[i]]<-plots[[i]]+
          ggplot2::labs(title=pt,subtitle=type[i],
                        x=expression("1-MEI("~paste(bold('MBD'[d])~")")),
                        y=expression("MBD("~paste(bold('MEI'[d])~")")))
      }else{
        plots[[i]]<-plots[[i]] +
          ggplot2::labs(title=pt,subtitle=type[i], x=expression("1-MEI("~paste(bold(widetilde(MBD)[t])~")")),
                                      y=expression("MBD("~paste(bold(widetilde(MEI)[t])~")")))
      }
    }

  }

  if(plotly){
    p<-plotly::subplot(plots,shareY = TRUE,shareX =TRUE)%>%
       plotly::layout(title = plot.title, yaxis = list(title="MBD(MEI)"), xaxis = list(title="1-MEI(MBD)"))
    print(p)
  }else{
    p<-do.call(gridExtra::grid.arrange, c(plots,nrow=1))
  }

  return(list(p=list(dimDG = plots[[1]], timeDG = plots[[2]], corrDG = plots[[3]], fullDG = p), out=out, color=color))
}
