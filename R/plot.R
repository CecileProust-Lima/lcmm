#' Plot of a fitted model
#' 
#' This function produces different plots (residuals, goodness-of-fit,
#' estimated link functions, estimated baseline risk/survival and posterior
#' probabilities distributions) of a fitted object of class hlme, lcmm,
#' multlcmm or Jointlcmm.
#' 
#' With \code{which="residuals"}, this function provides the marginal residuals
#' against the marginal predictions, the subject-specific residuals against the
#' subject-specific predictions, a normal QQ-plot with confidence bands for the
#' marginal residuals and a normal QQ-plot with confidence bands for the
#' subject-specific residuals.
#' 
#' With \code{which="postprob"}, the function provides the histograms of the
#' posterior class-membership probabilities stemmed from a \code{Jointlcmm},
#' \code{lcmm}, \code{hlme} or \code{multlcmm} object.
#' 
#' With \code{which="link"} or \code{which="linkfunction"}, the function
#' displays the estimated transformation(s) specified in the option \code{link}
#' of \code{lcmm} and \code{multlcmm} functions. It corresponds to the
#' (non)linear parameterized link estimated between the oberved longitudinal
#' outcome and the underlying latent process.
#' 
#' With \code{which="fit"}, the function provides the class-specific weighted
#' marginal and subject-specific mean predicted trajectories with time and the
#' class-specific weighted mean observed trajectories and their 95\% confidence
#' bounds. The predicted and observed class-specific values are weighted means
#' within each time interval; For each observation or prediction (in the
#' transformed scale if appropriate), the weights are the class-specific
#' (posterior with subject-specific or marginal otherwise) probabilities to
#' belong to the latent class.
#' 
#' With \code{which="baselinerisk"} or \code{which="hazard"}, the function
#' displays the estimated baseline risk functions for the time-to-event of
#' interest in each latent class.
#' 
#' With \code{which="survival"}, the function displays the estimated event-free
#' probabilities (survival functions) for the time-to-event of interest in each
#' latent class.
#' 
#' @param x an object inheriting from classes \code{hlme}, \code{lcmm},
#' \code{multlcmm} or \code{Jointlcmm}, representing respectively a fitted
#' latent class linear mixed model, a more general latent class mixed model or
#' a joint latent class model
#' @param which a character string indicating the type of plot to produce. For
#' \code{hlme} objects, are available "residuals", "postprob","fit". For
#' \code{lcmm} and \code{multlcmm} objects, are available "residuals",
#' "postprob", "link", "linkfunction", "fit".  For \code{Jointlcmm} objects,
#' are avaiable "residuals", "postprob", "link", "linkfunction", "fit",
#' "hazard", "baselinerisk", "survival". Default to "residuals"
#' @param var.time for \code{which="fit"} only, a character string containing
#' the name of the variable that corresponds to time in the longitudinal model.
#' @param break.times for \code{which="fit"} only, either a numeric vector
#' containing the cuts-off defining the time-intervals or an integer giving the
#' number of cut-offs. In the latter case, the cut-offs are placed at the
#' quantiles of the observed times distribution.
#' @param marg for \code{which="fit"} only, a logical indicating the type of
#' prediction. If \code{marg=TRUE} (the default), the marginal predictions are
#' provided. If \code{marg=FALSE}, the subject-specific predictions are
#' provided.
#' @param outcome for \code{which="fit"} and \code{multlcmm} objects only, the
#' outcome to consider.
#' @param event for \code{which="baselinerisk"} or \code{which="hazard"} only,
#' an integer corresponding to the numeric code (in the indicator variable) of
#' the event for which the baseline risk functions are to be plotted. By
#' default, the first event is considered.
#' @param subset for \code{which="fit"} only, a subset of the data used to
#' estimate the model, defining the data on which the fit is evaluated. By
#' default, all the data are used.
#' @param \dots other parameters to be passed through to plotting functions.
#' This includes graphical parameters described in par function and further arguments
#' legend (character or expression
#' to appear in the legend. If no legend should be added, \code{"legend"}
#' should be NULL. ),
#' legend.loc (keyword for the position of the legend from the list
#' \code{"bottomright"}, \code{"bottom"}, \code{"bottomleft"}, 
#' \code{"left"}, \code{"topleft"},\code{"top"}, \code{"topright"},
#' \code{"right"} and \code{"center"}. By default, the legend is located in
#' the top left of the plot. ) and
#' add (logical indicating if the curves
#' should be added to an existing plot. Default to FALSE.).
#' @param shades logical indicating if confidence intervals should be
#' represented with shades. Default to FALSE, confidence intervals are
#' represented as dotted lines.
#' 
#' @author Cecile Proust-Lima, Viviane Philipps and Benoit Liquet
#' @seealso \code{\link{hlme}}, \code{\link{lcmm}}, \code{\link{multlcmm}},
#' \code{\link{Jointlcmm}}
#' @examples
#' 
#' 
#' ###################### fit, residuals and postprob 
#' 
#' # estimation of the model
#' m<-lcmm(Y~Time*X1,mixture=~Time,random=~Time,classmb=~X2+X3,
#' subject='ID',ng=2,data=data_hlme,B=c(0.41,0.55,-0.18,-0.41,
#' -14.26,-0.34,1.33,13.51,24.65,2.98,1.18,26.26,0.97))  
#' 
#' # fit
#' plot(m,which="fit",marg=FALSE,var.time="Time",bty="n")
#' # residuals plot
#' plot(m)
#' # postprob plot
#' plot(m,which="postprob") 
#' 
#' 
#' ###################### fit, linkfunctions
#' 
#' #### Estimation of homogeneous mixed models with different assumed link
#' #### functions, a quadratic mean trajectory for the latent process with 
#' #### independent random intercept, slope and quadratic slope
#' #### (comparison of linear, Beta and 3 and 5 splines link functions)
#' \dontrun{
#' 
#' # linear link function
#' m10<-lcmm(Ydep2~Time+I(Time^2),random=~Time+I(Time^2),subject='ID',ng=1,
#'           data=data_lcmm,link="linear",
#'           B=c(-0.7454, -0.2031,  0.2715,  0.2916 , 0.6114, -0.0064,  0.0545,
#'               0.0128, 25.3795, 2.2371))
#'             
#' # Beta link function
#' m11<-lcmm(Ydep2~Time+I(Time^2),random=~Time+I(Time^2),subject='ID',ng=1,
#'           data=data_lcmm,link="beta",B=c(-0.9109, -0.0831,  0.5194,  0.1910 ,
#'           0.8984, -0.0179, -0.0636,  0.0045,  0.5514, -0.7692,  0.7037,  0.0899))
#'           
#' # fit 
#' par(mfrow=c(2,1),mar=c(4,4,1,1))
#' plot(m11,which="fit",var.time="Time",bty="l",ylim=c(-3,0))
#' plot(m11,which="fit",var.time="Time",marg=FALSE,bty="l",ylim=c(-3,0))
#' 
#' # I-splines with 3 equidistant nodes
#' m12<-lcmm(Ydep2~Time+I(Time^2),random=~Time+I(Time^2),subject='ID',ng=1,
#'           data=data_lcmm,link="3-equi-splines",B=c(-0.9272, -0.0753 , 0.5304, 
#'           0.1950,  0.9260, -0.0204, -0.0739 , 0.0059, -7.8369,  0.9228 ,-1.4689,
#'           2.0396,  1.8102))
#' 
#' # I-splines with 5 nodes, and interior nodes entered manually
#' m13<-lcmm(Ydep2~Time+I(Time^2),random=~Time+I(Time^2),subject='ID',ng=1,
#'           data=data_lcmm,link="5-manual-splines",intnodes=c(10,20,25),
#'           B=c(-0.9315, -0.0739 , 0.5254 , 0.1933,  0.9418, -0.0206, -0.0776,
#'           0.0064, -7.8645, 0.7470,  1.2080,  1.5537 , 1.7558 , 1.3386 , 1.0982))
#' 
#' # Plot of estimated different link functions:
#' # (applicable for models that only differ in the "link function" used. 
#' # Otherwise, the latent process scale is different and a rescaling
#' # is necessary)
#' plot(m10,which="linkfunction",bty="l")
#' plot(m11,which="linkfunction",bty="l",add=TRUE,col=2)
#' plot(m12,which="linkfunction",bty="l",add=TRUE,col=3)
#' plot(m13,which="linkfunction",bty="l",add=TRUE,col=4)
#' legend("topleft",legend=c("linear","beta","3-Isplines","5-Isplines"),
#' col=1:4,lty=1,bty='n')
#' }
#' 
#' 
#' ###################### fit, baselinerisk and survival
#' \dontrun{
#' #### estimation with 3 latent classes (ng=3) - see Jointlcmm 
#' #### help for details on the model
#' m3 <- Jointlcmm(fixed= Ydep1~Time*X1,mixture=~Time,random=~Time,
#' classmb=~X3,subject='ID',survival = Surv(Tevent,Event)~ X1+mixture(X2),
#' hazard="3-quant-splines",hazardtype="PH",ng=3,data=data_lcmm,
#' B=c(0.7576, 0.4095, -0.8232, -0.2737, 0, 0, 0, 0.2838, -0.6338, 
#' 2.6324, 5.3963, -0.0273, 1.3979, 0.8168, -15.041, 10.164, 10.2394, 
#' 11.5109, -2.6219, -0.4553, -0.6055, 1.473, -0.0383, 0.8512, 0.0389, 
#' 0.2624, 1.4982))
#' 
#' # fit
#' plot(m3,which="fit",var.time="Time",bty="l")
#' plot(m3,which="fit",var.time="Time",marg=FALSE,bty="l",ylim=c(0,15))
#' 
#' 
#' # Class-specific predicted baseline risk & survival functions in the 
#' # 3-class model retained (for the reference value of the covariates) 
#' plot(m3,which="baselinerisk",bty="l")
#' plot(m3,which="baselinerisk",ylim=c(0,5),bty="l")
#' plot(m3,which="survival",bty="l")
#' }
#' 
#' 
#' @name plot
#' 
#' @export
#' 
#' 
#' 
#' 
#' 
#' 
plot.hlme <- function(x,which="residuals",var.time,break.times,marg,subset,shades,...)
    {
        if(missing(x)) stop("The model should be specified")
        if(!inherits(x,"hlme")) stop("Use with 'hlme' objects only")

        if(!(which %in% c("residuals","postprob","fit"))) stop(paste("Argument 'which' should be one of:",paste(c("residuals","postprob","fit"),collapse=" ")))

        if(which=="residuals") .plotresid(x,...)

        if(which=="postprob") .plotpostprob(x,...)

        if(which=="fit")
            {
                if(missing(var.time))
                {
                    if(length(x$var.time))
                    {
                        var.time <- x$var.time
                    }
                    else
                    {
                        stop("Argument var.time should be specified")
                    }
                }
                if(missing(break.times)) break.times <- NULL
                if(missing(marg)) marg <- TRUE
                if(missing(subset)) subset <- NULL
                if(missing(shades)) shades <- FALSE
                #legend.loc?add?
                .plotfit(x,var.time=var.time,break.times,outcome=1,marg=marg,subset=subset,shades=shades,...)

            }

    }






#' @rdname plot
#' @export
plot.lcmm <- function(x,which="residuals",var.time,break.times,marg,subset,shades,...)
    {
        if(missing(x)) stop("The model should be specified")
        if(!inherits(x,"lcmm")) stop("Use with 'lcmm' objects only")
        
        if(!(which %in% c("residuals","postprob","link","linkfunction","fit"))) stop(paste("Argument 'which' should be one of:",paste(c("residuals","postprob","link","linkfunction","fit"),collapse=" ")))

        if(which=="residuals") .plotresid(x,...)

        if(which=="postprob") .plotpostprob(x,...)
           
        if(which %in% c("link","linkfunction")) .plotlinkfunction(x,...)
        

        if(which=="fit")
            {
                if(missing(var.time))
                {
                    if(length(x$var.time))
                    {
                        var.time <- x$var.time
                    }
                    else
                    {
                        stop("Argument var.time should be specified")
                    }
                }
                if(missing(break.times)) break.times <- NULL
                if(missing(marg)) marg <- TRUE
                if(missing(subset)) subset <- NULL
                if(missing(shades)) shades <- FALSE

                .plotfit(x,var.time=var.time,break.times,outcome=1,marg=marg,subset=subset,shades=shades,...)
            }
    }



#' @rdname plot
#' @export
plot.multlcmm <- function(x,which="residuals",var.time,break.times,marg,outcome,subset,shades,...)
    {
        if(missing(x)) stop("The model should be specified")
        if(!inherits(x,"multlcmm")) stop("Use with 'multlcmm' objects only")
        
        if(!(which %in% c("residuals","postprob","link","linkfunction","fit"))) stop(paste("Argument 'which' should be one of:",paste(c("residuals","postprob","link","linkfunction","fit"),collapse=" ")))

        if(which=="residuals") .plotresid(x,...)

        if(which=="postprob") .plotpostprob(x,...)

        if(which %in% c("link","linkfunction")) .plotlinkfunctionmult(x,...)
           
        if(which=="fit")
            {
                if(missing(var.time))
                {
                    if(length(x$var.time))
                    {
                        var.time <- x$var.time
                    }
                    else
                    {
                        stop("Argument var.time should be specified")
                    }
                }
                if(missing(break.times)) break.times <- NULL
                if(missing(marg)) marg <- TRUE
                if(missing(subset)) subset <- NULL
                if(missing(outcome)) outcome <- 1
                if(missing(shades)) shades <- FALSE
                
                .plotfit(x,var.time=var.time,break.times,outcome=outcome,marg=marg,subset=subset,shades=shades,...)
            }

    }




#' @rdname plot
#' @export
plot.Jointlcmm <- function(x,which="residuals",var.time,break.times,marg,event,subset,shades,...)
    {
        if(missing(x)) stop("The model should be specified")
        if(!inherits(x,"Jointlcmm")) stop("Use with 'Jointlcmm' objects only")
        
        if(!(which %in% c("residuals","postprob","link","linkfunction","fit","hazard","baselinerisk","survival"))) stop(paste("Argument 'which' should be one of:",paste(c("residuals","postprob","link","linkfunction","fit","hazard","baselinerisk","survival"),collapse=" ")))

        if(which=="residuals") .plotresid(x,...)

        if(which=="postprob") .plotpostprob(x,...)

        if(which %in% c("link","linkfunction")) .plotlinkfunction(x,...)

        if(which %in% c("hazard","baselinerisk"))
            {
                if(missing(event)) event <- 1
                .plotbaselinerisk(x,event=event,...) #legend.loc legend add
            }
           
        if(which=="survival") .plotsurvival(x,...)

        if(which=="fit")
            {
                if(missing(var.time))
                {
                    if(length(x$var.time))
                    {
                        var.time <- x$var.time
                    }
                    else
                    {
                        stop("Argument var.time should be specified")
                    }
                }
                if(missing(break.times)) break.times <- NULL
                if(missing(marg)) marg <- TRUE
                if(missing(subset)) subset <- NULL
                if(missing(shades)) shades <- FALSE
                
                .plotfit(x,var.time=var.time,break.times,outcome=1,marg=marg,subset=subset,shades=shades,...)
            }

    }




#' @rdname plot
#' @export
plot.mpjlcmm <- function(x,which,event,...)
{
    if(missing(x)) stop("The model should be specified")
    if(!inherits(x,"mpjlcmm")) stop("Use with 'mpjlcmm' objects only")
    
    if(which %in% c("residuals","link","linkfunction","fit")) stop(paste("For ",paste(c("residuals","link","linkfunction","fit"),collapse=" "),", please see the 'update' function",sep=""))
    
    if(!(which %in% c("hazard","baselinerisk","survival"))) stop("Argument 'which' should be one of hazard, baselinerisk or survival.")

    xx <- list(N=c(rep(0,9+x$nbevt)), conv=x$conv, predSurv=x$predSurv, ng=x$ng)
    class(xx) <- "Jointlcmm"
    
    plot.Jointlcmm(xx,which=which,event=event,...)
}





#' @rdname plot
#' @export
plot.externSurv <- function(x,which="hazard",event,...)
{
  if(missing(x)) stop("The model should be specified")
  if(!inherits(x,"externSurv")) stop("Use with 'externSurv' objects only")
  
  if(!(which %in% c("postprob","hazard","baselinerisk","survival"))) stop(paste("Argument 'which' should be one of:",paste(c("postprob","hazard","baselinerisk","survival"),collapse=" ")))
  
  if(which=="postprob") .plotpostprob(x,...)
  
  if(which %in% c("hazard","baselinerisk"))
  {
    if(missing(event)) event <- 1
    .plotbaselinerisk(x,event=event,...) #legend.loc legend add
  }
  
  if(which=="survival") .plotsurvival(x,...)
  
}



#' @rdname plot
#' @export
plot.externX <- function(x,which="postprob",event,...)
{
  if(missing(x)) stop("The model should be specified")
  if(!inherits(x,"externX")) stop("Use with 'externX' objects only")
  
  if(!(which %in% c("postprob"))) stop(paste("Argument 'which' should be one of:",paste(c("postprob"),collapse=" ")))
  
  if(which=="postprob") .plotpostprob(x,...)
  
}

