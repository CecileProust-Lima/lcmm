.plotpostprob <- function(x,...)
{
    if(x$ng==1)
        {
            cat("plot.postprob can only be used when ng > 1  \n")
        }
    else
        {
            if(x$conv==1|x$conv==2|x$conv==3) 
                {
                    dots <- list(...)
                    
                    dots <- dots[setdiff(names(dots),c("x","y","log"))] #ce qu'on ne veut pas changer

                    if(length(list(...)$main)) 
                        {
                            title1 <- as.character(eval(match.call()$main))
                            dots <- dots[setdiff(names(dots),"main")]
                            title1 <- rep(title1,length.out=x$ng)
                        }
                    else title1 <- paste("distribution of posterior probabilities in class",1:x$ng)
                    
                    if(length(list(...)$xlab)) 
                        {
                            xlab1 <- as.character(eval(match.call()$xlab))
                            dots <- dots[setdiff(names(dots),"xlab")]
                            xlab1 <- rep(xlab1,length.out=x$ng)
                        }
                    else xlab1 <- paste("class",1:x$ng)
                    
                    if(length(list(...)$probability)) 
                        {
                            prob1 <- as.character(eval(match.call()$probability))
                            dots <- dots[setdiff(names(dots),"probability")]
                        }
                    else prob1 <- TRUE

                    if(length(list(...)$col)) 
                        {
                            color <- as.character(eval(match.call()$col))
                            dots <- dots[setdiff(names(dots),"col")]
                            color <- rep(color,length.out=x$ng)
                        }
                    else color <- NULL
                    
                    if(length(list(...)$border)) 
                        {
                            border1 <- as.character(eval(match.call()$border))
                            dots <- dots[setdiff(names(dots),"border")]
                            border1 <- rep(border1,length.out=x$ng)
                        }
                    else border1 <- NULL

                    for(i in 1:x$ng)
                        {
                            do.call("hist",c(dots,list(x=x$pprob[,i+2],probability=prob1,xlab=xlab1[i],main=title1[i],col=color[i],border=border1[i])))   #va pas si col=1:ng
                        }
                }
            else
                {
                    cat("Output can not be produced since the program stopped abnormally.\n")
                }
        }
}


