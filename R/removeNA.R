## Function : removeNA
## Args :
##  - form : a list of formula
##  - data : a data frame containing all variables mentionned in form
## Returns : a list containing
##  - newdata : a data frame created from data, containing only the variables used in the formula, and without any missing values
##  - Y : the names of the response variables
##  - YlinesNA : the lines removed from each response variables
##  - nmes : the number of measures for each response variable
##  - X : the names of the covariates
##  - XlinesNA : the lines removed from data because of missing values in the covariates
##
## Caution : if several response variables are included in 'form', then the resulting data frame will contain one column (named outcome) grouping all these response variables. The irst lines will refer to the first response variables, a second block to the second, etc. All covariates will be repeated for each response variable.
## Example :
##  d <- data.frame(y1=c(3,6,4,2), y2=c(6.5,NA,9.5,4.5), x=c(0,0,1,0))
##  removeNA(form=list(y1+y2~x), data=d) gives
##   $newdata
##     outcome x
##  1      3.0 0   | first block referring to y1
##  2      6.0 0   |
##  3      4.0 1   |
##  4      2.0 0   |
##  11     6.5 0     | second block referring to y2
##  31     9.5 1     |
##  41     4.5 0     |

removeNA <- function(form, data)
{
    nform <- length(form)
    vars <- vector("list", length=nform) #covariates names
    nas <- vector("list", length=nform) # NA in covariates
    outcomes <- NULL # outcomes names
    naoutcomes <- list() # NA in each outcome
    term <- vector("list", length=nform)

    ## loop on the formulas : keep variables and lines with NA
    for(k in 1:nform)
    {
        ## covariates
        termsk <- delete.response(terms(form[[k]]))
        dk <- model.frame(termsk, data=data)
        vars[[k]] <- all.vars(termsk)
        nas[[k]] <- na.action(dk)
        term[[k]] <- terms(dk)

        ## outcomes
        if(attr(terms(form[[k]]), "response")==1)
        {
            y <- all.vars(form[[k]][[2]])
            outcomes <- c(outcomes, y)
            
            ny <- length(y)
            for(m in 1:ny)
            {
                naoutcomes <- c(naoutcomes, list(which(is.na(data[,y[m]]))))
            }           
        }
    }
    
    ## create newdata
    allvars <- unique(unlist(vars))
    linesNA <- unique(unlist(nas))
    
    if(is.null(outcomes))
    {
        ## keep same form as data
        
        if(is.null(linesNA))
        {
            newdata <- data[, allvars, drop=FALSE]
        }
        else
        {
            newdata <- data[-linesNA, allvars, drop=FALSE]
        }
        nmes <- nrow(newdata)
    }
    else
    {
        ## rbind the data, with one block per outcome
        ## a block has nrow(data)-length(linesNA) rows
        
        newdata <- NULL
        nmes <- rep(NA,length(outcomes))
        for(k in 1:length(outcomes))
        {
            nak <- unique(c(linesNA, naoutcomes[[k]]))
            if(!length(nak))
            {
                newdk <- data[, c(outcomes[k],allvars), drop=FALSE]
            }
            else
            {
                newdk <- data[-nak, c(outcomes[k],allvars), drop=FALSE]
            }
            colnames(newdk) <- c("outcome", allvars)
            nmes[k] <- nrow(newdk)

            newdata <- rbind(newdata, newdk)
        }
    }
    
    return(list(newdata=newdata, Y=outcomes, YlinesNA=naoutcomes, nmes=nmes, X=allvars, XlinesNA=linesNA, terms=term))    
}

