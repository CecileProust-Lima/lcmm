print.Diffepoce <- function(x,...){
	if (!inherits(x, "Diffepoce")) stop("use only with \"Diffepoce\" objects")
	cat("Difference in Expected Prognostic Observed Cross-Entropy (EPOCE) estimates", "\n")
	cat(" from the two following joint latent class models:", "\n")
	cat(" \n")
	cl <- x$call.JLCM1
	cl$B <- NULL
	dput(cl)
	cl <- x$call.JLCM1
	cl$B <- NULL
	dput(cl)
	cat(" \n")	
	if(x$new.data==TRUE) {
	cat("Difference in the Mean Prognostic Observed Log-likelihood (MPOL)", "\n")
	cat("      and its 95% tracking interval", "\n")


	}
	else{
	cat("Difference in the Cross-Validated Prognostic Observed Log-likelihood (CVPOL)", "\n")
	cat("      and its 95% tracking interval", "\n")
	}
	cat(" \n")
	print(x$DiffEPOCE,row.names=F,na.print=".")
	cat(" \n")
}
