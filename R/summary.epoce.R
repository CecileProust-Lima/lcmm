summary.epoce <- function(object,...){
	x <- object
	if (!inherits(x, "epoce")) stop("use only with \"epoce\" objects")
	cat("Expected Prognostic Observed Cross-Entropy (EPOCE) of the joint latent class model:", "\n")
	cat(" \n")
	cl <- x$call.Jointlcmm
	cl$B <- NULL
	dput(cl)
	cat(" \n")	
	if(x$new.data==TRUE) {
	cat("EPOCE estimator on new data: Mean Prognostic Observed Log-likelihood (MPOL)", "\n")
	}
	else{
	cat("EPOCE estimators on data used for estimation:", "\n")
	cat("    Mean Prognostic Observed Log-likelihood (MPOL)", "\n")
	cat("    and Cross-validated Prognostic Observed Log-likelihood (CVPOL)", "\n")
	cat("    (CVPOL is the bias-corrected MPOL obtained by approximated cross-validation)", "\n")
	}
	cat(" \n")
	print(x$EPOCE,na.print=".",row.names=F)
	cat(" \n")
}
