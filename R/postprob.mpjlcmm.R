#' @export
postprob.mpjlcmm <- function(x,threshold=c(0.7,0.8,0.9),...)
{
    if (!inherits(x, "mpjlcmm")) stop("use only with \"mpjlcmm\" objects")

    xx <- x
    class(xx) <- "Jointlcmm"

    postprob.Jointlcmm(xx,threshold=threshold,...)
 }
