\name{summarytable}
\alias{summarytable}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Summary of models
}
\description{
This function provides a table summarizing the results of different
models fitted by \code{hlme}, \code{lcmm}, \code{multlcmm} or \code{Jointlcmm}.
}
\usage{
summarytable(m1,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{m1}{
an object of class \code{hlme}, \code{lcmm}, \code{multlcmm} or \code{Jointlcmm}
}
  \item{...}{
    further arguments, in particular other objects of class \code{hlme}, \code{lcmm}, \code{multlcmm} or \code{Jointlcmm}
}
}

\value{
a matrix giving for each model: the number a latent classes, the log-likelihood, the
number of parameters, the BIC and the posterior probability of the latent
classes.
}

\author{
Cecile Proust-Lima, Viviane Philipps
}
\seealso{
\code{\link{summary}}, \code{\link{hlme}}, \code{\link{lcmm}}, \code{\link{multlcmm}}, \code{\link{Jointlcmm}}
}
