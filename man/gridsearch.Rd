\name{gridsearch}
\alias{gridsearch}
\title{Automatic grid search}
\description{
  This function provides an automatic grid search for latent class mixed
  models estimated with \code{hlme}, \code{lcmm}, \code{multlcmm} and
  \code{Jointlcmm} functions.
}
\usage{
gridsearch(m, rep, maxiter, minit)
}
\arguments{
  \item{m}{a call of \code{hlme}, \code{lcmm}, \code{multlcmm} or
    \code{Jointlcmm} corresponding to the model to estimate
}
  \item{rep}{the number of departures from random initial values
}
  \item{maxiter}{the number of iterations in the optimization algorithm
}

  \item{minit}{an object of class  \code{hlme}, \code{lcmm},
    \code{multlcmm} or \code{Jointlcmm} corresponding to the same model
    as specified in m except for the number of classes (it should be
    one). This object is used to generate random initial values}
}
\value{an object of class \code{hlme}, \code{lcmm}, \code{multlcmm} or
    \code{Jointlcmm} corresponding to the call specified in m.
}
\author{Cecile Proust-Lima and Viviane Philipps
}

\details{
The function permits the estimation of a model from a grid of random
initial values to reduce the odds of a convergence towards a local maximum.

The function was inspired by the emEM technique described in Biernacki
et al. (2003). It consists in:

1. randomly generating \code{rep} sets of initial values for \code{m} from
the estimates of \code{minit} (this is done internally using option
\code{B=random(minit)} \code{rep} times)

2. running the optimization algorithm for the model specified in \code{m}
from the \code{rep} sets of initial values with a maximum number of
iterations of \code{maxit} each time.

3. retaining the estimates of the random initialization that provides the
best log-likelihood after \code{maxiter} iterations.

4. running the optimization algorithm from these estimates for the final estimation.  
}
\references{
Proust-Lima C, Philipps V, Liquet B (2017). Estimation of Extended Mixed
Models Using Latent Classes and Latent Processes: the R package lcmm.
Journal of Statistical Software, 78(2), 1-56

Biernacki C, Celeux G, Govaert G (2003). Choosing Starting Values for
the EM Algorithm for Getting the Highest Likelihood in Multivariate
Gaussian Mixture models. Computational Statistics and Data Analysis, 41(3-4), 561-575.
}
\examples{
\dontrun{
# initial model with ng=1 for the random initial values
m1 <- hlme(Y ~ Time * X1, random =~ Time, subject = 'ID', ng = 1, 
      data = data_hlme)

# gridsearch with 10 iterations from 50 random departures
m2d <- gridsearch(rep = 50, maxiter = 10, minit = m1, hlme(Y ~ Time * X1,
      mixture =~ Time, random =~ Time, classmb =~ X2 + X3, subject = 'ID',
          ng = 2, data = data_hlme))
        }
}

