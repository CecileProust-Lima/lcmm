\name{plot.cuminc}
\alias{plot.cuminc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot of predicted cumulative incidences according to a profile of covariates
}
\description{
This function displays the predicted cause-specific cumulative incidences derived from a joint latent class model according to a profile of covariates. 
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
\method{plot}{cuminc}(x, profil, event, add = FALSE, legend, legend.loc = "topleft",  ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
an object of class \code{cuminc}
}
  \item{profil}{
an integer giving the profile number for which the cumulative incidences are to be plotted.
}
  \item{event}{
an integer giving the event indicator for which the cumulative incidence are to be plotted.
}
  \item{add}{
logical indicating if the curves should be added to an existing plot. Default to FALSE.
}
  \item{legend}{
character or expression to appear in the legend. If no legend should be added, \code{"legend"} should be NULL.
}
  \item{legend.loc}{
keyword for the position of the legend from the list \code{"bottomright"}, \code{"bottom"}, \code{"bottomleft"},
     \code{"left"}, \code{"topleft"},\code{"top"}, \code{"topright"}, \code{"right"} and
     \code{"center"}. By default, the legend is located in the top left of the plot.
}
  \item{\dots}{
    other parameters to be passed through to plotting functions
    }
}

\value{
returns NULL
}

\author{
Viviane Philipps and Cecile Proust-Lima
}

\seealso{
\code{\link{Jointlcmm}},\code{\link{plot.Jointlcmm}},\code{\link{cuminc}}
}
