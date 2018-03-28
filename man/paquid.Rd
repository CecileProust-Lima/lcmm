\name{paquid}
\alias{paquid}
\docType{data}
\title{
Longitudinal data on cognitive and physical aging in the elderly
}
\description{
The dataset consists in a subsample of the Paquid prospective cohort study. Repeated measures cognitive measures (MMSE, IST, BVRT psychometric tests), physical dependency (HIER) and depression sympatomatology (CESD) were collected over a maximum period of 20 years along with dementia information (age at dementia diagnosis, dementia diagnosis information). Time-independent socio-demographic information is also provided (CEP, male, age_init).
}
\usage{paquid}
\format{
  A data frame with 2250 observations over 500 subjects and 12 variables:
  \describe{
    \item{\code{ID}}{subject identification number}
    \item{\code{MMSE}}{score at the Mini-Mental State Examination (MMSE), a psychometric test of global cognitive functioning (integer in range 0-30)}
    \item{\code{BVRT}}{score at the Benton Visual Retention Test (BVRT), a psychometric test of spatial memory (integer in range 0-15)}
    \item{\code{IST}}{score at the Isaacs Set Test (IST) truncated at 15 seconds, a test of verbal memory (integer in range 0-40)}
    \item{\code{HIER}}{score of physical dependency (0=no dependency, 1=mild dependency, 2=moderate dependency, 3=severe dependency)}
    \item{\code{CESD}}{score of a short self-report scale CES-D designed to measure depressive symptomatology in the general population (integer in range 0-52)}
    \item{\code{age}}{age at the follow-up visit}
    \item{\code{dem}}{indicator of positive diagnosis of dementia}
    \item{\code{agedem}}{age at dementia diagnosis for \code{dem=1} and at last contact for \code{dem=0}}
    \item{\code{age_init}}{age at entry in the cohort}
    \item{\code{CEP}}{binary indicator of educational level (\code{CEP=1} for subjects who graduated from primary school; \code{CEP=0} otherwise)}
    \item{\code{male}}{binary indicator for gender (\code{male=1} for men; \code{male=0} for women)}
  }
}
\references{
Letenneur, L., Commenges, D., Dartigues, J. F., & Barberger-Gateau, P. (1994). Incidence
of dementia and Alzheimer's disease in elderly community residents of southwestern
France. International Journal of Epidemiology, 23 (6), 1256-61.
}
\examples{
summary(paquid)
}
\keyword{datasets}
