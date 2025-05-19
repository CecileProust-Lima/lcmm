# <img src="vignettes/lcmm.png" align="right" width=130 style="margin-right: 0px;vertical-align:middle"/> <span style="font-size:38px"> Extended Mixed Models Using Latent Classes and Latent Processes </span>

&nbsp;

<p align="justify">
The lcmm package implements various extensions of mixed models. It handles continuous (Gaussian or not) and ordinal outcomes with repeated measures and latent classes. A time-to-event can jointly be considered in a proportionnal hazard model.
All models are estimated with a maximum likelihood framework using a modified Marquardt-Levenberg algorithm.
The package also includes several predictions, visualization, and utility functions to conduct a complete statistical analysis.
</p>

A detailed companion paper is available in Journal of Statistical Software :

<p align="justify">
Proust-Lima C, Philipps V, Liquet B. Estimation of Extended Mixed Models
Using Latent Classes and Latent Processes: The R Package lcmm. Journal
of Statistical Software, Articles. 2017;78(2):1-56.
<https://doi.org/10.18637/jss.v078.i02>
</p>

<p align="justify">
And specific statistical models estimated are described in various statistical papers of the authors.
</p>

&nbsp;


## Install the package

The lcmm package needs version 3.5 or newer of the R software.
To install the released CRAN version of the package, use

``` r
install.packages("lcmm")
```

To get the most recent update, install it from github :

``` r
remotes::install_github("CecileProust-Lima/lcmm")
```

The lcmm package depends on other R package, namely :

- survival (>=2.37-2) for dealing with the survival outcomes
- parallel and doParallel for parallelizing some time consuming functions
- mvtnorm for generating random parameters
- spacefillr for the quasi Monte Carlo sequences
- marqLevAlg (>2.0) for the numerical optimization
- numDeriv for computing the Hessian


To run the examples proposed in this website, the following package are also needed :

- lattice
- NormPsy
- ggplot2
- ggpubr
- dplyr
- splines
- gridExtra
- ggalluvial

&nbsp;

## Documentation

<p align="justify">
This website is intended to help the lcmm users in their statistical analyses. It provides an [overview](articles/lcmm.html) of the package, several [vignettes](articles/index.html), a [FAQ](articles/usual_problems.html) page and the [help pages](reference/index.html) of all functions included in the lcmm package.
</p>

<p align="justify">
Further issues and questions about the use of the lcmm package are reported on the github issue page <https://github.com/CecileProust-Lima/lcmm/issues>.
Please check both opened and closed issues to make sure that the topic has not already been treated before creating a new issue. To report a bug, please provide a reproducible example.
</p>
