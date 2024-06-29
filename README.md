# SurvBal

## Overview
`SurvBal` enables the selection of microbiome balances in relation to censored survival or time-to-event outcomes, which are of considerable interest in many biomedical studies. The most commonly used survival models – the Cox proportional hazards and parametric survival (including accelerated failure time) models are included in the package, which are used in combination with step-wise selection procedures to identify the optimal associated balance of microbiome, i.e., the ratio of the geometric means of two groups of taxa’s relative abundances.

## Installation Guide

Please make sure you install the R package `zCompositions` with version number <= 1.4.0.1 before installing `SurvBal`: 

```
require(remotes) 
install_version("zCompositions", version = "1.4.0.1", repos = "http://cran.us.r-project.org")
```

Then you can install `SurvBal` from GitHub by:

```
# install.packages("devtools")
devtools::install_github("yinglia/SurvBal")
```

Please skip the recommended update about `zCompositions`.


## Resources


All functions in `SurvBal` are described in the manual: 

https://github.com/yinglia/SurvBal/blob/main/SurvBal_1.1.0.pdf

The vignette of `SurvBal` includes the instructions to run the method on the example dataset, with the input and output clearly described, which can be found at: 

https://yinglia.github.io/SurvBal-Vignette/

A Shiny app of `SurvBal`: 

https://yinglistats.shinyapps.io/shinyapp-survbal/
