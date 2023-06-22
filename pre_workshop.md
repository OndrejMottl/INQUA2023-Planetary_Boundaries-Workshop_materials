Pre-workshop
================

Pre-workshop instructions for participants.

## R and RStudio

[R](https://en.wikipedia.org/wiki/R_(programming_language)) is software for statistical analyses. [R-studio](https://posit.co/products/open-source/rstudio/) is software to help a user use the `R` (i.e. integrated development environment; IDE).

There are many guides on how to obtain and/or update R and RStudio, for example:

-   [R Basics for Paleoecologists](https://ckiahtipes.github.io/) by C.A. Kiahtipes, a previous part of the APD series of workshops.
-   [Install or Update R tutorial](https://jennhuck.github.io/workshops/install_update_R.html) by Jennifer Huck

## Packages

Packages are tools that contain a series of functions for a specific need. We need to make sure that all packages used throughout the workshop are installed on everyone’s computer.

### Install packages

Make a list of packages needed from CRAN

``` r
package_list <-
  c(
    "tidyverse", # general data wrangling and visualisation
    "pander", # nice tables
    "Bchron", # age-depth modelling
    "janitor", # string cleaning
    "remotes", # installing packages from GitHub
    "neotoma2", # access to the Neotoma database
    "mgcv", # GAM fitting
    "gratia" # GAM visualisation
  )
```

Install all packages from CRAN

``` r
lapply(
  package_list, utils::install.packages
)
```

Install packages from GitHub

``` r
# Install R-Ecopol
remotes::install_github("HOPE-UIB-BIO/R-Ecopol-package")
# Note that the package can be only installed on Windows machine. 
```

### Test if everything is set-up

The following code should produce `"Everything is good to go"`, not error (`"All required packages are not installed"`).

``` r
if (
  isTRUE(
    all(
      c(package_list, "REcopol") %in%
        as.data.frame(
          utils::installed.packages()
        )[, 1]
    )
  )
) {
  cat("Everything is good to go")
} else {
  warning("All required packages are not installed")
}
```
