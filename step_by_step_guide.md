Step-by-step guide
================

This workflow should serve as step-by-step guidance starting from downloading a dataset from Neotoma and processing it, to estimating ecosystem property (diversity in this case).

:warning: **This workflow is only meant as an example**: There are several additional steps for data preparation which should be done for any fossil pollen dataset from Neotoma!

See [**FOSSILPOL**](https://hope-uib-bio.github.io/FOSSILPOL-website/), an R-based modular workflow to process multiple fossil pollen records to create a comprehensive, standardised dataset compilation, ready for multi-record and multi-proxy analyses at various spatial and temporal scales.

## Install packages

Please follow the [pre-workshop instructions](/docs/pre_workshop.html) to make sure all packages are installed.

## Attach packages

``` r
library(tidyverse) # general data wrangling and visualisation ✨
library(pander) # nice tables 😍
library(neotoma2) # access to the Neotoma database 🌿
library(Bchron) # age-depth modelling 🕰️
library(mgcv) # GAM fitting 📈
library(marginaleffects) # predicting trends 📈
library(janitor) # string cleaning 🧹
library(here) # for working directory 🗺️
```

``` r
ggplot2::theme_set(
  ggplot2::theme_bw() +
    ggplot2::theme(
      axis.title = ggplot2::element_text(size = 25),
      axis.text = ggplot2::element_text(size = 15),
      strip.text = ggplot2::element_text(size = 15),
      panel.grid = ggplot2::element_blank()
    )
)
```

## Download a dataset from Neotoma

Here we have selected the **Lingua d’Oca** (ID = 52162) record by Feredico Di Rita.

Reference paper: Di Rita, F., A. Celant, and D. Magri. 2010. Holocene environmental instability in the wetland north of the Tiber delta (Rome, Italy): sea-lake-man interactions. *Journal of Paleolimnology* 44:51-67.

``` r
sel_dataset_download <-
  neotoma2::get_downloads(52162)
```

## Prepare the pollen counts

``` r
# get samples
sel_counts <-
  neotoma2::samples(sel_dataset_download)

# select only "pollen" taxa
sel_taxon_list_selected <-
  neotoma2::taxa(sel_dataset_download) %>%
  dplyr::filter(element == "pollen") %>%
  purrr::pluck("variablename")

# prepare taxa table
sel_counts_selected <-
  sel_counts %>%
  as.data.frame() %>%
  dplyr::mutate(sample_id = as.character(sampleid)) %>%
  tibble::as_tibble() %>%
  dplyr::select("sample_id", "value", "variablename") %>%
  # only include selected taxons
  dplyr::filter(
    variablename %in% sel_taxon_list_selected
  ) %>%
  # turn into the wider format
  tidyr::pivot_wider(
    names_from = "variablename",
    values_from = "value",
    values_fill = 0
  ) %>%
  # clean names
  janitor::clean_names()

head(sel_counts_selected)[, 1:5]
```

| sample_id | fagus | picea | rumex | salix |
|:---------:|:-----:|:-----:|:-----:|:-----:|
|  510392   |   1   |   1   |   1   |   1   |
|  510393   |   3   |   0   |   0   |   0   |
|  510394   |   6   |   0   |   1   |   0   |
|  510396   |   5   |   0   |   0   |   0   |
|  510391   |   8   |   0   |   0   |   0   |
|  510395   |   6   |   0   |   0   |   0   |

Here, we strongly advocate that attention should be paid to the selection of the ecological groups, the selection of depositional environments, as well as the harmonisation of the pollen taxa. However, that is not the subject of this workflow, but any analysis to be published needs careful preparation of the fossil pollen datasets!

We can now try to visualise the taxa per sample_id

``` r
sel_counts_selected %>%
  tibble::rowid_to_column("ID") %>%
  tidyr::pivot_longer(
    cols = -c(sample_id, ID),
    names_to = "taxa",
    values_to = "n_grains"
  ) %>%
  ggplot2::ggplot(
    mapping = ggplot2::aes(
      x = ID,
      y = n_grains,
      fill = taxa
    ),
  ) +
  ggplot2::geom_bar(
    stat = "identity",
    position = "fill"
  ) +
  ggplot2::labs(
    x = "sample_id",
    y = "proportion of pollen grains"
  ) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_blank(),
    legend.position = "none"
  )
```

![](step_by_step_guide_files/figure-gfm/count_vis-1.png)

## Preparation of the levels

### Sample depth

Extract depth for each level

``` r
sel_level <-
  neotoma2::samples(sel_dataset_download) %>%
  tibble::as_tibble() %>%
  dplyr::mutate(sample_id = as.character(sampleid)) %>%
  dplyr::distinct(sample_id, depth) %>%
  dplyr::relocate(sample_id)

head(sel_level)
```

| sample_id | depth |
|:---------:|:-----:|
|  510392   |   5   |
|  510393   |  10   |
|  510394   |  15   |
|  510396   |  20   |
|  510391   |  25   |
|  510395   |  30   |

### Age-depth modelling

We highly recommend recalculating the age-depth model ‘de Novo’ as different methodologies might have been applied to each record. However, age-depth modelling is a very complicated topic, which we will not dig into today. Just note that many parts of the age-depth modelling need attention (selection of chronology control points, using correct calibration curves, etc.).

``` r
# Here we only present a few of the important steps of preparation of the
#   chronology control table. There are many more potential issues, but
#   solving those is not the focus of this workflow.

# First, get the chronologies and check which we want to use used
sel_chron_control_table_download <-
  neotoma2::chroncontrols(sel_dataset_download)

# prepare the table
sel_chron_control_table <-
  sel_chron_control_table_download %>%
  # Here select the ID of one of the chronology
  dplyr::filter(chronologyid == 37228) %>%
  tibble::as_tibble() %>%
  # Here we calculate the error as the average of the age `limitolder` and
  #   `agelimityounger`
  dplyr::mutate(
    error = round((agelimitolder - agelimityounger) / 2)
  ) %>%
  # As Bchron cannot accept an error of 0, we need to replace the value with 1
  dplyr::mutate(
    error = replace(error, error == 0, 1),
    error = ifelse(is.na(error), 1, error)
  ) %>%
  # As Bchron cannot accept an thickness of 0, we need to replace the value with 1
  dplyr::mutate(
    thickness = ifelse(is.na(thickness), 1, thickness)
  ) %>%
  # We need to specify which calibration curve should be used for what point
  dplyr::mutate(
    curve = ifelse(as.data.frame(sel_dataset_download)["lat"] > 0, "intcal20", "shcal20"),
    curve = ifelse(chroncontroltype != "Radiocarbon", "normal", curve)
  ) %>%
  tibble::column_to_rownames("chroncontrolid") %>%
  dplyr::arrange(depth) %>%
  dplyr::select(
    chroncontrolage, error, depth, thickness, chroncontroltype, curve
  )

i_multiplier <- 0.1 # increase to 5

# Those are default values suggested by the Bchron package
n_iteration_default <- 10e3
n_burn_default <- 2e3
n_thin_default <- 8

# Let's multiply them by our i_multiplier
n_iteration <- n_iteration_default * i_multiplier
n_burn <- n_burn_default * i_multiplier
n_thin <- max(c(1, n_thin_default * i_multiplier))

# run Bchron
ad_model <-
  Bchron::Bchronology(
    ages = sel_chron_control_table$chroncontrolage,
    ageSds = sel_chron_control_table$error,
    positions = sel_chron_control_table$depth,
    calCurves = sel_chron_control_table$curve,
    positionThicknesses = sel_chron_control_table$thickness,
    iterations = n_iteration,
    burn = n_burn,
    thin = n_thin
  )
```

Visually check the age-depth models

``` r
plot(ad_model)
```

![](step_by_step_guide_files/figure-gfm/bchron_figure-1.png)

#### Predict ages

``` r
age_position <-
  Bchron:::predict.BchronologyRun(object = ad_model, newPositions = sel_level$depth)

age_uncertainties <-
  age_position %>%
  as.data.frame() %>%
  dplyr::mutate_all(., as.integer) %>%
  as.matrix()

colnames(age_uncertainties) <- sel_level$sample_id

# Let's take the median age of all possible ages (i.e. the estimated age
#   from each age-depth model run) as our default.
sel_level_predicted <-
  sel_level %>%
  dplyr::mutate(
    age = apply(
      age_uncertainties, 2,
      stats::quantile,
      probs = 0.5
    )
  )
```

### Visualisation of our data

Let’s now make a simple pollen diagram with proportions of the main pollen taxa (x-axis) against our estimated ages along depth (y-axis).

``` r
data_rownames <-
  sel_counts_selected %>%
  tibble::column_to_rownames("sample_id")

data_percentages <-
  (data_rownames / rowSums(data_rownames)) * 100

col_sum_non_zero <-
  colSums(data_percentages) > 0

data_filtered <-
  data_percentages %>%
  dplyr::select(
    dplyr::any_of(
      names(col_sum_non_zero)[col_sum_non_zero]
    )
  ) %>%
  tibble::rownames_to_column("sample_id") %>%
  dplyr::relocate(sample_id) %>%
  tibble::tibble() 

data_filtered %>%
  dplyr::inner_join(
    sel_level_predicted,
    by = dplyr::join_by(sample_id)
  ) %>%
  tidyr::pivot_longer(
    cols = -c(sample_id, depth, age),
    names_to = "taxa",
    values_to = "proportion_of_grains"
  ) %>%
  dplyr::group_by(taxa) %>%
  # Calculate the average proportion of grains
  dplyr::mutate(
    avg_prop = mean(proportion_of_grains)
  ) %>%
  # Only keep the main taxa (on average >= 1% pollen grains)
  dplyr::filter(avg_prop >= 1) %>%
  dplyr::ungroup() %>%
  ggplot2::ggplot(
    mapping = ggplot2::aes(
      y = age,
      x = proportion_of_grains,
      xmax = proportion_of_grains,
      xmin = 0,
      fill = taxa,
      col = taxa
    ),
  ) +
  ggplot2::geom_ribbon() +
  ggplot2::scale_y_continuous(trans = "reverse") +
  ggplot2::scale_x_continuous(breaks = c(0, 1)) +
  ggplot2::facet_wrap(~taxa, nrow = 1) +
  ggplot2::theme(
    legend.position = "none",
    panel.border = ggplot2::element_blank(),
    strip.background = ggplot2::element_blank(),
    strip.text = ggplot2::element_blank(),
    axis.ticks.x = ggplot2::element_blank(),
    axis.text.x = ggplot2::element_blank()
  ) +
  ggplot2::labs(
    x = "Proportion of grains",
    y = "Age (cal yr BP)"
  )
```

![](step_by_step_guide_files/figure-gfm/vis_data_with_ages-1.png)

## Estimation of ecosystem property

Now we will use our prepared fossil pollen data to estimate the diversity. We will use {REcopol} package, which has easy-to-use functions to analyse fossil pollen data. See package [website](https://hope-uib-bio.github.io/R-Ecopol-package/) for more information. Specifically, we will estimate rarefied values of [Hill numbers](https://esajournals.onlinelibrary.wiley.com/doi/abs/10.2307/1934352).

``` r
# function to estimate diversity
get_diversity <- function(data_source, round = TRUE, sel_method = "") {
  # helper function
  get_diversity_taxonomic <- function(data_matrix, sample_size) {
    hill0 <- function(data, sample_size) {
      data_sub <- data[data > 0]
      data_sum <- sum(data_sub)
      if (sample_size <= data_sum) {
        res <- sum(1 - exp(lchoose(data_sum - data_sub, sample_size) -
          lchoose(data_sum, sample_size)))
        return(res)
      } else {
        return(0)
      }
    }
    fk_hat <- function(data, sample_size) {
      data_sub <- data[data > 0]
      data_sum <- sum(data_sub)
      if (sample_size <= data_sum) {
        sub <- function(k) {
          sum(exp(lchoose(data_sub, k) + lchoose(data_sum -
            data_sub, sample_size - k) - lchoose(
            data_sum,
            sample_size
          )))
        }
        res <- sapply(1:sample_size, sub)
        return(res)
      } else {
        return(0)
      }
    }
    hill1 <- function(data, sample_size) {
      data_sub <- data[data > 0]
      data_sum <- sum(data_sub)
      if (sample_size <= data_sum) {
        k <- 1:sample_size
        res <- exp(-sum(k / sample_size * log(k / sample_size) *
          fk_hat(data_sub, sample_size)))
        return(res)
      } else {
        return(0)
      }
    }
    hill2 <- function(data, sample_size) {
      data_sub <- data[data > 0]
      data_sum <- sum(data_sub)
      if (sample_size <= data_sum) {
        res <- 1 / (1 / sample_size + (1 - 1 / sample_size) * sum(data_sub *
          (data_sub - 1) / data_sum / (data_sum - 1)))
        return(res)
      } else {
        return(0)
      }
    }
    est_n_0 <-
      sapply(sample_size, function(n) {
        apply(data_matrix, 1, hill0, sample_size = n)
      })
    est_n_1 <-
      sapply(sample_size, function(n) {
        apply(data_matrix, 1, hill1, sample_size = n)
      })
    est_n_2 <-
      sapply(sample_size, function(n) {
        apply(data_matrix, 1, hill2, sample_size = n)
      })
    hill_diversity <-
      cbind(est_n_0, est_n_1, est_n_2, est_n_1 -
        est_n_2, est_n_2 / est_n_1, est_n_1 / est_n_0) %>%
      as.data.frame() %>%
      tibble::as_tibble() %>%
      dplyr::mutate(sample_id = row.names(data_matrix)) %>%
      tibble::column_to_rownames("sample_id") %>%
      purrr::set_names(
        nm = c(
          "n0",
          "n1", "n2", "n1_minus_n2", "n2_divided_by_n1", "n1_divided_by_n0"
        )
      )
    return(hill_diversity)
  }

  if (
    isTRUE(round)
  ) {
    data_matrix <-
      data_source %>%
      tibble::column_to_rownames("sample_id") %>%
      dplyr::mutate_all(., .f = floor) %>%
      as.matrix() %>%
      round()
  } else {
    data_matrix <-
      data_source %>%
      tibble::column_to_rownames("sample_id") %>%
      as.matrix()
  }

  sample_size <-
    apply(data_matrix, 1, sum) %>%
    floor() %>%
    min()

  div <-
    get_diversity_taxonomic(
      data_matrix = data_matrix,
      sample_size = sample_size
    )

  res <-
    div %>%
    tibble::rownames_to_column("sample_id") %>%
    dplyr::relocate(sample_id) %>%
    as.data.frame() %>%
    tibble::as_tibble()
  return(res)
}
```

``` r
data_diversity <-
  get_diversity(
    data_source = sel_counts_selected,
    sel_method = "taxonomic"
  )

head(data_diversity)
```

| sample_id |  n0   |  n1   |  n2   | n1_minus_n2 | n2_divided_by_n1 |
|:---------:|:-----:|:-----:|:-----:|:-----------:|:----------------:|
|  510392   | 25.17 | 12.76 | 8.095 |    4.664    |      0.6344      |
|  510393   | 24.39 | 11.66 | 6.795 |    4.867    |      0.5826      |
|  510394   | 24.33 | 11.13 | 6.707 |    4.426    |      0.6024      |
|  510396   | 20.29 | 7.68  | 3.952 |    3.728    |      0.5146      |
|  510391   | 26.27 | 13.12 | 7.43  |    5.694    |      0.5661      |
|  510395   | 22.78 | 11.89 | 7.615 |    4.279    |      0.6402      |

Table continues below

| n1_divided_by_n0 |
|:----------------:|
|      0.5068      |
|      0.4782      |
|      0.4576      |
|      0.3784      |
|      0.4995      |
|      0.522       |

Now we can fit a temporal trend as a GAM model.

``` r
data_to_fit <-
  dplyr::inner_join(
    data_diversity,
    sel_level_predicted,
    by = "sample_id"
  )

mod_n0 <-
  mgcv::gam(
    n0 ~ s(age, k = 25, bs = "tp"),
    data = data_to_fit,
    method = "REML",
    family = mgcv::tw(link = "log")
  )

summary(mod_n0)
#> 
#> Family: Tweedie(p=1.01) 
#> Link function: log 
#> 
#> Formula:
#> n0 ~ s(age, k = 25, bs = "tp")
#> 
#> Parametric coefficients:
#>             Estimate Std. Error t value Pr(>|t|)    
#> (Intercept)   2.9714     0.0221   134.5   <2e-16 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Approximate significance of smooth terms:
#>          edf Ref.df     F p-value    
#> s(age) 4.667  5.828 14.31  <2e-16 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> R-sq.(adj) =  0.596   Deviance explained = 60.1%
#> -REML = 178.17  Scale est. = 0.58249   n = 64

# mgcv::gam.check(mod_n0, k.sample = 10e3, k.rep = 1e3)
```

Now we can visualise the results.

``` r
age_dummy <-
  tibble::tibble(
    age = seq(
      from = min(data_to_fit$age),
      to = max(data_to_fit$age),
      length.out = 100
    )
  )

data_predicted <-
  marginaleffects::predictions(
    model = mod_n0,
    newdata = age_dummy
  ) %>%
  tibble::as_tibble() %>%
  dplyr::rename(
    fit = estimate,
    lwr = conf.low,
    upr = conf.high,
    sd_error = std.error
  ) %>%
  dplyr::select(
    dplyr::all_of(
      c(
        "age",
        "fit",
        "sd_error",
        "lwr",
        "upr"
      )
    )
  )

data_predicted %>%
  ggplot2::ggplot(
    mapping = ggplot2::aes(
      x = age,
      y = fit
    )
  ) +
  ggplot2::geom_ribbon(
    mapping = ggplot2::aes(
      ymin = lwr,
      ymax = upr
    ),
    alpha = 0.2
  ) +
  ggplot2::geom_line() +
  geom_point(
    data = data_to_fit,
    mapping = ggplot2::aes(
      y = n0
    )
  ) +
  ggplot2::scale_x_continuous(
    trans = "reverse"
  ) +
  ggplot2::labs(
    y = "Hill's N0",
    x = "Age (cal yr BP)"
  )
```

![](step_by_step_guide_files/figure-gfm/plot_gam-1.png)
