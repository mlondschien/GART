
<!-- README.md is generated from README.Rmd. Please edit that file -->

# GART

<!-- badges: start -->
<!-- badges: end -->

Guided Adversarial Robust Transfer (GART) Learning aims to effectively
bridges the realms of transfer learning and distributional robustness
prediction models when dealing with a limited amount of target data and
a diverse range of source models. By leveraging the source mixing
assumption, GART is designed to learn valuable knowledge that may be
present in different yet potentially related auxiliary samples, and
achieve a faster convergence rate than the model fitted with the target
data.

## Installation

You can install the development version of GART from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("xinxiong0238/GART")
```

## Example

The following is a basic example to run GART. We first generate the
training data (a small number of target data and 4 large but
heterogenerous source data) as well as a validation data set. For
illustration purpose, the validation data share the same generation
mechanism as the training target data. Then `GART` function is called to
estimate GART parameter and validate model performance. Five benchmarks
are also included (i.e., target only estimator, source mixture, maximin,
transLasso and transGLM).

``` r
library(GART)
data_sim = simu_data()
data = data_sim$GART_est_input
data_valid = data_sim$GART_eval_input

fit = GART(data, is_benchmark = T, is_valid = T, data_valid = data_valid)
#> Estimate source coef...
#> Estimate GART coef...
#> Add benchmark methods (maximin, transLasso, transGLM, etc)...
#> Valid model performance..
```
