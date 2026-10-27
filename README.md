# Re-analysis of Gade et al. (2021)
Shaheed Azaad

- [Overview](#overview)
- [Reproducing Gade et al.’s original
  results](#reproducing-gade-et-als-original-results)
  - [Model with interaction](#model-with-interaction)
  - [Model without interaction](#model-without-interaction)
    - [Model comparison](#model-comparison)
    - [Fixed effects](#fixed-effects)
  - [Summary](#summary)
- [After correccting](#after-correccting)
  - [Model with interaction](#model-with-interaction-1)
  - [Model without interaction](#model-without-interaction-1)
    - [Model comparison](#model-comparison-1)
    - [Fixed effects](#fixed-effects-1)
    - [Summary](#summary-1)

# Overview

This document contains (1) a minimal reproduction of Gade et al.’s
(2021; Journal of cognition) meta-analysis, and (2) a re-analysis based
after fixing a critical bug in their control flow.

The main issue is this: the code contains an `if` statement that is
supposed to fit models with condition random effects if each study
contributes more than one effect. However, there is a bug in the code
that means only models without condition random effects were fit,
contrary to what is reported in the paper.

To re-run this QMD, place it inside Gade et al.’s OSF folder, found
here: https://osf.io/ukjq4/

``` r
library(tidyr)
library("brms")
```

    Loading required package: Rcpp

    Loading 'brms' package (version 2.23.0). Useful instructions
    can be found by typing help('brms'). A more detailed introduction
    to the package is available through vignette('brms_overview').


    Attaching package: 'brms'

    The following object is masked from 'package:stats':

        ar

``` r
library("ggplot2")
library("bayesplot")
```

    This is bayesplot version 1.16.0

    - Online documentation and vignettes at mc-stan.org/bayesplot

    - bayesplot theme set to bayesplot::theme_default()

       * Does _not_ affect other ggplot2 plots

       * See ?bayesplot_theme_set for details on theme setting


    Attaching package: 'bayesplot'

    The following object is masked from 'package:brms':

        rhat

``` r
temp <- read.table("Data R2/temp.txt", header=T, sep="\t") # the datafile created by the authors' processing code
options(mc.cores = 6)
set.seed(123) # the original code has no seed so numbers wont perfectly reproduce
```

# Reproducing Gade et al.’s original results

The below code fits models from Gade et al.’s script `analysis_code.R`.
I only run the main analyses here (not the moderator tests).

I have made *one* modification to the code. I have added print
statements to the relevant code blocks. So, if the model without random
effects runs, we will see “WITHOUT CONDITION RANDOM EFFECTS” printed in
the output, and if it does, we will see “WITH CONDITION RANDOM EFFECTS”
printed.

## Model with interaction

``` r
prior <- get_prior(
  dm ~ 1 +
    dominance +
    language_transition +
    dominance * language_transition +
    (Condition | Study),
  temp,
  family = gaussian(),
  autocor = NULL,
  sparse = NULL,
  internal = TRUE
)

if (length(temp$dm == 4 * length(unique(temp$Study)))) {
  print("WITHOUT CONDITION RANDOM EFFECTS")  
  test <- brm(
    formula = dm ~ 1 +
      dominance +
      language_transition +
      dominance * language_transition +
      (1 | Study),
    data = temp,
    family = gaussian(),
    inits = 0,

    prior = c(
      set_prior("student_t(3, 847.5, 215)", class = "Intercept"),
      set_prior("normal(0,10)", class = "b"),
      set_prior("student_t(3, 0, 215)", class = "sd")
    ),
    warmup = 1000,
    iter = 5000,
    chains = 4,
    silent = 2,
    refresh = 0,
    control = list(adapt_delta = 0.95, max_treedepth = 15), 
    open_progress = FALSE
  )
} else if (length(temp$dm != 4 * length(unique(temp$Study)))) {
  print("WITH CONDITION RANDOM EFFECTS")  
  test <- brm(
    formula = dm ~ 1 +
      dominance +
      language_transition +
      dominance * language_transition +
      (Condition | Study),
    data = temp,
    family = gaussian(),
    inits = 0,

    prior = c(
      set_prior("student_t(3, 847.5, 215)", class = "Intercept"),
      set_prior("normal(0,10)", class = "b"),
      set_prior("student_t(3, 0, 215)", class = "sd")
    ),
    warmup = 1000,
    iter = 5000,
    chains = 4,
    silent = 2,
    refresh = 0,
    cores = 4,
    control = list(adapt_delta = 0.95, max_treedepth = 15)
  )
}
```

    [1] "WITHOUT CONDITION RANDOM EFFECTS"

    Warning: Argument 'inits' is deprecated. Please use argument 'init' instead.

``` r
modelwinteraction <- test

# fit criteria

modelwinteraction <- add_criterion(
  modelwinteraction,
  c("loo", "waic", "bayes_R2")
)
```

    Warning: 
    39 (4.1%) p_waic estimates greater than 0.4. We recommend trying loo instead.

``` r
modelwinteraction$criteria$waic
```


    Computed from 16000 by 956 log-likelihood matrix.

              Estimate   SE
    elpd_waic  -5729.5 31.0
    p_waic        65.0  4.7
    waic       11459.0 61.9

    39 (4.1%) p_waic estimates greater than 0.4. We recommend trying loo instead. 

``` r
modelwinteraction$criteria$loo
```


    Computed from 16000 by 956 log-likelihood matrix.

             Estimate   SE
    elpd_loo  -5730.2 31.0
    p_loo        65.8  4.8
    looic     11460.5 62.0
    ------
    MCSE of elpd_loo is 0.1.
    MCSE and ESS estimates assume MCMC draws (r_eff in [0.4, 2.0]).

    All Pareto k estimates are good (k < 0.7).
    See help('pareto-k-diagnostic') for details.

``` r
mean(modelwinteraction$criteria$bayes_R2)
```

    [1] 0.8105422

``` r
#model fit
pp_check(modelwinteraction, type = "scatter_avg") +
  geom_abline(intercept = 0, slope = 1, color = "red", lty = 4) # y vs. yrep
```

    Using all posterior draws for ppc type 'scatter_avg' by default.

![](gade_files/figure-commonmark/unnamed-chunk-2-1.png)

## Model without interaction

``` r
if (length(temp$dm == 4 * length(unique(temp$Study)))) {
  print("WITHOUT CONDITION RANDOM EFFECTS")      
  test <- brm(
    formula = dm ~ 1 + dominance + language_transition + (1 | Study),
    data = temp,
    family = gaussian(),
    inits = 0,

    prior = c(
      set_prior("student_t(3, 847.5, 215)", class = "Intercept"),
      set_prior("normal(0,10)", class = "b"),
      set_prior("student_t(3, 0, 215)", class = "sd")
    ),
    warmup = 1000,
    iter = 5000,
    chains = 4,
    silent = 2,
    refresh = 0,
    control = list(adapt_delta = 0.95, max_treedepth = 15), 
    open_progress = FALSE
  )
} else if (length(temp$dm != 4 * length(unique(temp$Study)))) {
  print("WITH CONDITION RANDOM EFFECTS")      
  test <- brm(
    formula = dm ~ 1 + dominance + language_transition + (Condition | Study),
    data = temp,
    family = gaussian(),
    inits = 0,
    prior = c(
      set_prior("student_t(3, 847.5, 215)", class = "Intercept"),
      set_prior("normal(0,10)", class = "b"),
      set_prior("student_t(3, 0, 215)", class = "sd")
    ),
    warmup = 1000,
    iter = 5000,
    chains = 4,
    silent = 2,
    refresh = 0,
    cores = 4,
    control = list(adapt_delta = 0.95, max_treedepth = 15)
  )
}
```

    [1] "WITHOUT CONDITION RANDOM EFFECTS"

    Warning: Argument 'inits' is deprecated. Please use argument 'init' instead.

``` r
modelwointeraction <- test

# fit criteria

modelwointeraction <- add_criterion(
  modelwointeraction,
  c("loo", "waic", "bayes_R2")
)
```

    Warning: 
    38 (4.0%) p_waic estimates greater than 0.4. We recommend trying loo instead.

``` r
modelwointeraction$criteria$waic
```


    Computed from 16000 by 956 log-likelihood matrix.

              Estimate   SE
    elpd_waic  -5728.9 30.9
    p_waic        63.9  4.6
    waic       11457.8 61.8

    38 (4.0%) p_waic estimates greater than 0.4. We recommend trying loo instead. 

``` r
modelwointeraction$criteria$loo
```


    Computed from 16000 by 956 log-likelihood matrix.

             Estimate   SE
    elpd_loo  -5729.6 31.0
    p_loo        64.6  4.7
    looic     11459.2 61.9
    ------
    MCSE of elpd_loo is 0.1.
    MCSE and ESS estimates assume MCMC draws (r_eff in [0.4, 2.5]).

    All Pareto k estimates are good (k < 0.7).
    See help('pareto-k-diagnostic') for details.

``` r
mean(modelwointeraction$criteria$bayes_R2)
```

    [1] 0.8104332

``` r
#model fit
pp_check(modelwinteraction, type = "scatter_avg") +
  geom_abline(intercept = 0, slope = 1, color = "red", lty = 4) # y vs. yrep
```

    Using all posterior draws for ppc type 'scatter_avg' by default.

![](gade_files/figure-commonmark/unnamed-chunk-3-1.png)

### Model comparison

``` r
loo(modelwinteraction, modelwointeraction)
```

    Output of model 'modelwinteraction':

    Computed from 16000 by 956 log-likelihood matrix.

             Estimate   SE
    elpd_loo  -5730.2 31.0
    p_loo        65.8  4.8
    looic     11460.5 62.0
    ------
    MCSE of elpd_loo is 0.1.
    MCSE and ESS estimates assume MCMC draws (r_eff in [0.4, 2.0]).

    All Pareto k estimates are good (k < 0.7).
    See help('pareto-k-diagnostic') for details.

    Output of model 'modelwointeraction':

    Computed from 16000 by 956 log-likelihood matrix.

             Estimate   SE
    elpd_loo  -5729.6 31.0
    p_loo        64.6  4.7
    looic     11459.2 61.9
    ------
    MCSE of elpd_loo is 0.1.
    MCSE and ESS estimates assume MCMC draws (r_eff in [0.4, 2.5]).

    All Pareto k estimates are good (k < 0.7).
    See help('pareto-k-diagnostic') for details.

    Model comparisons:
                       elpd_diff se_diff
    modelwointeraction  0.0       0.0   
    modelwinteraction  -0.6       1.1   

### Fixed effects

``` r
fixef(modelwinteraction)
```

                                    Estimate Est.Error       Q2.5      Q97.5
    Intercept                     894.926160 22.393281 851.478159 939.078397
    dominance                      -4.421014  2.899920 -10.110029   1.213794
    language_transition            26.996371  2.896898  21.288224  32.595299
    dominance:language_transition  -3.102481  2.886172  -8.785487   2.476238

``` r
fixef(modelwointeraction)
```

                          Estimate Est.Error      Q2.5      Q97.5
    Intercept           892.961124 23.936248 846.71888 941.401027
    dominance            -4.429385  2.873447 -10.07973   1.144334
    language_transition  27.001267  2.900509  21.25408  32.707197

## Summary

As in the paper the model without an interaction term fit better. And
the credible intervals for switch cost asymmetry and language dominance
cross zero.

But as the output shows, there is a bug in the control flow. We can see
that the “WITHOUT CONDITION RANDOM EFFECTS” is printed both times. In
other words, the model without condition random effects was run.

The issue is in the `if` statements. The below code appears to be
written to check if the length of the `dm` (RT) column equals the number
of unique studies. If these don’t match, it would mean that studies
contribute multiple effects.

However, a misplaced parenthesis means that what it *actually* does is
compare each value of `dm` to the number of unique studies, and then
returns the length of that vector. In this case, 956:

``` r
length(temp$dm != 4 * length(unique(temp$Study)))
```

    [1] 956

The correct version is this, which checks whether there are more RTs
than studies. It evaluates to `TRUE` as we’d expect.

``` r
length(temp$dm) > 4 * length(unique(temp$Study))
```

    [1] TRUE

In other words, none of these models had random effects for `condition`,
which meant that between-condition variability went into the residual.
This is a problem not only for the models themselves, but for the paper,
which states that there *were* condition random effects.

# After correccting

Here, I re-run the analyses by removing the `if` statement since we know
that studies contribute multiple conditions. Note that I do not think
this is the best way to analyse the data either, since effects are not
weighted by their precision as they should be in a meta-analysis. This
re-analysis is just to show the consequence of the bug.

The original code added random effects per condition with
`(Condition|Study)` but this doesn’t match the article’s description for
conditions *within* studies.

So I have changed the relevant code to match, adding the following
terms: `(1 | Study) + (1 | Study:Condition)`.

## Model with interaction

``` r
  test <- brm(
    formula = dm ~ 1 +
      dominance +
      language_transition +
      dominance * language_transition +
      (1 | Study) + (1 | Study:Condition),
    data = temp,
    family = gaussian(),
    inits = 0,

    prior = c(
      set_prior("student_t(3, 847.5, 215)", class = "Intercept"),
      set_prior("normal(0,10)", class = "b"),
      set_prior("student_t(3, 0, 215)", class = "sd")
    ),
    warmup = 1000,
    iter = 5000,
    chains = 4,
    silent = 2,
    refresh = 0,
    control = list(adapt_delta = 0.95, max_treedepth = 15), 
    open_progress = FALSE
  )
```

    Warning: Argument 'inits' is deprecated. Please use argument 'init' instead.

``` r
modelwinteraction <- test

# fit criteria

modelwinteraction <- add_criterion(
  modelwinteraction,
  c("loo", "waic", "bayes_R2")
)
```

    Warning: Found 7 observations with a pareto_k > 0.7 in model
    'modelwinteraction'. We recommend to set 'moment_match = TRUE' in order to
    perform moment matching for problematic observations.

    Warning: 
    118 (12.3%) p_waic estimates greater than 0.4. We recommend trying loo instead.

``` r
modelwinteraction$criteria$waic
```


    Computed from 16000 by 956 log-likelihood matrix.

              Estimate   SE
    elpd_waic  -4986.2 37.5
    p_waic       209.1 14.9
    waic        9972.5 75.0

    118 (12.3%) p_waic estimates greater than 0.4. We recommend trying loo instead. 

``` r
modelwinteraction$criteria$loo
```


    Computed from 16000 by 956 log-likelihood matrix.

             Estimate   SE
    elpd_loo  -4998.1 38.4
    p_loo       221.0 15.7
    looic      9996.2 76.7
    ------
    MCSE of elpd_loo is NA.
    MCSE and ESS estimates assume MCMC draws (r_eff in [0.4, 1.1]).

    Pareto k diagnostic values:
                             Count Pct.    Min. ESS
    (-Inf, 0.7]   (good)     949   99.3%   311     
       (0.7, 1]   (bad)        5    0.5%   <NA>    
       (1, Inf)   (very bad)   2    0.2%   <NA>    
    See help('pareto-k-diagnostic') for details.

``` r
mean(modelwinteraction$criteria$bayes_R2)
```

    [1] 0.9668558

``` r
#model fit
pp_check(modelwinteraction, type = "scatter_avg") +
  geom_abline(intercept = 0, slope = 1, color = "red", lty = 4) # y vs. yrep
```

    Using all posterior draws for ppc type 'scatter_avg' by default.

![](gade_files/figure-commonmark/unnamed-chunk-9-1.png)

## Model without interaction

``` r
  test <- brm(
    formula = dm ~ 1 + dominance + language_transition +
      (1 | Study) + (1 | Study:Condition),
    data = temp,
    family = gaussian(),
    inits = 0,
    prior = c(
      set_prior("student_t(3, 847.5, 215)", class = "Intercept"),
      set_prior("normal(0,10)", class = "b"),
      set_prior("student_t(3, 0, 215)", class = "sd")
    ),
    warmup = 1000,
    iter = 5000,
    chains = 4,
    silent = 2,
    refresh = 0,
    control = list(adapt_delta = 0.95, max_treedepth = 15), 
    open_progress = FALSE
  )
```

    Warning: Argument 'inits' is deprecated. Please use argument 'init' instead.

``` r
modelwointeraction <- test

# fit criteria

modelwointeraction <- add_criterion(
  modelwointeraction,
  c("loo", "waic", "bayes_R2")
)
```

    Warning: Found 5 observations with a pareto_k > 0.7 in model
    'modelwointeraction'. We recommend to set 'moment_match = TRUE' in order to
    perform moment matching for problematic observations.

    Warning: 
    124 (13.0%) p_waic estimates greater than 0.4. We recommend trying loo instead.

``` r
modelwointeraction$criteria$waic
```


    Computed from 16000 by 956 log-likelihood matrix.

              Estimate   SE
    elpd_waic  -4990.2 37.3
    p_waic       208.4 14.8
    waic        9980.5 74.7

    124 (13.0%) p_waic estimates greater than 0.4. We recommend trying loo instead. 

``` r
modelwointeraction$criteria$loo
```


    Computed from 16000 by 956 log-likelihood matrix.

             Estimate   SE
    elpd_loo  -5001.9 38.0
    p_loo       220.0 15.4
    looic     10003.7 75.9
    ------
    MCSE of elpd_loo is NA.
    MCSE and ESS estimates assume MCMC draws (r_eff in [0.4, 1.1]).

    Pareto k diagnostic values:
                             Count Pct.    Min. ESS
    (-Inf, 0.7]   (good)     951   99.5%   295     
       (0.7, 1]   (bad)        4    0.4%   <NA>    
       (1, Inf)   (very bad)   1    0.1%   <NA>    
    See help('pareto-k-diagnostic') for details.

``` r
mean(modelwointeraction$criteria$bayes_R2)
```

    [1] 0.9665481

``` r
#model fit
pp_check(modelwinteraction, type = "scatter_avg") +
  geom_abline(intercept = 0, slope = 1, color = "red", lty = 4) # y vs. yrep
```

    Using all posterior draws for ppc type 'scatter_avg' by default.

![](gade_files/figure-commonmark/unnamed-chunk-10-1.png)

### Model comparison

``` r
loo(modelwinteraction, modelwointeraction)
```

    Output of model 'modelwinteraction':

    Computed from 16000 by 956 log-likelihood matrix.

             Estimate   SE
    elpd_loo  -4998.1 38.4
    p_loo       221.0 15.7
    looic      9996.2 76.7
    ------
    MCSE of elpd_loo is NA.
    MCSE and ESS estimates assume MCMC draws (r_eff in [0.4, 1.1]).

    Pareto k diagnostic values:
                             Count Pct.    Min. ESS
    (-Inf, 0.7]   (good)     949   99.3%   311     
       (0.7, 1]   (bad)        5    0.5%   <NA>    
       (1, Inf)   (very bad)   2    0.2%   <NA>    
    See help('pareto-k-diagnostic') for details.

    Output of model 'modelwointeraction':

    Computed from 16000 by 956 log-likelihood matrix.

             Estimate   SE
    elpd_loo  -5001.9 38.0
    p_loo       220.0 15.4
    looic     10003.7 75.9
    ------
    MCSE of elpd_loo is NA.
    MCSE and ESS estimates assume MCMC draws (r_eff in [0.4, 1.1]).

    Pareto k diagnostic values:
                             Count Pct.    Min. ESS
    (-Inf, 0.7]   (good)     951   99.5%   295     
       (0.7, 1]   (bad)        4    0.4%   <NA>    
       (1, Inf)   (very bad)   1    0.1%   <NA>    
    See help('pareto-k-diagnostic') for details.

    Model comparisons:
                       elpd_diff se_diff
    modelwinteraction   0.0       0.0   
    modelwointeraction -3.8       3.1   

### Fixed effects

``` r
fixef(modelwinteraction)
```

                                    Estimate Est.Error       Q2.5       Q97.5
    Intercept                     891.406415 23.622617 846.545321 938.7722999
    dominance                      -4.744846  1.251358  -7.197089  -2.2793836
    language_transition            29.006128  1.261203  26.539962  31.4649575
    dominance:language_transition  -3.333520  1.256104  -5.842025  -0.8469377

``` r
fixef(modelwointeraction)
```

                          Estimate Est.Error       Q2.5      Q97.5
    Intercept           891.503094 23.608871 844.803875 938.298371
    dominance            -4.766423  1.251439  -7.183035  -2.355282
    language_transition  29.004656  1.257380  26.537596  31.455579

``` r
library(emmeans)
```

    Welcome to emmeans.
    Caution: You lose important information if you filter this package's results.
    See '? untidy'

``` r
grid <- list(dominance = c(-1, 1), language_transition = c(-1, 1))

emm <- emmeans(modelwinteraction, ~ language_transition | dominance, at = grid)
emm
```

    dominance = -1:
     language_transition emmean lower.HPD upper.HPD
                      -1    864       819       911
                       1    928       883       976

    dominance =  1:
     language_transition emmean lower.HPD upper.HPD
                      -1    861       815       908
                       1    912       865       958

    Point estimate displayed: median 
    HPD interval probability: 0.95 

``` r
sc <- contrast(emm, method = list(switch_cost = c(-1, 1)))
sc
```

    dominance = -1:
     contrast    estimate lower.HPD upper.HPD
     switch_cost     64.7      57.6      71.4

    dominance =  1:
     contrast    estimate lower.HPD upper.HPD
     switch_cost     51.3      44.5      58.7

    Point estimate displayed: median 
    HPD interval probability: 0.95 

``` r
contrast(sc, method = list(asymmetry = c(1, -1)), by = NULL)
```

     contrast  estimate lower.HPD upper.HPD
     asymmetry     13.3      3.69      23.6

    Point estimate displayed: median 
    HPD interval probability: 0.95 

### Summary

The fit flips: the model *with* interaction fits better now. And the
relevant coefficients’ credible intervals no longer include zero. The
posterior predictive checks look a lot better, too.

The coefficients are, of course, not the actual RT differences because
of the -1/+1 coding. The emmeans output suggests that the asymmetry is
around 13ms.
