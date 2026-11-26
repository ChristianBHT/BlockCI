# CausalBlocks: causal discovery between blocks of data

R package for causal discovery between **blocks** of data. A block is a matrix. Two settings are covered:

- the same samples with different variables (e.g. several measurement blocks on the same objects),
- the same variables with different samples (environments).

## Aims

1. A viable test for **conditional independence (CI) between blocks** of data (matrices).
2. A **causal discovery** method built on that CI test.
3. An **R package** that makes both usable.

## Status

Early development. The current work is running simulations and tinkering with the
method to get the best possible testing results (type I error, power, speed). Expect
changes in function names, arguments and defaults. Preliminary simulation results are
summarised [below](#preliminary-simulation-results); details are in `LOG.md`.

The idea and the method are the developers' own. The code is written mainly with
Claude Code and is studied and checked by the developers.

## Installation

```r
# from the package folder
devtools::install()            # or devtools::load_all() during development
library(CausalBlocks)
```

Required: `multiblock`, `pls`. Optional: `grf` (`method = "RF"`), `xgboost` >= 2.0
(`method = "XGB"`), `FNN` (`method = "KNN"`), `energy` (`stat = "dcor"`).

## The tests in short

All tests take a formula `Y ~ X | Z1 + Z2 + ...` and a named list of matrices, and test
H0: Y is independent of X given (Z1, Z2, ...). Each name after the bar is one
conditioning block.

| Function | Idea | p-value from |
|---|---|---|
| `Block.CCI.test()` | Does X improve the held-out prediction of Y beyond Z? `T = RMSEP(Y ~ Z + X) - RMSEP(Y ~ Z)` | `nperm` refits on a resampled X (rotation, permutation or conditional permutation) |
| `Block.CCI.fast()` | The same statistic, one fit per split instead of a Monte Carlo loop | resampled t-test over `nsplit` splits (experimental) |
| `Block.resid.test()` | Remove Z from X and from Y by regression; are the residual matrices still dependent? | permutation (RV, first canonical correlation, distance correlation) or a Gaussian approximation (weighted GCM) |

The regression is chosen with `method`: `"PLS"`, `"SOPLS"`, `"ROSA"`, `"RF"`
(multi-output random forest, grf), `"XGB"` (multi-target XGBoost) or `"KNN"` (k
nearest neighbours; new, not part of the simulation results below). In the residual
test RF and XGB first remove the linear part of Z with PLS (`linear_first = TRUE`).
Many conditioning variables can be reduced first with `reduce = "pca"` or `"pls"`.

## Testing under different scenarios

The examples use one small generator so that each scenario can be run as it stands.
X and Y depend only on the *observed* Z blocks plus their own noise, so with
`effect = 0` the null holds exactly; `effect > 0` adds a direct effect X -> Y.

```r
make_blocks <- function(n, design = c("linear", "nonlinear", "highdim"), effect = 0,
                        q = 4, hetero = FALSE) {
  design <- match.arg(design)
  if (design == "highdim") {                     # 200 Z variables driven by 3 factors
    L <- matrix(rnorm(3 * 200), 3, 200)
    Z <- matrix(rnorm(n * 3), n, 3) %*% L + matrix(rnorm(n * 200), n, 200)
    S <- scale(Z %*% t(L))                       # the 3 directions of the observed Z
  } else {
    Z <- matrix(rnorm(n * 8), n, 8)              # Z1, Z2 with 4 columns each
    S <- Z
  }
  Fm <- if (design == "nonlinear") scale(cbind(S, S^2, sin(2 * S), S[, 1] * S[, 2])) else S
  mu <- function() Fm %*% matrix(rnorm(ncol(Fm) * q, sd = sqrt(2 / ncol(Fm))), ncol(Fm), q)
  sd_e <- if (hetero) sqrt(0.2 + 0.8 * as.numeric(scale(Z %*% rnorm(ncol(Z))))^2) else 1
  X <- mu() + matrix(rnorm(n * q), n, q) * sd_e
  Y <- mu() + matrix(rnorm(n * q), n, q) * sd_e + effect * X
  h <- ncol(Z) / 2
  list(Y = Y, X = X, Z1 = Z[, 1:h], Z2 = Z[, -(1:h)])
}
```

### 1. Linear relations

Every test works here. PLS is the fastest and has the best power of the CCI tests;
the residual test with PLS is faster still.

```r
set.seed(1)
d0 <- make_blocks(500, "linear")               # H0 true
d1 <- make_blocks(500, "linear", effect = 0.5) # H0 false

Block.CCI.test(Y ~ X | Z1 + Z2, d0, method = "PLS", nperm = 99)   # large p-value
Block.CCI.test(Y ~ X | Z1 + Z2, d1, method = "PLS", nperm = 99)   # small p-value
Block.CCI.test(Y ~ X | Z1, d0, method = "PLS", nperm = 99)        # Z2 missing: H0 false

Block.resid.test(Y ~ X | Z1 + Z2, d0, method = "PLS")
Block.resid.test(Y ~ X | Z1 + Z2, d1, method = "PLS")
```

### 2. Noise whose level depends on Z (heteroskedastic)

Permutation tests of the residual RV coefficient reject true nulls far too often
here (about 0.5 in our simulations). The weighted GCM statistic (`stat = "wgcm"`) and
the CCI test hold the level.

```r
dh <- make_blocks(500, "linear", hetero = TRUE)
Block.resid.test(Y ~ X | Z1 + Z2, dh, method = "PLS")                 # too liberal
Block.resid.test(Y ~ X | Z1 + Z2, dh, method = "PLS", stat = "wgcm")  # holds the level
Block.CCI.test(Y ~ X | Z1 + Z2, dh, method = "PLS", nperm = 99)
```

### 3. Nonlinear relations

A linear regression cannot remove a nonlinear effect of Z, so the linear methods reject
true nulls. Use a nonlinear learner in the CCI test. The residual tests are not
reliable here, even with a random forest: the part of the shared nonlinear signal that
the forest does not learn stays in both residual matrices.

```r
dn <- make_blocks(500, "nonlinear")
Block.CCI.test(Y ~ X | Z1 + Z2, dn, method = "RF", nperm = 99,
               method_args = list(num.trees = 500, num.threads = NULL))   # recommended
Block.resid.test(Y ~ X | Z1 + Z2, dn, method = "PLS")                      # rejects H0 wrongly
Block.resid.test(Y ~ X | Z1 + Z2, dn, method = "RF", stat = "wgcm")        # often too liberal
```

CCI with RF holds the level up to n = 500 and becomes somewhat liberal at n = 1000
(see the results below). Keep the default `resampling = "rotation"`: in our first runs
permutation and conditional permutation made CCI with RF strongly liberal.

### 4. Many conditioning variables (high-dimensional Z)

With 200 Z variables and linear relations the residual test with PLS is fast, holds the
level and has full power. The CCI test needs n well above the number of Z variables:
its null distribution rotates X against all of Z, and when Z has about as many columns
as there are rows the rotated X equals X and the test has no power.

```r
dz0 <- make_blocks(500, "highdim")           # linear in the 3 factor directions
dz1 <- make_blocks(500, "highdim", effect = 0.5)
Block.resid.test(Y ~ X | Z1 + Z2, dz0, method = "PLS")
Block.resid.test(Y ~ X | Z1 + Z2, dz1, method = "PLS")
Block.CCI.test(Y ~ X | Z1 + Z2, dz1, method = "PLS", reduce = "pca", nperm = 99)
```

`reduce = "pca"` replaces Z by its principal component scores (the components above the
Marchenko-Pastur noise edge are kept with `reduce_ncomp = "auto"`); `reduce = "pls"`
uses PLS scores (in the residual test: RF or XGB on the scores of the PLS step).

### 5. Nonlinear and high-dimensional: a limit case

When X and Y share a nonlinear function of a few directions hidden among many noisy Z
variables, no regression-based test can be trusted at moderate n (Shah and Peters
2020): the regression cannot be learned well enough. In our simulations every residual
test rejects almost all true nulls, also with dimension reduction, and CCI with RF on
reduced Z holds the level only at n = 200 (0.11-0.21 at n = 500-1000). Treat results in
this setting with caution.

### Which test when

| Situation | Suggested test |
|---|---|
| Linear relations, homoskedastic noise | `Block.resid.test(method = "PLS")` (fast) or `Block.CCI.test(method = "PLS")` |
| Noise level depends on Z | `Block.resid.test(stat = "wgcm")` or `Block.CCI.test(method = "PLS")` |
| Nonlinear relations | `Block.CCI.test(method = "RF")`, rotation |
| Many Z variables, linear | `Block.resid.test(method = "PLS")`; CCI only when n is well above ncol(Z) |
| Many Z variables, nonlinear | no reliable test yet; CCI with RF and `reduce` at small n |

## Preliminary simulation results

Script: `inst/benchmarks/simulation_study_noise.R` (one function per scenario, e.g.
`study_t3_noise(c("CCI_PLS", "CCI_RF"), effect = 0, n = 500, n_tests = 100)`).
Blocks: Y and X with 4 columns, Z1 and Z2 with 4 columns each (high-dimensional: X and
Y with 10 columns, Z1 and Z2 with 100 each). Each X and Y column has variance 2
explained by Z plus noise of variance 1. Noise types: Gaussian, t3, Laplace, skewed
(centred exponential), Poisson counts, 5 % outliers, heteroskedastic (sd depends on
Z), correlated columns. Power: `Y <- Y + 0.5 * X`. 100 datasets per cell, level 0.05:
the Monte Carlo standard error of a rejection rate near 0.05 is about 0.02, so rates
above about 0.09 are clearly too high. CCI tests: `nperm = 99`, RF with 500 trees.
All results here predate two changes: the CCI tests used a training fraction
`p = 0.75` (now 0.5), and conditional permutation, also the default of the residual
tests with RV, was within k-means groups (now local kNN).
Residual tests: 5 folds, `B = 199`.

**Linear design, type I error** (median over the 8 noise types, worst in brackets)

| Test | n = 200 | n = 500 | n = 1000 | Power (median), n = 200 / 500 / 1000 |
|---|---|---|---|---|
| CCI, PLS | 0.03 (0.05) | 0.04 (0.08) | 0.04 (0.09) | 0.99 / 1.00 / 1.00 |
| CCI, ROSA | 0.04 (0.08) | 0.03 (0.13) | 0.05 (0.08) | 0.94 / 1.00 / 1.00 |
| CCI, SO-PLS (4 noise types) | 0.04 (0.05) | 0.07 (0.07) | 0.05 (0.07) | 0.47 / 0.54 / 0.64 |
| CCI, RF | 0.05 (0.08) | 0.07 (0.11) | 0.08 (0.17) | 0.56 / 0.85 / 0.97 |
| Residual, PLS, RV | 0.08 (0.42) | 0.08 (0.51) | 0.06 (0.55) | 1.00 / 1.00 / 1.00 |
| Residual, PLS + RF, RV | 0.11 (0.41) | 0.06 (0.53) | 0.07 (0.56) | 1.00 / 1.00 / 1.00 |
| Residual, PLS + RF, wGCM | 0.04 (0.07) | 0.03 (0.07) | 0.05 (0.08) | 1.00 / 1.00 / 1.00 |

The worst cases for the RV residual tests are all heteroskedastic noise; for CCI with
RF they are skewed noise.

**Nonlinear design, type I error** (median over the 8 noise types, worst in brackets)

| Test | n = 200 | n = 500 | n = 1000 | Power (median), n = 200 / 500 / 1000 |
|---|---|---|---|---|
| CCI, RF | 0.07 (0.11) | 0.09 (0.15) | 0.16 (0.25) | 0.84 / 0.96 / 1.00 |
| Residual, PLS, RV (Gaussian only) | 0.73 | 0.99 | 1.00 | – |
| Residual, PLS + RF, RV | 0.71 (0.82) | 0.92 (0.95) | 0.99 (1.00) | – |
| Residual, PLS + RF, wGCM | 0.18 (0.30) | 0.49 (0.53) | 0.73 (0.85) | – |

**High-dimensional designs** (Gaussian noise, 200 Z variables from a 3-factor model)

| Design, test | Type I, n = 200 / 500 / 1000 | Power, n = 200 / 500 / 1000 |
|---|---|---|
| linear, residual PLS | 0.06 / 0.06 / 0.02 | 1.00 / 1.00 / 1.00 |
| linear, residual PLS on PCA scores | 0.11 / 0.06 / 0.02 | 1.00 / 1.00 / 1.00 |
| linear, CCI PLS | 0.06 / 0.04 / 0.03 | 0.03 / 0.59 / 0.99 |
| linear, CCI PLS on PCA scores | 0.04 / 0.05 / 0.01 | 0.04 / 1.00 / 1.00 |
| nonlinear, CCI RF on PCA scores | 0.06 / 0.13 / 0.14 | 0.03 / 0.61 / 0.91 |
| nonlinear, CCI RF on PLS scores | 0.05 / 0.11 / 0.21 | 0.01 / 0.78 / 0.97 |
| nonlinear, residual RF (PCA or PLS scores), RV | 1.00 / 1.00 / 1.00 | – |
| nonlinear, residual RF on PCA scores, wGCM | 0.61 / 0.88 / 0.97 | – |

The low CCI power at n = 200 comes from the rotation: with 200 Z variables and 200 rows
the rotated X equals X.

**Resampling of X in the CCI test** (first results: linear design, Gaussian noise, type I error)

| Test | rotation, n = 500 / 1000 | permutation | conditional permutation |
|---|---|---|---|
| CCI, PLS | 0.04 / 0.02 | 0.04 / 0.01 | 0.02 / 0.01 |
| CCI, RF | 0.06 / 0.08 | 0.29 / 0.42 | 0.33 / 0.46 |

Permuting X breaks its relation to Z. Under the null the observed X still carries
information about Z that the forest has not fully learned, while the permuted X does
not, so the null distribution is shifted and the test rejects too often. Conditional
permutation (then within k-means groups of Z; since replaced by a local kNN
permutation, not yet studied) did not help here. The remaining noise types
and the nonlinear design are running.

**Time per test** (median, one test, all cores for RF): CCI with PLS 1-3 s, ROSA
5-12 s, RF 14-29 s, SO-PLS 4-46 s; residual tests under 2 s.

**Summary so far.** With linear relations all tests work, and the residual test with
the weighted GCM statistic is the fastest test that holds the level under every noise
type. With nonlinear relations only the CCI test with a random forest comes close to
the nominal level; it degrades slowly as n grows, while the residual tests break down.
High-dimensional nonlinear conditioning remains an open problem.

## Block CCI test: details

`Block.CCI.test(Y ~ X | V, data, method = "PLS")` is a block version of the
Computational Conditional Independence (CCI) test, with the same structure as the CCI
package.

- **Statistic:** `T = RMSEP(Y ~ V + X) - RMSEP(Y ~ V)` on held-out data (a fraction
  `p` of the rows is used for training). Negative `T` means X improves prediction
  beyond V.
- **Regression (`method`):** `"PLS"` (default, `pls::plsr` on the cbind-ed blocks),
  `"SOPLS"` (`multiblock::sopls`, each V block kept as a separate block), `"ROSA"`
  (`multiblock::rosa`, every block kept separate), `"XGB"` (multi-target XGBoost,
  experimental), `"RF"` (multi-output random forest, `grf::multi_regression_forest`) or
  `"KNN"` (the mean Y of the k nearest training rows, k chosen by inner CV or set with
  `method_args = list(k = 10)`; Euclidean distance on standardised columns).
  Method-specific options go in `method_args`.
- **Null distribution:** `nperm` Monte Carlo CV instances on a resampled X: rotation
  (`resampling = "rotation"`, default; keeps the linear relation between X and V),
  permutation (`"permutation"`) or the greedy local permutation of Runge (2018)
  (`"conditional"`): each row gets the X of one of its `k_perm` (default 5) nearest
  neighbours in V, using each neighbour once where possible.
- **Components (SOPLS, PLS, ROSA):** chosen by inner cross-validation on the training
  rows (`ncomp = "auto"`) or given. For SOPLS `c(a_1, ..., a_k, b)`: components for
  each of the k V blocks, then b for X. For PLS and ROSA `c(a, b)`: components of
  Y ~ V and Y ~ V + X.
- **p-value:** share of null statistics at or below the observed one
  (`pvalue = "empirical"`), or from a normal fit to the null statistics
  (`pvalue = "parametric"`).
- **Observed splits (`npval`, default 1):** the observed statistic depends on the
  random train/test split. With `npval > 1` (a whole number up to `nperm`,
  experimental) it is computed on `npval` splits, each is compared with the same null
  statistics and the p-values are averaged. The p-value then depends less on the
  split, but the average is not uniform under the null: in linear Gaussian
  simulations with PLS, 2 or more made the test far too conservative and lowered its
  power at level 0.05. Keep the default for testing. The simulation results above
  were made with a single split.
- **Cross-fitting (`crossfit = TRUE`, experimental, off by default):** the models are
  fitted on each half of the rows and predict the other, one statistic is computed
  from the errors of all rows, and every null draw is cross-fitted the same way
  (about twice the fits). With linear Gaussian data and PLS (n = 500) it held the
  level (0.039) and had more power than one 50/50 split (0.30 against 0.24 at effect
  0.10, 0.68 against 0.57 at 0.15; a second run at 0.10 showed a smaller gain, 0.235
  against 0.208, with level 0.038). With a random forest it rejected true nulls too
  often (0.10-0.13), so do not use it with RF on the present evidence.

```r
set.seed(1)
data <- sim_1(n = 500, p = 4)   # list(Y, X, Z1, Z2); Y _||_ X | Z1, Z2
res <- Block.CCI.test(Y ~ X | Z1 + Z2, data, nperm = 99)
res
plot(res)                        # null distribution with the observed statistic
Block.CCI.test(Y ~ X | Z1 + Z2, data, nperm = 99, pvalue = "parametric")
Block.CCI.test(Y ~ X | Z1 + Z2, data, nperm = 99, method = "XGB",
               method_args = list(nrounds = 300, params = list(max_depth = 4, eta = 0.05)))
```

With few null draws the empirical p-value cannot go below `1/(nperm + 1)`; the
parametric p-value is finer, at the price of assuming normality. SO-PLS has low power in
our simulations (one component per V block by default). The data generators
`sim_1()` to `sim_5()` (in `R/scenarios.R`) give further null scenarios: latent common
causes (`sim_2()`), a chain (`sim_3()`), nearly collinear X (`sim_4()`) and a nonlinear,
heteroskedastic case (`sim_5()`).

**Adding a method.** A method is a *learner* function
`function(Ys, X, Vl, train, test, opts)` that fits Y ~ V and Y ~ V + X on the training
rows and returns `list(yh0, yh1, a, b)`: the test-row predictions of both models and
their complexity. The interface is described at the top of `R/learners.R`. A learner can
be passed directly as `method`, or made available by name by adding it to
`.cci_methods` in that file.

## Fast version

`Block.CCI.fast()` avoids the Monte Carlo loop. For each of `nsplit` train/test splits
it fits the model once on the real X and once on a rotated copy X* (the fit of Y ~ V is
reused), compares the two on the same held-out rows, and combines the splits with a
resampled t-test. It takes the same `method`, `resampling`, `method_args` and `reduce`
arguments as `Block.CCI.test()`.

```r
Block.CCI.fast(Y ~ X | Z1 + Z2, data)                                  # plain t across splits (default)
Block.CCI.fast(Y ~ X | Z1 + Z2, data, correction = "nadeau-bengio")    # valid but low power
Block.CCI.fast(Y ~ X | Z1, data, nsplit = 20, method = "PLS")
Block.CCI.fast(Y ~ X | Z1 + Z2, data, method = "RF",
               method_args = list(num.trees = 500, num.threads = NULL))
```

This version is experimental. In small simulations (n = 300, SO-PLS) the default (no
correction) had much better power than Nadeau-Bengio but rejected slightly too often
under the null (0.05-0.08); Nadeau-Bengio was valid but had low power. Use the Monte
Carlo test for a final p-value until the fast version has been calibrated.

## Residual test: details

`Block.resid.test(Y ~ X | Z, data, method = "PLS")` removes Z from X and from Y with a
regression and tests whether the two residual matrices are still dependent.

```r
Block.resid.test(Y ~ X | Z1 + Z2, data)                    # PLS (default), RV, conditional permutation
Block.resid.test(Y ~ X | Z1 + Z2, data, method = "ROSA")   # ROSA, Z blocks kept separate
Block.resid.test(Y ~ X | Z1 + Z2, data, method = "RF",     # PLS, then multi-output random forest
                 method_args = list(num.trees = 1000, num.threads = NULL))
Block.resid.test(Y ~ X | Z1 + Z2, data, stat = "wgcm")     # weighted GCM
```

- **Residuals:** with `resid = "cv"` (default) every row's residual comes from a model
  that did not use that row (`nfolds` folds, shared by the X and the Y regression). PLS
  and ROSA choose the number of components by the CV error over the same folds.
- **PLS first for RF and XGB** (`linear_first = TRUE`, default): in each fold a PLS
  regression on the training rows, then the tree learner on its training residuals.
  Tree ensembles approximate linear effects poorly; without this step a random forest
  left part of a linear effect of Z in both residual matrices and rejected a true null
  in every run of a small test.
- **Statistics:** `stat = "rv"` (RV coefficient), `"cc1"` (first canonical
  correlation), `"dcor"` (distance correlation), calibrated by permutation
  (`perm = "conditional"`, default: the greedy local permutation among the `k_perm`
  nearest neighbours in Z, see `permute_X_cond()`); or
  `"wgcm"`, the multivariate weighted generalised covariance measure with fixed weights
  (mWGCM.fix; Scheidegger, Hörrmann and Bühlmann 2022, JMLR 23, Appendix A.1): for every
  pair of residual columns and every weight `w_k(Z)` (the constant and
  `sign(Z_d - quantile)` at `k0` quantiles per column of Z) a studentised mean of the
  weighted products, the maximum of their absolute values, and a p-value from a
  Gaussian approximation (`B` draws).
- **Dimension reduction:** `reduce = "pca"` (Z replaced by PCA scores, fitted on all
  rows) or `"pls"` (RF / XGB on the PLS scores of the linear step).
- **Adding a method:** a *residual learner* `function(M, Zl, folds, opts)` returning
  `list(res, complexity)`; pass it as `method` or add it to `.resid_methods` in
  `R/resid_learners.R`.

`sim_scale()` simulates blocks of any size. In earlier simulations
(`inst/benchmarks/scale_timing.R`, `scale_power.R`) the PLS residual test handled one
million rows in about 5 to 7 minutes and thousands of variables in under a minute.

## Development

Open `BlockCI.Rproj` in RStudio.
`devtools::document()` rebuilds docs and NAMESPACE (roxygen2), `devtools::test()` runs tests,
`devtools::install()` installs.

- `R/Block.CCI.test.R`: the block CCI test (Monte Carlo)
- `R/cci_test_gen.R`, `R/cci_null_dist.R`: one statistic and the null distribution
- `R/learners.R`: the regression methods (SOPLS, PLS, ROSA, XGB, RF) and how to add one
- `R/fast.R`: the fast resampled-t version (`Block.CCI.fast()`)
- `R/Block.resid.test.R`, `R/resid_learners.R`: the residual test and its regression methods
- `R/reduce.R`: dimension reduction of the conditioning blocks
- `R/utils.R`: resampling of X, statistics, weighted GCM
- `R/plot.R`: plot method for test results
- `R/scenarios.R`, `R/sim_scale.R`: data generators
- `inst/benchmarks/simulation_study_noise.R`: the simulation studies above
- `tests/testthat/`: unit tests
- `LOG.md`: project log
