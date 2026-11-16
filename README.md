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
changes in function names, arguments and defaults. Results so far are in `LOG.md`.

The idea and the method are the developers' own. The code is written mainly with
Claude Code and is studied and checked by the developers.

## Block CCI test

`Block.CCI.test(Y ~ X | V, data, method = "SOPLS")` tests whether Y is independent of
X given V. It is a block version of the Computational Conditional Independence (CCI)
test, with the same structure as the CCI package.

- **Statistic:** `T = RMSEP(Y ~ V + X) - RMSEP(Y ~ V)` on held-out data. Negative `T`
  means X improves prediction beyond V.
- **Regression (`method`):** `"SOPLS"` (default, `multiblock::sopls`, each V block kept
  as a separate block), `"PLS"` (`pls::plsr` on the cbind-ed blocks), `"XGB"`
  (multi-target XGBoost, experimental) or `"RF"` (multi-output random forest,
  `grf::multi_regression_forest`). Method-specific options go in `method_args`.
- **Null distribution:** R Monte Carlo CV instances on a resampled X: rotation
  (`resampling = "rotation"`, default), permutation (`"permutation"`) or permutation
  within clusters of similar V (`"conditional"`).
- **Components (SOPLS, PLS):** chosen by inner cross-validation on the training rows
  (`ncomp = "auto"`) or given. For SOPLS `c(a_1, ..., a_k, b)`: components for each
  of the k V blocks, then b for X. For PLS `c(a, b)`: components of Y ~ V and Y ~ V + X.
  Several V blocks are written `Y ~ X | Z1 + Z2`.
- **p-value:** share of null statistics at or below the observed one
  (`pvalue = "empirical"`), or from a normal fit to the null statistics
  (`pvalue = "parametric"`).

## Examples

The data generators are in `R/scenarios.R` (`sim_1()` to `sim_5()` are exported). The
package needs the `multiblock` and `pls` packages; `xgboost` and `grf` are needed only
for `method = "XGB"` and `method = "RF"`:

```r
devtools::install()           # from the package folder, or devtools::load_all()
library(CausalBlocks)
```

### A simple null case

`sim_1()` makes four matrices Y, X, Z1 and Z2 where X and Y both depend on Z1
and Z2 and nothing else, so Y is independent of X given (Z1, Z2) but not given Z1 alone.

```r
set.seed(1)
data <- sim_1(n = 500, p = 4)   # list(Y, X, Z1, Z2), p variables per block

# Y _||_ X | Z1, Z2 is true: we expect a large p-value
res <- Block.CCI.test(Y ~ X | Z1 + Z2, data, nperm = 99)
res
plot(res)                        # null distribution with the observed statistic

# Y _||_ X | Z1 is false (Z2 is missing): we expect a small p-value
Block.CCI.test(Y ~ X | Z1, data, nperm = 99)

# the same with plain PLS
Block.CCI.test(Y ~ X | Z1, data, nperm = 99, method = "PLS")
```

### Parametric p-value

With few null draws the empirical p-value cannot go below `1/(R + 1)`. A normal
approximation of the null distribution gives a finer p-value, at the price of
assuming normality:

```r
res_par <- Block.CCI.test(Y ~ X | Z1 + Z2, data, nperm = 99, pvalue = "parametric")
c(empirical = res_par$p.empirical, parametric = res_par$p.parametric)
plot(res_par)                    # adds the fitted normal curve
```

### Scenario with mixed signals

`sim_2()` generates blocks Y (4 variables), X (4), Z1 (5) and Z2 (5) from latent
common causes. Y is independent of X given (Z1, Z2) by construction.

```r
d2 <- sim_2(n = 500)
Block.CCI.test(Y ~ X | Z1 + Z2, d2, nperm = 99)   # null is true
Block.CCI.test(Y ~ X | Z1,      d2, nperm = 99)   # null is false
```

`sim_3()` (chain `X <- Z1 -> Z2 -> Y`), `sim_4()` (nearly collinear X) and `sim_5()`
(nonlinear, heteroskedastic) work the same way. The nonlinear `sim_5()` is known to
fail for the linear methods: they reject a true null too often (see `LOG.md`).
`sim_6()` and `sim_7()` are still being worked on.

### Nonlinear methods (XGB, RF)

`method = "XGB"` uses multi-target XGBoost (`xgboost` >= 2.0,
`multi_strategy = "multi_output_tree"`). Each model is trained with a fixed number of
boosting rounds (`nrounds`, default 1000) on all the training rows, without early
stopping, and only predicts the held-out rows. It is experimental: in our simulations
it rejects true nulls too often on nonlinear data (see `LOG.md`).

`method = "RF"` uses a multi-output random forest (`grf::multi_regression_forest()`).
Each forest is grown on all the training rows with splits that reduce the error summed
over the Y columns. Options in `method_args` are passed to grf (defaults here:
`num.trees = 2000`, `num.threads = 1`). Not yet checked in simulations.

```r
Block.CCI.test(Y ~ X | Z1 + Z2, data, nperm = 49, method = "XGB",
               method_args = list(nrounds = 300, params = list(max_depth = 4, eta = 0.05)))
Block.CCI.test(Y ~ X | Z1 + Z2, data, nperm = 99, method = "RF",
               method_args = list(num.trees = 500, min.node.size = 10, num.threads = NULL))
```

The rotation or permutation of X does not condition on V nonlinearly, so on strongly
nonlinear scenarios these tests can still reject a true null too often.

### Adding a method

A method is a *learner* function `function(Ys, X, Vl, train, test, opts)` that fits
Y ~ V and Y ~ V + X on the training rows and returns `list(yh0, yh1, a, b)`: the
test-row predictions of both models and their complexity. The interface is described
at the top of `R/learners.R`. A learner can be passed directly as `method`, or made
available by name by adding it to `.cci_methods` in that file.

## Fast version

`Block.CCI.fast()` avoids the Monte Carlo loop. For each of `nsplit` train/test splits
it fits the model once on the real X and once on a rotated copy X*, compares the two
on the same held-out rows, and combines the splits with a resampled t-test. It takes
the same `method`, `resampling` and `method_args` arguments as `Block.CCI.test()`. With
SO-PLS it takes about 1-2 seconds where the Monte Carlo test takes about 10 seconds or
more.

```r
Block.CCI.fast(Y ~ X | Z1 + Z2, data)                                  # plain t across splits (default)
Block.CCI.fast(Y ~ X | Z1 + Z2, data, correction = "nadeau-bengio")    # valid but low power
Block.CCI.fast(Y ~ X | Z1, data, nsplit = 20, method = "PLS")
```

This version is experimental. In small simulations (n = 300, SO-PLS) the default (no
correction) had much better power than Nadeau-Bengio but rejected slightly too often
under the null (0.05-0.08); Nadeau-Bengio was valid but had low power. Like the other
tests it fails on the nonlinear scenario. Use the Monte Carlo test for a final
p-value until the fast version has been calibrated (details in `LOG.md`).

## Residual-based tests

Three tests share one idea: remove Z from X and from Y with a regression, then test
by permutation whether the residuals are still dependent. They differ only in the
regression, and all take `formula` and `data`:

```r
pls_ci_test(Y ~ X | Z1 + Z2, data)       # PLS
rosa_ci_test(Y ~ X | Z1 + Z2, data)      # ROSA (Z blocks kept separate)
multixg_ci_test(Y ~ X | Z1 + Z2, data)   # multi-target XGBoost, cross-fitted (needs more work)
```

`sim_scale()` simulates blocks of any size for large-N and high-dimensional studies.
In our simulations (`inst/benchmarks/scale_timing.R`, `scale_power.R`, results in `LOG.md`)
`pls_ci_test()` handles one million rows (about 5 to 7 minutes) and thousands of
variables (under a minute); `rosa_ci_test()` is about 2.5 times slower. Type I error is
close to nominal for Gaussian, Laplace, skewed and heavy-tailed noise, but when the
noise level depends on Z use `perm = "conditional"`, and in high dimension use
cross-validated residuals (`resid = "cv"`).

`multixg_ci_test()` predicts each fold from a model trained on the other folds
(`nfolds = 5`), so no row is predicted by a model that has seen it. Experimental: in a
small simulation (`inst/benchmarks/multixg_power.R`) it rejected true nulls far too
often when Z affects X and Y in a way that is hard to learn from a few hundred rows
(see `LOG.md`).

## Development

Open `BlockCI.Rproj` in RStudio.
`devtools::document()` rebuilds docs and NAMESPACE (roxygen2), `devtools::test()` runs tests,
`devtools::install()` installs. Requires the `multiblock` and `pls` packages.

- `R/Block.CCI.test.R`: the block CCI test (Monte Carlo)
- `R/cci_test_gen.R`, `R/cci_null_dist.R`: one statistic and the null distribution
- `R/learners.R`: the regression methods (SOPLS, PLS, XGB, RF) and how to add one
- `R/fast.R`: the fast resampled-t version (`Block.CCI.fast()`)
- `R/plot.R`: plot method for test results
- `R/scenarios.R`: data generators for null scenarios
- `inst/examples/`: demo and benchmark scripts
- `inst/benchmarks/`: slow simulations for the SO-PLS test (type I error and power);
  `simulation_study.R` is a scratch script for the current simulation runs
- `LOG.md`: project log
