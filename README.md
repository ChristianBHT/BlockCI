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

## SO-PLS conditional independence test

`SOPLS.CCI.test(Y ~ X | V, data)` tests whether Y is independent of X given V. It is an
SO-PLS version of the Computational Conditional Independence (CCI) test, with the
same structure as the CCI package.

- **Statistic:** `T = RMSEP(Y ~ V + X) - RMSEP(Y ~ V)` on held-out data, from one
  `multiblock::sopls` fit. Negative `T` means X improves prediction beyond V.
- **Null distribution:** R Monte Carlo CV instances on a resampled X: rotation
  (`method = "rotation"`, default) or permutation (`method = "permutation"`).
- **Components:** chosen by inner cross-validation (`ncomp = "auto"`) or given as
  `c(a_1, ..., a_k, b)`: components for each of the k V blocks, then b for X.
  Several V blocks are written `Y ~ X | Z1 + Z2` and are kept as separate blocks.
- **p-value:** share of null statistics at or below the observed one
  (`pvalue = "empirical"`), or from a normal fit to the null statistics
  (`pvalue = "parametric"`).

## Examples

The data generators are in `R/scenarios.R` (`?sim_simple`). The package needs the
`multiblock` package:

```r
devtools::install()           # from the package folder, or devtools::load_all()
library(CausalBlocks)
```

### A simple null case

`sim_simple()` makes four matrices Y, X, Z1 and Z2 where X and Y both depend on Z1
and Z2 and nothing else, so Y is independent of X given (Z1, Z2) but not given Z1 alone.

```r
set.seed(1)
data <- sim_simple(n = 500)     # list(Y, X, Z1, Z2), 4 variables per block

# Y _||_ X | Z1, Z2 is true: we expect a large p-value
res <- SOPLS.CCI.test(Y ~ X | Z1 + Z2, data, nperm = 99)
res
plot(res)                        # null distribution with the observed statistic

# Y _||_ X | Z1 is false (Z2 is missing): we expect a small p-value
SOPLS.CCI.test(Y ~ X | Z1, data, nperm = 99)
```

### Parametric p-value

With few null draws the empirical p-value cannot go below `1/(R + 1)`. A normal
approximation of the null distribution gives a finer p-value, at the price of
assuming normality:

```r
res_par <- SOPLS.CCI.test(Y ~ X | Z1 + Z2, data, nperm = 99, pvalue = "parametric")
c(empirical = res_par$p.empirical, parametric = res_par$p.parametric)
plot(res_par)                    # adds the fitted normal curve
```

### Scenario with mixed signals

`sim_1()` generates blocks Y (4 variables), X (4), Z1 (5) and Z2 (5) from latent
common causes. Y is independent of X given (Z1, Z2) by construction.

```r
d1 <- sim_1(n = 500)
SOPLS.CCI.test(Y ~ X | Z1 + Z2, d1, nperm = 99)   # null is true
SOPLS.CCI.test(Y ~ X | Z1,      d1, nperm = 99)   # null is false
```

`sim_2()` (chain `X <- Z1 -> Z2 -> Y`), `sim_3()` (nearly collinear X) and `sim_4()`
(nonlinear, heteroskedastic) work the same way. The nonlinear `sim_4()` is known to
fail for the linear SO-PLS test: it rejects a true null too often (see `LOG.md`).

## Fast version

`SOPLS.CCI.fast()` avoids the Monte Carlo loop. For each of `nsplit` train/test splits
it fits SO-PLS once on the real X and once on a rotated copy X*, compares the two on
the same held-out rows, and combines the splits with a resampled t-test. It takes
about 1-2 seconds where the Monte Carlo test takes about 10 seconds or more.

```r
SOPLS.CCI.fast(Y ~ X | Z1 + Z2, data)                       # Nadeau-Bengio corrected (default)
SOPLS.CCI.fast(Y ~ X | Z1 + Z2, data, correction = "none")  # more power, can be liberal
SOPLS.CCI.fast(Y ~ X | Z1, data, nsplit = 20)
```

This version is experimental. In a small simulation (n = 300) the default correction
was valid but had low power, and `correction = "none"` had much better power but
rejected slightly too often under the null. Use the Monte Carlo test for a final
p-value until the fast version has been calibrated (details in `LOG.md`).

## Development

Open `BlockCI.Rproj` in RStudio.
`devtools::document()` rebuilds docs and NAMESPACE (roxygen2), `devtools::test()` runs tests,
`devtools::install()` installs. Requires the `multiblock` package.

- `R/socci.R`: the SO-PLS CI test (Monte Carlo)
- `R/fast.R`: the fast resampled-t version
- `R/plot.R`: plot method for test results
- `R/scenarios.R`: data generators for null scenarios
- `inst/examples/`: demo and benchmark scripts
- `inst/benchmarks/`: slow simulations for the SO-PLS test (type I error and power);
  `simulation_study.R` is a scratch script for the current simulation runs
- `LOG.md`: project log
