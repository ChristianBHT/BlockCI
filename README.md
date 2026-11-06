# CausalBlocks: causal discovery between blocks of data

R package. A "block" is a matrix; two settings are covered:

This project aims at three things. First is to create a viable test for conditional independence between blocks of data (matrices). 
Secondly, based on CI testing, create a causal discovery method. 
Thirdly, create a R-package. 
## SO-PLS conditional independence test

`SOPLS.CCI.test(Y ~ X | V, data)` tests Y independent of X given V, modelled on the CCI package.
Statistic: T = RMSEP(Y ~ V + X) - RMSEP(Y ~ V) on held-out data, from one `multiblock::sopls` fit.
Null distribution: R Monte Carlo CV instances on rotated X (`method = "rotation"`, default) or
permuted X (`method = "permutation"`). Components are chosen by inner CV (`ncomp = "auto"`) or given as `c(a, b)`.
p-value = share of null statistics <= T_obs.

```r
res <- SOPLS.CCI.test(Y ~ X | V, list(Y = Y, X = X, V = V), nperm = 99)
res
```

## Development

Open `Causal Discovery between blocks.Rproj` in RStudio.
`devtools::document()` rebuilds docs and NAMESPACE (roxygen2), `devtools::test()` runs tests,
`devtools::install()` installs. Requires the `multiblock` package.

- `inst/examples/`: demo and benchmark scripts
- `inst/benchmarks/`: slow simulations for the SO-PLS test (type I error and power)
- `LOG.md`: project log
