# Project log: Causal Discovery between blocks of data

## 2026-11-05
- Project created. Language: base R.
- Blocks can mean (A) same samples, different variables, or (B) same variables, different samples.
- Added simulators (`R/simulate.R`), baseline methods (`R/discover.R`), evaluation (`R/evaluate.R`), `run_demo.R`, `benchmark.R`, `README.md`.
- Results: setting A over 30 seeds (n=1000, alpha=0.05): precision 0.84, recall 0.74, SHD 1.13. Setting B (ICP) on demo seed: precision 1.0, recall 0.4.
- Project moved to `Documents\Causal Discovery between blocks`.

- Read thesis section 3.1.4 and Article 1 summary (CCI test, pp. 36 and 39-41).
- Added SO-PLS version of CCI: `R/socci.R` (`SOPLS.CCI.test`, `sopls_perm_test`, `sopls_test_gen`), same structure as the CCI package. Uses `multiblock::sopls` (installed from CRAN); one fit gives M0 (V only, "a,0") and M1 (V then X, "a,b"); T = M1 - M0 (RMSEP); p = share of T* <= T_obs. Resampling of X: rotation (default, linear V-adjustment, Langsrud-style) or permutation.
- Checks in `test_socci.R` (40 sims, R = 49, n = 200, linear sim): type I error 0.03-0.07 for both resampling methods. Power depends strongly on `ncomp`: effect 0.6 gives 0.10-0.17 with ncomp=c(3,3) but 0.65-0.80 with ncomp=c(5,2) (V has 5 columns, so 3 comps left V information unexplained, which X then picks up).

- Added `ncomp = "auto"` (now the default): inside every MCCV split, (a, b) chosen by 5-fold CV on the training part only (a = best V-only RMSECV in 1..10; b = best added X comps in 0..5, b = 0 gives T = 0). `max_ncomp` and `cv_segments` control it. Also fixed `rotate_X` for the no-V case.
- Results with auto (30 sims, R = 39, n = 200): H0 rejection 0.03 (rotation) / 0.00 (permutation); effect 0.6: 0.57 / 0.63; effect 1.5: 0.83 / 0.83. Unconditional (Y ~ X) check OK (`test_unconditional.R`). Cost: 30 sims x 6 settings took ~5 min, so each test is slow (CV in every split).

- Converted to an R package named `CausalBlocks` (rename in DESCRIPTION if wanted): DESCRIPTION, roxygen2 docs and NAMESPACE, `.Rbuildignore`, `.gitignore`, testthat tests (15 expectations pass). Scripts moved to `inst/examples/` and `inst/benchmarks/`. Installed roxygen2 and testthat from CRAN. Note: `R/simulate.R`, `R/discover.R`, `R/evaluate.R` were missing from the project folder after the move and were recreated from the earlier versions (now with `stats::`/`utils::` prefixes and roxygen docs).

- Assessed publishability (honest view: plausible for chemometrics/computational-statistics journals only with large type I/power study, nonlinear/non-Gaussian robustness of rotation test, competitor comparison, downstream PC results, real data).
- Added `sim_h0_scenario()` (`R/scenarios.R`): 4 scenarios, blocks Y(4) X(4) Z1(5) Z2(5), Y _||_ X | Z1,Z2 true by construction (X, Y generated from the OBSERVED Z blocks; first draft used latent factors and was wrong): (1) common causes linear Gaussian, (2) chain X<-Z1->Z2->Y, (3) X almost a linear function of Z (noise sd 0.1, collinear), (4) nonlinear + Laplace + heteroskedastic.
- Type I error study `inst/benchmarks/h0_scenarios.R`: n = 200/500/1000/2000, both resampling methods, ncomp auto, 30 reps, nperm 49 (timing: one test takes ~6 s at n = 200 and ~125 s at n = 2000 per core). Results to `inst/benchmarks/h0_scenarios_results.csv`.
- RESULTS (done 2026-11-06, rejection rate at 0.05, 30 reps per cell, same for rotation and permutation): scenarios 1-3 reject 0-7% (mean p 0.83-0.94, so conservative; mean chosen X comps 0.03-0.3, i.e. mostly b = 0 -> T = 0 in both observed and null -> ties counted as extreme). Scenario 4 (nonlinear): rejects 97-100% at every n, mean b about 2.1. Linear SO-PLS on V cannot capture the nonlinear Z -> Y mean, X carries that information, and rotation (linear V-adjustment) cannot repair it. Rotation and permutation were indistinguishable, including scenario 3 (strong X-Z dependence), so no advantage of rotation shown yet.
- RESULTS `inst/benchmarks/h0_conditioning_sets.R` (done 2026-11-06; 100 reps, n = 500, rotation, rejection rate at 0.05): given Z1,Z2 (H0 true): scen 1: 0.01, scen 2: 0.00, scen 3: 0.02, scen 4: 0.99 (fails). Given Z1 only: scen 2 (H0 true): 0.00; scen 1/3/4 (H0 false, power): 1.00 / 0.99 / 1.00. Both conditioning sets behave as expected in the linear scenarios; the nonlinear scenario is rejected under both, so its Z1-only power is not informative.

## 2026-11-06 (later)
- Repo now lives in `Documents\GitHub\BlockCI` (GitHub Desktop); first commit message written. R files are edited there, not in the old project folder.
- Speed: the test is slow because each of the 161 MCCV instances runs a 5-fold CV `sopls` over all (a, b) pairs plus a final fit (~1000 fits at the old defaults). Options noted: fixed `ncomp`, smaller `max_ncomp`/`cv_segments`, parallel null draws (not implemented), fewer `nperm`. Defaults of `sopls_test_gen` were lowered by the author (`max_ncomp = c(6, 10)`, `cv_segments = 3`).
- Multiple V blocks: `Y ~ X | Z1 + Z2` now keeps Z1 and Z2 as separate `sopls` blocks (before: cbind-ed). `ncomp = "auto"` currently gives each V block one component and selects only b (X components) by inner CV (the author's current setting in `R/socci.R`); explicit `ncomp = c(a_1, ..., a_k, b)`. Rotation still adjusts X for all V columns together. A per-block CV search for the V components was tried (n = 300, H0: chose 8 V comps and b = 0, p = 0.93; real effect: 7 and 4, p = 0.025) but costs `max_ncomp[1]^k x (max_ncomp[2] + 1)` combinations (15 s vs 3.6 s for 40 fits). One component per block alone looked low-powered in one run (4 variables per block leaves V signal unexplained; effect not detected, p = 0.275). Not yet checked in a proper simulation.
- Fixed leftovers from the `tr`/`te` -> `train`/`test` rename in `sopls_test_gen`.
- `pvalue = c("empirical", "parametric")` in `SOPLS.CCI.test`: parametric = `pnorm(T_obs, mean(T*), sd(T*))` (normal null). Both stored as `p.empirical` / `p.parametric`. Example (R = 39): 0.025 empirical vs 0.011 parametric.
- `plot()` method for test results (`R/plot.R`): null histogram + observed T line (no shaded tail); fitted normal curve for parametric p. Added `graphics` to Imports.
- Fast version `SOPLS.CCI.fast()` (`R/fast.R`): fit on real X and rotated X* with the same split, paired per-row squared-error comparison, repeated over `nsplit` splits. Added `split`/`losses` arguments to `sopls_test_gen`. About 1 s per test (vs ~15 s).
  - Single split with a normal approximation over test rows: badly anti-conservative (H0 rejection 0.26-0.33 at 0.05, n = 300 and 1000, nrot 1 or 5), because it ignores training-fit/split variability.
  - K splits with Nadeau-Bengio corrected t (default): valid but low power (n = 300, K = 8, 60 sims: reject 0.00 / 0.03 / 0.10 at effect 0 / 0.15 / 0.5). `correction = "none"` (plain t across splits): 0.08 / 0.23 / 0.58, slightly liberal. Monte Carlo test remains the safe option for a final p-value.
  - Only n = 300 with V = Z1, Z2 and 60 sims checked; n = 1000, nrot > 1 and the H0 scenarios not run.
- Package repair (README task): `R/scenarios.R` had undefined `f1`/`f2`/`p` in `sim_2`/`sim_3` and a broken local `lap` in `sim_4`; fixed, with `sim_simple`, `sim_1`-`sim_4` and `sim_h0_scenario` (wrapper used by the benchmark scripts) now documented and exported. `R/simulate.R` had been overwritten by a scratch simulation script, which broke package loading and removed `sim_multiview`/`sim_environments`; the scratch script is kept unchanged as `inst/benchmarks/simulation_study.R` and the two simulators were restored from the earlier version. Regenerated `NAMESPACE` and `man/` with roxygen2; installs into a temporary library and the 15 testthat expectations pass.

## Next
- Fast test: calibrate the variance between "none" and Nadeau-Bengio (or use more rotations/splits) and check it on the four H0 scenarios and n = 1000/2000.
- Decide the default for V components (one per block vs CV), then run a type I / power study for it.
- Nonlinear V: add polynomial/interaction expansion of V (like CCI's `poly`, `degree`) to both the SO-PLS fit and the rotation adjustment; rerun scenario 4.
- Break ties in the p-value (randomised) so the test is not conservative when b = 0.
- Use the SO-PLS CI test inside the block-level PC search (replace Fisher-z in `pcor_test`).
- (to be decided)
