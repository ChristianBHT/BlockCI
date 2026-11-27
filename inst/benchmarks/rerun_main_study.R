# Rerun of the main simulation study (Section 4.2 of paper/disposition.tex)
# ----------------------------------------------------------------------------------------
# The data generators, test methods and settings are those of simulation_study_noise.R
# (only its definitions are read; none of the study calls in that file are run). Every
# cell (design, noise, method, n, effect) is saved in its own file, with the seed
# n + 10 * effect. The data seeds do not depend on the method, so all methods and parts
# use the SAME datasets in a given design, noise, n and effect and can be compared pair
# by pair.
#
# RESUMABLE: a cell whose file is complete and free of errors is skipped. A job can be
# stopped at any time and started again; only the cell it was working on is lost.
# Nothing in the old result folders is touched.
#
# Run from the project folder, after installing the current package (devtools::install()).
#
#   Terminal:   Rscript inst/benchmarks/rerun_main_study.R                  # plan only
#               Rscript inst/benchmarks/rerun_main_study.R cci_linear       # run a part
#               Rscript inst/benchmarks/rerun_main_study.R all              # everything
#               Rscript inst/benchmarks/rerun_main_study.R resid_linear threads=1
#
#   R console or RStudio job script:
#               source("inst/benchmarks/rerun_main_study.R")
#               rerun_plan()
#               run_rerun("cci_linear")
#
# To run parts at the same time, start one R process per part and give each threads=1
# (the forests otherwise use all cores in every process). Do not start the same part
# twice. Keep the computer from sleeping during long runs.
#
# Parts, in order of priority:
#   cci_linear      CCI with PLS, ROSA, RF; linear design
#   cci_nonlinear   CCI with RF; nonlinear design
#   cci_sopls       CCI with SO-PLS; linear design, four noise types (as before)
#   resid_linear    residual tests (RV with PLS and RF, weighted GCM); linear design
#   resid_nonlinear residual tests; nonlinear design (RV with PLS: Gaussian noise only)
#   highdim         high-dimensional Z, linear and nonlinear
#   resampling      CCI with permutation and conditional permutation (type I error)
# The weighted GCM is not affected by the code changes since the old runs; it is cheap
# and is rerun so that all tests share the same datasets.

library(CausalBlocks)

rerun_root <- "simulation_results/rerun_main"

# ---- definitions from the main simulation script ---------------------------------------

if (!file.exists("inst/benchmarks/simulation_study_noise.R"))
  stop("run this script from the project folder (BlockCI)")
sim <- new.env()
for (e in parse("inst/benchmarks/simulation_study_noise.R", keep.source = FALSE))
  if (is.call(e) && identical(e[[1]], as.name("<-"))) eval(e, sim)

# ---- the cells -------------------------------------------------------------------------

rerun_cells <- function() {
  noise <- c("gaussian", "t3", "laplace", "skewed", "poisson", "outlier", "hetero", "correlated")
  grid <- function(part, method, design, noise, n = c(200, 500, 1000), effect = c(0, 0.5))
    expand.grid(part = part, method = method, design = design, noise = noise, n = n,
                effect = effect, stringsAsFactors = FALSE)
  cells <- rbind(
    grid("cci_linear", c("CCI_PLS", "CCI_ROSA", "CCI_RF"), "linear", noise),
    grid("cci_nonlinear", "CCI_RF", "nonlinear", noise),
    grid("cci_sopls", "CCI_SOPLS", "linear", noise[1:4]),
    grid("resid_linear", c("RES_PLS", "RES_RF", "RES_RF_wgcm"), "linear", noise),
    grid("resid_nonlinear", c("RES_RF", "RES_RF_wgcm"), "nonlinear", noise),
    grid("resid_nonlinear", "RES_PLS", "nonlinear", "gaussian"),
    grid("highdim", c("RES_PLS", "RES_PLS_pca", "CCI_PLS", "CCI_PLS_pca"), "highdim_linear",
         "gaussian", effect = c(0, 0.2, 0.5)),
    grid("highdim", c("RES_RF_pls", "RES_RF_pca", "RES_RF_pca_wgcm", "CCI_RF_pca", "CCI_RF_pls"),
         "highdim", "gaussian"),
    grid("resampling", c("CCI_PLS_perm", "CCI_PLS_cond", "CCI_RF_perm", "CCI_RF_cond"), "linear",
         c("gaussian", "hetero", "outlier", "skewed"), n = c(500, 1000), effect = 0),
    grid("resampling", c("CCI_PLS_perm", "CCI_PLS_cond", "CCI_RF_perm", "CCI_RF_cond"), "nonlinear",
         c("hetero", "outlier", "skewed"), n = c(500, 1000), effect = 0))
  # within a part: small n first and type I error before power, so that complete
  # tables for the smaller sample sizes are ready early
  cells <- cells[order(match(cells$part, unique(cells$part)), cells$n, cells$effect), ]
  cells$study <- ifelse(grepl("^highdim", cells$design), paste0("highdim_", cells$noise),
                        paste0(cells$noise, "_noise"))
  rownames(cells) <- NULL
  cells
}

cell_file <- function(cell, root) {
  file.path(root, cell$part, sprintf("%s_%s_%s_n%d_effect%s.rds", cell$study, cell$design,
                                     cell$method, as.integer(cell$n), format(cell$effect)))
}

# complete = readable, the requested number of tests, no failed test
cell_done <- function(file, n_tests) {
  if (!file.exists(file)) return(FALSE)
  x <- tryCatch(readRDS(file), error = function(e) NULL)
  !is.null(x) && nrow(x$results) == n_tests && !anyNA(x$results$p)
}

.parts <- function(parts) {
  known <- unique(rerun_cells()$part)
  if (identical(parts, "all")) return(known)
  bad <- setdiff(parts, known)
  if (length(bad)) stop("unknown part: ", paste(bad, collapse = ", "),
                        "; available: all, ", paste(known, collapse = ", "))
  parts
}

# ---- plan: what is left, with a rough time estimate --------------------------------------

# Median seconds per test in the OLD result files, by design, method and n (older code,
# training fraction 0.75, varying computer load: a rough guide only).
old_seconds <- function() {
  files <- list.files("simulation_results", "[.]rds$", recursive = TRUE, full.names = TRUE)
  files <- files[!grepl("/(rerun[^/]*|followup[^/]*|test|sim_noise_results|sim_summaries)/", files) &
                 grepl("_results[^/]*/", files)]
  s <- do.call(rbind, lapply(files, function(f)
    tryCatch(readRDS(f)$summary[, c("design", "method", "n", "median_sec")],
             error = function(e) NULL)))
  if (is.null(s)) return(NULL)
  stats::aggregate(median_sec ~ design + method + n, s, stats::median)
}

rerun_plan <- function(n_tests = 100, root = rerun_root) {
  cells <- rerun_cells()
  cells$done <- vapply(seq_len(nrow(cells)),
                       function(i) cell_done(cell_file(cells[i, ], root), n_tests), TRUE)
  old <- old_seconds()
  key <- function(d) paste(d$design, d$method, d$n)
  sec <- if (is.null(old)) rep(NA_real_, nrow(cells)) else old$median_sec[match(key(cells), key(old))]
  cells$hours <- sec * n_tests / 3600
  out <- do.call(rbind, lapply(split(cells, factor(cells$part, unique(cells$part))), function(s)
    data.frame(part = s$part[1], cells = nrow(s), done = sum(s$done),
               hours_left = round(sum(s$hours[!s$done]), 1))))
  rownames(out) <- NULL
  print(out, row.names = FALSE)
  cat(sprintf("total: %d cells, %d done, about %.0f hours left in one process\n",
              nrow(cells), sum(cells$done), sum(cells$hours[!cells$done])))
  cat("(hours: rough, from the median times of the old runs; NA if there are none)\n")
  invisible(out)
}

# ---- running -----------------------------------------------------------------------------

# Stop if the installed package differs from the code in R/: the tests are run with the
# installed package.
check_installed <- function() {
  src <- new.env()
  for (f in list.files("R", "[.][Rr]$", full.names = TRUE)) sys.source(f, src, keep.source = FALSE)
  inst <- asNamespace("CausalBlocks")
  code <- function(f) paste(c(deparse(formals(f), control = NULL),
                              deparse(body(f), control = NULL)), collapse = "\n")
  differs <- Filter(function(name) {
    f <- get(name, src)
    is.function(f) && (!exists(name, inst, inherits = FALSE) ||
                         !identical(code(f), code(get(name, inst))))
  }, ls(src, all.names = TRUE))
  if (length(differs))
    stop("the installed CausalBlocks differs from the code in R/ (", paste(utils::head(differs, 5), collapse = ", "),
         if (length(differs) > 5) ", ...", "). Restart R and run devtools::install().",
         call. = FALSE)
}

# the Git commit of the code ("+changes": R/ or the simulation script differ from it)
code_version <- function() {
  git <- function(...) tryCatch(suppressWarnings(system2("git", c(...), stdout = TRUE, stderr = FALSE)),
                                error = function(e) character(0))
  head <- git("rev-parse", "--short", "HEAD")
  if (!length(head)) return(NA_character_)
  paste0(head[1], if (length(git("status", "--porcelain", "--", "R",
                                 "inst/benchmarks/simulation_study_noise.R"))) "+changes")
}

# parts: names of parts, or "all". n_tests: datasets per cell. threads: threads per
# forest (NULL = all cores). root: result folder.
run_rerun <- function(parts, n_tests = 100, threads = NULL, root = rerun_root) {
  parts <- .parts(parts)
  check_installed()
  cells <- rerun_cells()
  cells <- cells[cells$part %in% parts, ]
  cells <- cells[order(match(cells$part, parts)), ]           # parts in the order asked for
  version <- code_version()
  cat(sprintf("%d cells in %s; code %s; results in %s\n", nrow(cells),
              paste(parts, collapse = ", "), version, root))
  failed <- character(0)
  start <- Sys.time()
  for (i in seq_len(nrow(cells))) {
    cell <- cells[i, ]
    file <- cell_file(cell, root)
    if (cell_done(file, n_tests)) next
    cat(sprintf("\n[%s] cell %d of %d (%s)\n", format(Sys.time(), "%d %b %H:%M"), i,
                nrow(cells), cell$part))
    args <- list(cell$method, effect = cell$effect, n = cell$n, n_tests = n_tests,
                 seed = cell$n + 10 * cell$effect,
                 settings = if (is.null(threads)) list() else list(threads = threads))
    out <- if (grepl("^highdim", cell$design)) {
      do.call(sim$study_highdim, c(args, list(noise = cell$noise, nonlinear = cell$design == "highdim")))
    } else {
      do.call(sim[[paste0("study_", cell$noise, "_noise")]], c(args, list(design = cell$design)))
    }
    out$settings$code_version <- version
    out$settings$finished <- Sys.time()
    # write to a temporary file first: an interrupted save must not look like a result
    dir.create(dirname(file), showWarnings = FALSE, recursive = TRUE)
    tmp <- paste0(file, ".tmp")
    saveRDS(out, tmp)
    if (file.exists(file)) file.remove(file)
    if (!file.rename(tmp, file)) stop("could not save ", file)
    if (out$summary$errors > 0) {
      failed <- c(failed, basename(file))
      # every test failed: something is wrong with the setup, not with a dataset
      if (out$summary$errors == n_tests) stop("all tests failed in ", basename(file), ": ",
                                              out$results$error[1])
    }
  }
  cat(sprintf("\nfinished %s after %.1f hours\n", paste(parts, collapse = ", "),
              as.numeric(difftime(Sys.time(), start, units = "hours"))))
  if (length(failed))
    cat("cells with failed tests (they are run again next time):\n", paste(" ", failed, collapse = "\n"), "\n")
  invisible(failed)
}

# One table with the summary rows of all finished cells.
rerun_summary <- function(root = rerun_root) {
  files <- list.files(root, "[.]rds$", recursive = TRUE, full.names = TRUE)
  if (!length(files)) stop("no results in ", root)
  out <- do.call(rbind, lapply(files, function(f) {
    x <- readRDS(f)
    cbind(part = basename(dirname(f)), x$summary, code_version = x$settings$code_version)
  }))
  out <- out[order(out$part, out$study, out$design, out$method, out$effect, out$n), ]
  rownames(out) <- NULL
  out
}

# ---- command line --------------------------------------------------------------------

# runs only when THIS file is the script given to Rscript, not when it is source()d
if (any(grepl("^--file=.*rerun_main_study[.]R$", commandArgs()))) local({
  args <- commandArgs(TRUE)
  opt <- grepl("=", args)
  value <- function(name, default) {
    a <- grep(paste0("^", name, "="), args[opt], value = TRUE)
    if (length(a)) sub("^[^=]*=", "", a[1]) else default
  }
  unknown <- setdiff(sub("=.*", "", args[opt]), c("threads", "n_tests", "root"))
  if (length(unknown)) stop("unknown option: ", paste(unknown, collapse = ", "))
  n_tests <- as.integer(value("n_tests", 100))
  root <- value("root", rerun_root)
  threads <- value("threads", NULL)
  if (!any(!opt)) {
    rerun_plan(n_tests, root)
    cat("Give one or more parts, or 'all', to run.\n")
  } else {
    run_rerun(args[!opt], n_tests = n_tests, root = root,
              threads = if (is.null(threads)) NULL else as.integer(threads))
  }
})
