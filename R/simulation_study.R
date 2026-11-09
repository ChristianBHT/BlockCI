library(multiblock)
library(parallel)
library(pbapply)

# small tests
data <- sim_1()
test <- SOPLS.CCI.test(Y ~ X | Z1 + Z2, data = data)
test
test <- SOPLS.CCI.test(Y ~ X | Z1, data = data)
test

#simulation 1

## ---------------- Design ----------------
n_grid   <- c(100, 250, 500)
nsim     <- 500          # replikasjoner per celle (MC-SE ved p=0.05 er ca. 0.01)
alpha    <- 0.05
effects  <- c(null = 0, power = 0.3)   # effekt-størrelse for H1; tilpass til sim_1

conditioning <- list(
  "Z1+Z2" = Y ~ X | Z1 + Z2,   # korrekt betingingsmengde
  "Z1"    = Y ~ X | Z1         # utelater Z2 (misspesifisert)
)

design <- expand.grid(
  n = n_grid,
  effect = names(effects),
  cond = names(conditioning),
  stringsAsFactors = FALSE
)

## ---------------- Én replikasjon ----------------
run_one <- function(i, n, effect, formula) {
  data <- sim_1(n = n, effect = effect)
  tryCatch(
    SOPLS.CCI.test(formula, data = data, progress = FALSE)$p.value,
    error = function(e) NA_real_
  )
}

## ---------------- Klynge ----------------
cl <- makeCluster(max(1, detectCores() - 1))
on.exit(stopCluster(cl), add = TRUE)

clusterEvalQ(cl, library(multiblock))
fns <- ls(globalenv(), all.names = TRUE)
fns <- fns[vapply(fns, function(x) is.function(get(x, globalenv())), NA)]
clusterExport(cl, varlist = fns)
clusterExport(cl, varlist = c("run_one", "effects", "conditioning"))
clusterSetRNGStream(cl, 123)    # reproduserbare, uavhengige strømmer

## ---------------- Kjør alle celler ----------------
results <- vector("list", nrow(design))
pvals   <- vector("list", nrow(design))

for (k in seq_len(nrow(design))) {
  d <- design[k, ]
  message(sprintf("Celle %d/%d: n=%d, %s, betinger på %s",
                  k, nrow(design), d$n, d$effect, d$cond))
  
  p <- unlist(pblapply(seq_len(nsim), function(i) {
    run_one(i, n = d$n, effect = effects[[d$effect]],
            formula = conditioning[[d$cond]])
  }, cl = cl))
  
  ok  <- !is.na(p)
  rej <- sum(p[ok] < alpha)
  ci  <- binom.test(rej, sum(ok))$conf.int   # eksakt 95 % KI
  
  pvals[[k]]   <- p
  results[[k]] <- data.frame(
    d,
    nsim_ok = sum(ok),
    rejection_rate = rej / sum(ok),
    ci_low = ci[1], ci_high = ci[2]
  )
}

res <- do.call(rbind, results)
res$metric <- ifelse(res$effect == "null", "type I error", "power")

saveRDS(list(summary = res, pvalues = pvals, design = design),
        file = "simulation_study.rds")
write.csv(res, "simulation_study_summary.csv", row.names = FALSE)

print(res[, c("metric", "cond", "n", "rejection_rate", "ci_low", "ci_high", "nsim_ok")],
      row.names = FALSE, digits = 3)
