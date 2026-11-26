# 2. Simulation results: residual (RES) and CCI tests ------------------------------------
library(ggplot2)

base <- "simulation_results"     # run from the project folder

# ---- read the summaries of all saved studies ----
dirs  <- list.dirs(base, recursive = FALSE)
dirs  <- dirs[grepl("_(residual|cci)_results$", dirs)]
files <- unlist(lapply(dirs, list.files, pattern = "\\.rds$", full.names = TRUE))
results <- do.call(rbind, lapply(files, function(f) readRDS(f)$summary))

# ---- tidy labels ----
method_labels <- c(CCI_PLS = "PLS", CCI_SOPLS = "SO-PLS", CCI_ROSA = "ROSA", CCI_RF = "RF",
                   RES_PLS = "PLS (RV)", RES_RF = "PLS + RF (RV)", RES_RF_wgcm = "PLS + RF (wGCM)")
noise_levels  <- c("gaussian", "t3", "laplace", "skewed", "poisson", "outlier", "hetero", "correlated")

results$noise  <- sub("_noise$", "", results$study)
results$noise  <- factor(results$noise, levels = intersect(noise_levels, unique(results$noise)))
results$family <- ifelse(grepl("^CCI", results$method), "CCI test", "Residual test")
results$method_label <- factor(method_labels[results$method], levels = unique(method_labels))
# type I error clearly too high: more than 2 Monte Carlo s.e. above the nominal 5 %
results$too_high <- results$effect == 0 &
  results$rejection_rate > 0.05 + 2 * sqrt(0.05 * 0.95 / results$tests)
results <- results[order(results$noise, results$n, results$method_label), ]

# ---- tables: one row per noise type and n, one column per test ----
make_table <- function(d, eff) {
  d <- d[d$effect == eff, ]
  d$cell <- sprintf("%.2f%s", d$rejection_rate, ifelse(d$too_high, "*", ""))
  w <- reshape(d[, c("noise", "n", "method", "cell")], idvar = c("noise", "n"),
               timevar = "method", direction = "wide")
  names(w) <- sub("^cell\\.", "", names(w))
  cols <- intersect(names(method_labels), names(w))
  w <- w[order(w$noise, w$n), c("noise", "n", cols)]
  w[is.na(w)] <- "-"                                   # cells not run (yet)
  names(w) <- c("Noise", "n", paste0(ifelse(grepl("^CCI", cols), "CCI: ", "RES: "), method_labels[cols]))
  rownames(w) <- NULL
  w
}
type1_table <- make_table(results, 0)
power_table <- make_table(results, 0.5)

print(knitr::kable(type1_table, align = c("l", rep("r", ncol(type1_table) - 1)),
                   caption = "Type I error (H0 true), nominal 0.05. * = more than 2 MC s.e. above 0.05"))
print(knitr::kable(power_table, align = c("l", rep("r", ncol(power_table) - 1)),
                   caption = "Power (effect = 0.5)"))

dir.create(file.path(base, "sim_summaries"), showWarnings = FALSE)
write.csv(results,     file.path(base, "sim_summaries/all_results.csv"), row.names = FALSE)
write.csv(type1_table, file.path(base, "sim_summaries/type1_table.csv"), row.names = FALSE)
write.csv(power_table, file.path(base, "sim_summaries/power_table.csv"), row.names = FALSE)

# ---- heatmaps ----
row_lab <- paste0(results$noise, ", n = ", results$n)
results$row <- factor(row_lab, levels = rev(unique(row_lab)))      # first noise type at the top

heat_theme <- theme_minimal(base_size = 11) +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 30, hjust = 1),
        strip.text = element_text(face = "bold"),
        plot.title = element_text(face = "bold"))

# Type I error: diverging around the nominal 0.05 (blue = conservative, red = too many rejections)
d0 <- results[results$effect == 0, ]
p_type1 <- ggplot(d0, aes(method_label, row, fill = rejection_rate)) +
  geom_tile(colour = "#fcfcfb", linewidth = 1) +
  geom_text(aes(label = sprintf("%.2f%s", rejection_rate, ifelse(too_high, "*", "")),
                fontface = ifelse(too_high, "bold", "plain")),
            size = 3.2, colour = "#2b2b29") +
  scale_fill_gradient2(low = "#3987e5", mid = "#f0efec", high = "#e34948", midpoint = 0.05,
                       limits = c(0, 0.2), oob = scales::squish, name = "Rejection\nrate") +
  facet_grid(. ~ family, scales = "free_x", space = "free_x") +
  labs(title = "Type I error (H0 true), nominal level 0.05",
       subtitle = "Colour capped at 0.20; * and bold = more than 2 Monte Carlo s.e. above 0.05",
       x = NULL, y = NULL) +
  heat_theme

# Power: one hue, light (low) to dark (high)
d1 <- results[results$effect == 0.5, ]
p_power <- ggplot(d1, aes(method_label, row, fill = rejection_rate)) +
  geom_tile(colour = "#fcfcfb", linewidth = 1) +
  geom_text(aes(label = sprintf("%.2f", rejection_rate), colour = rejection_rate > 0.6), size = 3.2) +
  scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = "#2b2b29"), guide = "none") +
  scale_fill_gradient(low = "#cde2fb", high = "#0d366b", limits = c(0, 1), name = "Power") +
  facet_grid(. ~ family, scales = "free_x", space = "free_x") +
  labs(title = "Power (effect = 0.5)", x = NULL, y = NULL) +
  heat_theme

print(p_type1)
print(p_power)

h <- function(d) 1.5 + 0.28 * length(unique(d$row))           # height grows with the number of rows
ggsave(file.path(base, "sim_summaries/heatmap_type1.png"), p_type1, width = 9, height = h(d0), dpi = 200)
ggsave(file.path(base, "sim_summaries/heatmap_power.png"), p_power, width = 9, height = h(d1), dpi = 200)
