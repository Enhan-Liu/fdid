# Reproducible plots for the continuous-G validation diagnostics.
# Run from the package root with:
#   Rscript inst/plot_validation_results.R

args <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", grep("^--file=", args, value = TRUE)[1])
if (is.na(file_arg) || !nzchar(file_arg)) {
  script_dir <- "inst"
} else {
  script_dir <- dirname(normalizePath(file_arg))
}
pkg_root <- normalizePath(file.path(script_dir, ".."))
out_dir <- file.path(pkg_root, "inst", "validation_figures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

read_result <- function(...) {
  read.csv(file.path(pkg_root, ...), stringsAsFactors = FALSE)
}

write_result <- function(x, name) {
  utils::write.csv(x, file.path(out_dir, name), row.names = FALSE)
}

short_scenario <- function(x) {
  out <- x
  out[out == "heavy_tail_interactions"] <- "Heavy tail"
  out[out == "heteroskedastic_interactions"] <- "Heterosked."
  out[out == "mixture_kink"] <- "Mixture/kink"
  out
}

short_density <- function(x) {
  out <- x
  out[out == "residual_kde"] <- "Residual KDE"
  out[out == "location_scale"] <- "Location-scale"
  out
}

plot_png <- function(filename, width = 1800, height = 1100, res = 180,
                     pointsize = 18) {
  grDevices::png(file.path(out_dir, filename), width = width, height = height,
                 res = res, pointsize = pointsize, type = "cairo")
}

cols <- c("#2c7fb8", "#f03b20", "#31a354", "#756bb1", "#636363")
density_cols <- c("Residual KDE" = "#2c7fb8", "Location-scale" = "#f03b20")

# Figure 1: oracle RMSE convergence in the smooth DGP.
oracle <- read_result("inst", "oracle_signal_validation_results", "oracle_summary.csv")
oracle_plot_data <- data.frame(
  n = rep(oracle$n, 3L),
  target = rep(c("Level curve", "Derivative curve", "Scalar derivative"),
               each = nrow(oracle)),
  rmse = c(oracle$level_rmse, oracle$derivative_rmse, oracle$scalar_rmse),
  coverage = c(oracle$level_coverage, oracle$derivative_coverage,
               oracle$scalar_coverage)
)
write_result(oracle_plot_data, "oracle_rmse_plot_data.csv")

plot_png("oracle_rmse.png")
mat <- cbind(oracle$level_rmse, oracle$derivative_rmse, oracle$scalar_rmse)
matplot(oracle$n, mat, type = "b", pch = c(16, 17, 15), lwd = 2,
        col = cols[1:3], xlab = "Sample size n", ylab = "RMSE",
        main = "Oracle-signal validation: RMSE falls with sample size",
        cex = 1.2, cex.lab = 1.15, cex.axis = 1.05, cex.main = 1.15)
legend("topright", legend = c("Level curve", "Derivative curve", "Scalar derivative"),
       col = cols[1:3], pch = c(16, 17, 15), lwd = 2, bty = "n", cex = 1.0)
grid()
grDevices::dev.off()

# Figure 2: stress-test diagnostics at n = 1000.
stress <- read_result("inst", "oracle_signal_stress_results",
                      "stress_estimated_summary.csv")
stress1000 <- stress[stress$n == 1000, ]
stress1000$scenario_label <- short_scenario(stress1000$scenario)
stress1000$density_label <- short_density(stress1000$density_method)
stress1000 <- stress1000[order(stress1000$scenario_label, stress1000$density_label), ]
stress_plot_data <- stress1000[, c(
  "scenario", "scenario_label", "density_method", "density_label",
  "flex_level_rmse", "flex_derivative_rmse",
  "incremental_scalar_rmse", "incremental_scalar_coverage"
)]
write_result(stress_plot_data, "stress_diagnostics_plot_data.csv")

scens <- c("Heavy tail", "Heterosked.", "Mixture/kink")
dens <- c("Residual KDE", "Location-scale")
make_mat <- function(value) {
  out <- matrix(NA_real_, nrow = length(dens), ncol = length(scens),
                dimnames = list(dens, scens))
  for (i in seq_len(nrow(stress1000))) {
    out[stress1000$density_label[i], stress1000$scenario_label[i]] <- stress1000[[value]][i]
  }
  out
}

plot_png("stress_diagnostics.png", width = 2100, height = 1300)
layout(matrix(c(1, 2, 3), nrow = 1))
op <- par(mar = c(7, 5, 3, 1), cex = 1.05)
barplot(make_mat("flex_level_rmse"), beside = TRUE, col = density_cols,
        las = 2, ylim = c(0, 0.13), ylab = "RMSE",
        main = "Level curve", cex.names = 0.8, cex.axis = 0.9,
        cex.lab = 1.0, cex.main = 1.05)
legend("topright", legend = dens, fill = density_cols, bty = "n", cex = 0.85)
barplot(make_mat("flex_derivative_rmse"), beside = TRUE, col = density_cols,
        las = 2, ylim = c(0, 0.35), ylab = "RMSE",
        main = "Derivative curve", cex.names = 0.8, cex.axis = 0.9,
        cex.lab = 1.0, cex.main = 1.05)
barplot(make_mat("incremental_scalar_coverage"), beside = TRUE,
        col = density_cols, las = 2, ylim = c(0, 1.05), ylab = "Coverage",
        main = "Scalar coverage", cex.names = 0.8, cex.axis = 0.9,
        cex.lab = 1.0, cex.main = 1.05)
abline(h = 0.95, lty = 2, lwd = 2, col = "#555555")
par(op)
grDevices::dev.off()

# Figure 3: nuisance-oracle ablation for the scalar derivative at n = 1000.
abl <- read_result("inst", "nuisance_oracle_ablation_results",
                   "nuisance_oracle_ablation_summary.csv")
abl1000 <- abl[abl$n == 1000, ]
sel <- rbind(
  transform(abl1000[abl1000$layer == "estimated_outcome_oracle_density", ],
            component = "Outcome ML only"),
  transform(abl1000[abl1000$layer == "oracle_outcome_estimated_density" &
                      abl1000$density_method == "residual_kde", ],
            component = "Density only: residual KDE"),
  transform(abl1000[abl1000$layer == "oracle_outcome_estimated_density" &
                      abl1000$density_method == "location_scale", ],
            component = "Density only: location-scale"),
  transform(abl1000[abl1000$layer == "estimated_all" &
                      abl1000$density_method == "residual_kde", ],
            component = "Estimated all: residual KDE"),
  transform(abl1000[abl1000$layer == "estimated_all" &
                      abl1000$density_method == "location_scale", ],
            component = "Estimated all: location-scale")
)
sel$scenario_label <- short_scenario(sel$scenario)
component_order <- c("Outcome ML only", "Density only: residual KDE",
                     "Density only: location-scale",
                     "Estimated all: residual KDE",
                     "Estimated all: location-scale")
sel <- sel[order(sel$scenario_label, match(sel$component, component_order)), ]
ablation_plot_data <- sel[, c("scenario", "scenario_label", "component",
                              "scalar_vs_oracle_rmse",
                              "derivative_vs_oracle_rmse",
                              "level_vs_oracle_rmse",
                              "scalar_coverage")]
write_result(ablation_plot_data, "ablation_scalar_plot_data.csv")

abl_mat <- matrix(NA_real_, nrow = length(component_order), ncol = length(scens),
                  dimnames = list(component_order, scens))
for (i in seq_len(nrow(sel))) {
  abl_mat[sel$component[i], sel$scenario_label[i]] <- sel$scalar_vs_oracle_rmse[i]
}

plot_png("ablation_scalar_rmse.png", width = 2100, height = 950)
op <- par(mar = c(7, 5, 3, 2), cex = 1.05)
barplot(abl_mat, beside = TRUE, col = cols, las = 2,
        ylim = c(0, max(abl_mat, na.rm = TRUE) * 1.25),
        ylab = "Scalar RMSE vs oracle-all",
        main = "Nuisance-oracle ablation: scalar derivative",
        cex.names = 0.85, cex.axis = 0.9, cex.lab = 1.0, cex.main = 1.05)
legend("topleft", legend = component_order, fill = cols, bty = "n", cex = 0.8)
par(op)
grDevices::dev.off()

message("Saved validation figures and plot data to: ", out_dir)
