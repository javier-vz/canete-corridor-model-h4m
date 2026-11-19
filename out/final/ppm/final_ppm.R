# final_ppm.R  (JCAA-80 final recompute, H4m)
# Copy of multi_scale_ppm_v2.R: ONLY the input/output file names change (see diff below).
# Input column E holds E' (the final affinity field); the summary header lines still say "v2"
# because they are copied verbatim; the file names say H4m.

suppressPackageStartupMessages({
  library(spatstat)
  library(raster)
  library(dplyr)
})
set.seed(42)

cat("[load] Reading fields_H4m_for_ppm.csv...\n")
df <- read.csv("fields_H4m_for_ppm.csv")
cat("[load]", nrow(df), "cells,",
    sum(df$is_site), "sites\n")

CELL <- 250
CRS_UTM <- 32718

xmin <- min(df$X) - CELL/2
xmax <- max(df$X) + CELL/2
ymin <- min(df$Y) - CELL/2
ymax <- max(df$Y) + CELL/2

cat("[window] Building observation window from valid cells...\n")
hull_idx <- chull(df$X, df$Y)
hull_idx <- rev(hull_idx)   # anti-clockwise for spatstat
Omega_win <- owin(poly = list(x = df$X[hull_idx],
                                y = df$Y[hull_idx]))
cat("[window] Area:", round(area.owin(Omega_win) / 1e6, 1), "km^2\n")

cat("[ppp] Building point pattern of 84 sites...\n")
sites_df <- df[df$is_site == 1, ]
site_ppp <- ppp(x = sites_df$X, y = sites_df$Y,
                 window = Omega_win)
cat("[ppp]", npoints(site_ppp), "sites in window\n")

cat("[raster] Rasterising F, L, E to 250 m grid...\n")
build_im <- function(df, col_name) {
  r <- rasterFromXYZ(df[, c("X", "Y", col_name)],
                      res = c(CELL, CELL),
                      crs = CRS(paste0("+init=epsg:", CRS_UTM)),
                      digits = 0)
  m <- as.matrix(r)
  m <- m[nrow(m):1, ]
  im(mat = m, xrange = c(xmin, xmax), yrange = c(ymin, ymax))
}
F_im <- build_im(df, "F")
L_im <- build_im(df, "L")
E_im <- build_im(df, "E")

scales <- c(250, 500, 750, 1000, 1500, 2000, 2750, 3500, 5000)
cat("\n[scales]", length(scales), "scales:",
    paste(scales, collapse = " "), "m\n")

smooth_and_fit <- function(covar_im, radius_m, covar_name) {
  if (radius_m <= 250) {
    smoothed <- covar_im
  } else {
    smoothed <- blur(covar_im, sigma = radius_m / 2,
                       normalise = TRUE, bleed = FALSE)
  }
  fit <- ppm(site_ppp, trend = ~cov,
              covariates = list(cov = smoothed), rbord = 0)
  data.frame(covar = covar_name, radius_m = radius_m,
             AIC = AIC(fit), deviance = fit$maxlogpl * (-2),
             coef_slope = coef(fit)[2], stringsAsFactors = FALSE)
}

results <- data.frame()
covariates <- list(F = F_im, L = L_im, E = E_im)
for (cname in names(covariates)) {
  cat("\n[fit] Covariate:", cname, "\n")
  for (r in scales) {
    cat("  scale r =", r, "m ... ")
    result <- tryCatch({
      smooth_and_fit(covariates[[cname]], r, cname)
    }, error = function(e) {
      cat("ERROR:", conditionMessage(e), "\n")
      data.frame(covar = cname, radius_m = r, AIC = NA,
                 deviance = NA, coef_slope = NA)
    })
    results <- rbind(results, result)
    cat("AIC =", round(result$AIC, 2), "\n")
  }
}

cat("\n[null] Fitting null (constant intensity)...\n")
fit_null <- ppm(site_ppp, trend = ~1, rbord = 0)
null_AIC <- AIC(fit_null)
cat("[null] AIC =", round(null_AIC, 2), "\n")

cat("\n=====================================================\n")
cat(" SUMMARY: multi-scale PPM AIC (v2 model, river-fix applied)\n")
cat("=====================================================\n")
for (cname in names(covariates)) {
  sub <- results[results$covar == cname, ]
  best_idx <- which.min(sub$AIC)
  cat(sprintf("\n  %s:  best r = %d m, AIC = %.2f, dAIC = %.2f\n",
              cname, sub$radius_m[best_idx], sub$AIC[best_idx],
              null_AIC - sub$AIC[best_idx]))
}
best_row <- results[which.min(results$AIC), ]
cat(sprintf("\n  OVERALL BEST: %s @ r = %d m, AIC = %.2f, dAIC = %.2f\n",
            best_row$covar, best_row$radius_m, best_row$AIC,
            null_AIC - best_row$AIC))

results$null_AIC <- null_AIC
results$delta_AIC_vs_null <- null_AIC - results$AIC
write.csv(results, "ppm_H4m_multi_scale_aic.csv", row.names = FALSE)

sink("ppm_H4m_summary.txt")
cat("# Multi-scale PPM analysis (Lewis 2024) -- MODEL V2 with river-buffer fix\n")
cat("# Lower Canete Valley, 84 archaeological sites, 12,826-cell Omega\n")
cat("# Canonical CA v2: eps=0.15, kappa=0.03, seed=42, river_buffer_m=125\n\n")
cat("Null model (constant intensity): AIC =", round(null_AIC, 2), "\n\n")
for (cname in names(covariates)) {
  cat("## Covariate:", cname, "\n")
  sub <- results[results$covar == cname, ]
  print(sub[, c("radius_m", "AIC", "delta_AIC_vs_null")], row.names = FALSE)
  cat("\n")
}
cat("## Overall best model\n")
print(best_row, row.names = FALSE)
sink()

cat("\n[done] Results saved to:\n")
cat("       ppm_H4m_multi_scale_aic.csv\n")
cat("       ppm_H4m_summary.txt\n")
