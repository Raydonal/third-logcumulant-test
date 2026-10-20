# =============================================================================
# article_code.R
#
# R code for the article "Discriminating between the Weibull and log-normal
# distributions with a test based on the third log-cumulant".
#
# The blocks follow the order of the text:
#   1. Functions: T3 statistic, RML rule and Kappenman's ratio
#   2. Section 3: constants and Figure 1 (log-cumulant diagram)
#   3. Section 4: correction factor, exact variances, Table 1 (Monte Carlo
#      calibration) and minimum sample size
#   4. Section 5: simulation study of the PCS of the three methods (Table 2
#      and Figure 2)
#   5. Section 6: application to the air-conditioning, ball-bearing and
#      carbon-fibre data, and Figure 3 (the three data sets in the
#      log-cumulant diagram)
#   6. Outputs: .xlsx workbook and LaTeX rows of Tables 1 and 2
#
# Notation of the text: Weibull with shape beta and scale lambda, log-normal
# with parameters mu and sigma, gamma with shape rho, significance level alpha.
#
# Seeds: in all simulations (Table 1, minimum sample size and Table 2), Monte
# Carlo replicate r uses set.seed(SEED + r - 1), with SEED = 123. The first
# replicate uses seed 123, the second uses 124, and so on up to the last one.
# The sequence is the same for every n and every configuration, so any
# replicate of any cell of the tables can be reproduced on its own (see
# generate_samples). The generator is set to Mersenne-Twister with normal
# deviates by inversion, the R default since version 3.6.0.
#
# Requirements: R (tested with version 4.3.3) and, optionally, the openxlsx
# package (without it, the tables are written as CSV files). The whole script
# runs in a few minutes.
#
# Outputs, written to the working directory: article_results.xlsx,
# fig_lc_diagram.pdf, fig_pcs.pdf, fig_data_diagram.pdf and latex_tables.tex
# =============================================================================

SEED    <- 123                       # seed of the first Monte Carlo replicate
RNGkind(kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")
ALPHA   <- 0.05
Z_ALPHA <- qnorm(ALPHA)              # -1.6449: standard normal critical value

# Cut-off of Kappenman (1988): the log-normal is chosen when the ratio r
# exceeds the cut-off. The treatment of non-integer values of 0.05 n and
# 0.2 n is described in kappenman_ratio.
KAPPENMAN_CUTOFF <- 0.7477


# =============================================================================
# 1. FUNCTIONS
# The first three take a matrix Y of log-data, with one sample per row
# (n = ncol(Y)), and return one value per sample. The last two generate the
# Monte Carlo replicates, with one seed per replicate.
# =============================================================================

# T3 statistic of eq. (7): bias-corrected sample third log-cumulant,
# standardised by the null variance 6 kappa2^3 with kappa2 estimated.
stat_T3 <- function(Y) {
  n  <- ncol(Y)
  Yc <- Y - rowMeans(Y)                  # deviations from the mean of each row
  k2 <- rowMeans(Yc^2)                   # sample second log-cumulant
  k3 <- rowMeans(Yc^3)                   # sample third log-cumulant
  cn <- n^2 / ((n - 1) * (n - 2))        # bias-correction factor, eq. (4)
  cn * sqrt(n) * k3 / sqrt(6 * k2^3)
}

# Difference between the maximised log-likelihoods of the Weibull and the
# log-normal (RML rule: the Weibull is chosen when the difference is
# positive). The log-normal MLE has a closed form. For the Weibull, the shape
# solves the profile equation
#   h(k) = sum(x^k log x)/sum(x^k) - 1/k - mean(log x) = 0
# by Newton's method, for all rows at once, and the scale then follows in
# closed form.
loglik_diff <- function(Y, iter = 50) {
  n    <- ncol(Y)
  ybar <- rowMeans(Y)
  s2   <- rowMeans((Y - ybar)^2)
  llLN <- -rowSums(Y) - 0.5 * n * log(2 * pi * s2) - 0.5 * n
  ymax <- Y[cbind(seq_len(nrow(Y)), max.col(Y, ties.method = "first"))]
  Yd   <- Y - ymax                       # avoids overflow in exp(k * y)
  k    <- 1.2825 / sqrt(s2)              # starting value: sd of log X = 1.2825/k
  for (it in seq_len(iter)) {
    W  <- exp(k * Yd)
    sw <- rowSums(W)
    A  <- rowSums(W * Y) / sw            # mean of y weighted by x^k
    B  <- rowSums(W * Y * Y) / sw
    h  <- A - 1 / k - ybar               # profile equation
    hp <- B - A^2 + 1 / k^2              # its derivative, always positive
    kn <- k - h / hp
    k  <- ifelse(kn > 0, kn, k / 2)      # guards against negative steps
  }
  W   <- exp(k * Yd)
  llW <- n * log(k) - n * (k * ymax + log(rowMeans(W))) + (k - 1) * rowSums(Y) - n
  # largest residual of the profile equation: checks convergence in every row
  residual <- max(abs(rowSums(W * Y) / rowSums(W) - 1 / k - ybar))
  structure(llW - llLN, residual = residual)
}

# Kappenman's (1988) ratio, r = (A3 - A2) / (A2 - A1), computed from the
# ordered logarithms. A1 and A3 are the means of the smallest 5% and of the
# largest 5% of the values, and A2 is the mean of the values that remain
# after discarding 20% on each side. When 0.05 n and 0.2 n are not integers,
# we use ceiling(0.05 n) observations in each tail and discard floor(0.2 n)
# on each side. Under the log-normal, r is close to 1. Under the Weibull,
# whose left tail is longer on the log scale, r falls below 1.
kappenman_ratio <- function(Y) {
  n   <- ncol(Y)
  S   <- t(apply(Y, 1, sort))
  if (nrow(Y) == 1) S <- matrix(S, nrow = 1)
  k5  <- ceiling(0.05 * n)
  k20 <- floor(0.20 * n)
  A1  <- rowMeans(S[, 1:k5, drop = FALSE])
  A3  <- rowMeans(S[, (n - k5 + 1):n, drop = FALSE])
  A2  <- rowMeans(S[, (k20 + 1):(n - k20), drop = FALSE])
  (A3 - A2) / (A2 - A1)
}

# Sample generation. Replicate r uses set.seed(SEED + r - 1), whatever the
# value of n, and its sample consists of the first n observations generated
# from that seed. Returns replicates r_start, ..., r_start + M - 1, one per
# row, on the original scale.
generate_samples <- function(M, n, generator, r_start = 1) {
  X <- matrix(0, M, n)
  for (i in seq_len(M)) {
    set.seed(SEED + (r_start + i - 1) - 1)   # seed of replicate r = r_start + i - 1
    X[i, ] <- generator(n)
  }
  X
}

# T3 in replicates 1, ..., M for each sample size in ns, computed in blocks of
# replicates to save memory. Replicate r uses set.seed(SEED + r - 1) and
# generates max(ns) observations, and its sample of size n consists of the
# first n of them. Since the generator produces the observations in sequence,
# these are exactly the observations that set.seed(SEED + r - 1) followed by
# the generation of n observations would produce. Replicate r therefore uses
# the same seed for every n (the check after Table 1 confirms this), and the
# samples of different sizes from the same replicate are nested. By
# Proposition 4.2, it suffices to simulate the standard log-normal, LN(0, 1),
# under H0 and the Weibull(1, 1), that is, the exponential distribution, under
# the alternative. Returns an M x length(ns) matrix.
simulate_T3 <- function(ns, M, family, block = 10000) {
  generator <- if (family == "LN") function(n) rlnorm(n, meanlog = 0, sdlog = 1)
               else                function(n) rweibull(n, shape = 1, scale = 1)
  out <- matrix(0, M, length(ns), dimnames = list(NULL, ns))
  for (start in seq(1, M, by = block)) {
    m <- min(block, M - start + 1)
    Y <- log(generate_samples(m, max(ns), generator, r_start = start))
    for (k in seq_along(ns))
      out[start:(start + m - 1), k] <- stat_T3(Y[, seq_len(ns[k]), drop = FALSE])
  }
  out
}


# =============================================================================
# 2. SECTION 3: CONSTANTS AND FIGURE 1
# =============================================================================

zeta3        <- -psigamma(1, 2) / 2                       # zeta(3)
gamma1_W     <- psigamma(1, 2) / psigamma(1, 1)^1.5       # -2 zeta(3)/(pi^2/6)^(3/2)
gamma1_gamma <- function(rho) psigamma(rho, 2) / psigamma(rho, 1)^1.5
constants <- data.frame(
  quantity = c("gamma1 (Weibull)", "gamma1 (gamma, rho = 1/2)",
               "gamma1 (gamma, rho = 2)", "gamma1 (gamma, rho = 5)"),
  value = c(gamma1_W, gamma1_gamma(0.5), gamma1_gamma(2), gamma1_gamma(5)))
print(constants, digits = 5)

# Figure 1: Weibull on the curve kappa3 = gamma1_W kappa2^(3/2), log-normal on
# the half-line kappa3 = 0 and gamma on the curve (psi''(rho), psi'(rho))
col_W <- "#1F4E79"; col_LN <- "#8B1A1A"; col_G <- "#2E7D32"
pdf("fig_lc_diagram.pdf", width = 6.2, height = 5)
par(mar = c(4.3, 4.6, 1, 1))
k3W <- seq(-4, 0, length.out = 400)
k2W <- (k3W / gamma1_W)^(2 / 3)
rho <- exp(seq(log(0.82), log(400), length.out = 800))
plot(NA, xlim = c(-4.1, 0.6), ylim = c(0, 4), las = 1,
     xlab = expression(kappa[3] ~~ "(third log-cumulant)"),
     ylab = expression(kappa[2] ~~ "(second log-cumulant)"))
abline(h = 0, col = "grey80")
lines(k3W, k2W, col = col_W, lwd = 2)
segments(0, 0, 0, 4, col = col_LN, lwd = 2)
lines(psigamma(rho, 2), psigamma(rho, 1), col = col_G, lwd = 2, lty = 2)
points(psigamma(c(2, 3), 2), psigamma(c(2, 3), 1), pch = 19, col = col_G)
text(psigamma(2, 2), psigamma(2, 1), expression(rho == 2), pos = 4, col = col_G, cex = 0.8)
# the label of rho = 3 sits to the right of the log-normal half-line so that it does not cross it
text(0.02, psigamma(3, 1), expression(rho == 3), pos = 4, col = col_G, cex = 0.8)
points(psigamma(1, 2), psigamma(1, 1), pch = 1, cex = 1.6, lwd = 2)   # exponential
text(psigamma(1, 2), psigamma(1, 1), "Exp", pos = 3, cex = 0.8)
text(-2.6, 0.5, expression(kappa[3] / kappa[2]^{3/2} == -1.1395), col = col_W, cex = 0.9)
text(0.3, 3.5, expression(kappa[3] == 0), col = col_LN, cex = 0.9)
legend("topleft", c("Weibull", "Log-normal", "Gamma"), col = c(col_W, col_LN, col_G),
       lty = c(1, 1, 2), lwd = 2, bty = "n")
dev.off()


# =============================================================================
# 3. SECTION 4: BIAS CORRECTION, CALIBRATION AND SAMPLE SIZE
# =============================================================================

# Correction factor n^2/((n-1)(n-2)) and exact null variances of the
# standardised statistic without and with the correction. Under normality,
# the sample skewness has variance 6(n-2)/((n+1)(n+3)) (Kendall and Stuart,
# 1977).
correction_factor <- function(n) n^2 / ((n - 1) * (n - 2))
var_uncorrected   <- function(n) n * (n - 2) / ((n + 1) * (n + 3))
var_corrected     <- function(n) correction_factor(n)^2 * var_uncorrected(n)
correction <- data.frame(n = c(20, 25, 50, 100),
                         factor = correction_factor(c(20, 25, 50, 100)),
                         var_uncorrected = var_uncorrected(c(20, 25, 50, 100)),
                         var_corrected = var_corrected(c(20, 25, 50, 100)))
print(correction, digits = 4)

# Table 1: for each n, the 5% quantile of the null distribution of T3 (to be
# compared with -1.6449), and the size and power of the rule "reject H0 if
# T3 <= -1.6449", that is, of the test as it is used in practice. The power
# is the PCS under the Weibull. We also keep the size of the statistic
# without the bias correction (T3 divided by the factor) and the null
# tabulations for n = 23 and n = 30, used for the Monte Carlo p-values in
# Section 6 (the one for n = 69 comes from the minimum-sample-size scan
# below). No seed is set here, because each replicate sets its own inside
# simulate_T3 (seeds 123 to 100122, the same for every n and for both
# families).
ns_tab1 <- c(20, 23, 25, 30, 40, 50, 75, 100, 150, 200, 300, 500)
M1      <- 1e5
T0 <- simulate_T3(ns_tab1, M1, "LN")      # under H0 (log-normal), one column per n
T1 <- simulate_T3(ns_tab1, M1, "W")       # under the Weibull
tab1 <- NULL
for (k in seq_along(ns_tab1)) {
  n  <- ns_tab1[k]
  t0 <- T0[, k]; t1 <- T1[, k]
  q05 <- quantile(t0, ALPHA, names = FALSE)
  tab1 <- rbind(tab1, data.frame(
    n = n, q005 = q05,
    size = mean(t0 <= Z_ALPHA),
    power = mean(t1 <= Z_ALPHA),
    size_uncorrected = mean(t0 / correction_factor(n) <= Z_ALPHA)))
}
null_tabs <- list("23" = T0[, "23"], "30" = T0[, "30"])
rm(T0, T1)
print(tab1, digits = 4)

# Check of the seed scheme: replicates of n = 23 regenerated on their own,
# from their seed alone, reproduce the tabulated values.
for (r in c(1, 2, M1)) {
  set.seed(SEED + r - 1)
  y <- log(rlnorm(23, meanlog = 0, sdlog = 1))
  stopifnot(isTRUE(all.equal(stat_T3(matrix(y, nrow = 1)), null_tabs[["23"]][r])))
}
cat("Seeds checked: replicates 1, 2 and", M1, "of n = 23 were regenerated on their own\n")

# Minimum sample size for power 0.90, 0.95 and 0.99 at the 5% level. For each
# n from 60 to 135, the size and power of the rule T3 <= -1.6449 are
# estimated with the same 1e5 replicates as in Table 1 (seeds 123 to 100122),
# and the smallest n whose power reaches the target is taken.
M2      <- M1
ns_scan <- 60:135
P0 <- simulate_T3(ns_scan, M2, "LN")
P1 <- simulate_T3(ns_scan, M2, "W")
power_scan <- data.frame(n = ns_scan, size = colMeans(P0 <= Z_ALPHA),
                         power = colMeans(P1 <= Z_ALPHA))
# null tabulation for n = 69 (same seeds), used for the Monte Carlo p-value of
# the carbon-fibre data in Section 6
null_tabs[["69"]] <- P0[, "69"]
rm(P0, P1)
min_n <- do.call(rbind, lapply(c(0.90, 0.95, 0.99), function(a) {
  n_min <- min(power_scan$n[power_scan$power >= a])
  # the range of n must contain the point at which the power crosses the target
  stopifnot(is.finite(n_min), n_min > min(ns_scan))
  data.frame(target_power = a, n_min = n_min)
}))
print(min_n)


# =============================================================================
# 4. SECTION 5: SIMULATION STUDY OF THE PCS (TABLE 2 AND FIGURE 2)
# The data are generated on the original scale, with R = 5000 samples per
# configuration, and the three methods are applied to the same samples:
#   T3        chooses the Weibull when T3 <= -1.6449
#   RML       chooses the Weibull when l_W - l_LN > 0
#   Kappenman chooses the Weibull when r <= 0.7477
# The PCS is the proportion of choices of the family that generated the data.
# All configurations use the same seeds per replicate (123 to 5122). Since
# the samples with different parameters are, replicate by replicate,
# transformations x -> b x^c of one another, and the three rules are
# invariant under these transformations (Proposition 4.2), the rows with the
# same n and different beta (or sigma) coincide. The largest difference
# between them is printed at the end.
# =============================================================================
ns_tab2 <- c(20, 30, 40, 50, 75, 100, 200, 300, 500)
R2      <- 5000
configs <- rbind(data.frame(family = "Weibull",    parameter = c(0.5, 1, 2, 4)),
                 data.frame(family = "Log-normal", parameter = c(0.5, 1)))
tab2 <- NULL
for (i in seq_len(nrow(configs))) {
  for (n in ns_tab2) {
    fam <- configs$family[i]; par_ <- configs$parameter[i]
    X <- if (fam == "Weibull") {
      generate_samples(R2, n, function(n) rweibull(n, shape = par_, scale = 1))  # beta = par_, lambda = 1
    } else {
      generate_samples(R2, n, function(n) rlnorm(n, meanlog = 0, sdlog = par_)) # mu = 0, sigma = par_
    }
    Y <- log(X)
    d <- loglik_diff(Y)
    chooses_W <- cbind(T3        = stat_T3(Y) <= Z_ALPHA,
                       RML       = as.vector(d) > 0,
                       Kappenman = kappenman_ratio(Y) <= KAPPENMAN_CUTOFF)
    correct <- if (fam == "Weibull") chooses_W else !chooses_W
    tab2 <- rbind(tab2, data.frame(family = fam, parameter = par_, n = n,
                                   PCS_T3 = mean(correct[, "T3"]),
                                   PCS_RML = mean(correct[, "RML"]),
                                   PCS_Kappenman = mean(correct[, "Kappenman"]),
                                   newton_residual = attr(d, "residual")))
  }
  cat(sprintf("Table 2: %s with parameter %.1f done\n", fam, par_))
}
print(tab2, digits = 3)
cat("Largest Newton residual in Table 2:", max(tab2$newton_residual), "\n")

# Checks: (i) identical PCS across the parameters of the same family for the
# same n; (ii) under the log-normal with sigma = 1 and n = 30, the PCS of T3
# equals that obtained from the first 5000 replicates of the null tabulation
# of Table 1 (same seeds).
cfg_diff <- aggregate(cbind(PCS_T3, PCS_RML, PCS_Kappenman) ~ family + n, data = tab2,
                      FUN = function(v) max(v) - min(v))
cat("Largest PCS difference across parameters (same family and same n):",
    max(cfg_diff[, c("PCS_T3", "PCS_RML", "PCS_Kappenman")]), "\n")
stopifnot(isTRUE(all.equal(
  tab2$PCS_T3[tab2$family == "Log-normal" & tab2$parameter == 1 & tab2$n == 30],
  mean(null_tabs[["30"]][1:R2] > Z_ALPHA))))
cat("PCS of T3 (log-normal, n = 30) equals that of the first 5000 replicates of Table 1\n")

# Figure 2: PCS of the three methods as a function of n, with beta = 1 for the
# Weibull data and sigma = 1 for the log-normal data. With the same seeds, the
# other parameter values give identical curves (Proposition 4.2).
cols       <- c(T3 = "#8B1A1A", RML = "#1F4E79", Kappenman = "#2E7D32")
pchs       <- c(T3 = 19, RML = 17, Kappenman = 15)
leg_labels <- c(expression(tilde(T)[3]), "RML", "Kappenman")
pdf("fig_pcs.pdf", width = 9.2, height = 4.3)
par(mfrow = c(1, 2), mar = c(4.3, 4.6, 2.2, 1))
for (fam in c("Weibull", "Log-normal")) {
  sub <- tab2[tab2$family == fam, ]
  plot(NA, xlim = range(ns_tab2), ylim = if (fam == "Weibull") c(0.4, 1) else c(0.7, 1),
       log = "x", xaxt = "n", las = 1, main = paste(fam, "data"),
       xlab = "n (sample size)", ylab = "PCS (prob. of correct selection)")
  axis(1, at = c(20, 50, 100, 200, 300, 500), gap.axis = 0.25)   # gap.axis: keeps the 300 label
  if (fam == "Log-normal") abline(h = 1 - ALPHA, lty = 2)
  s <- sub[sub$parameter == 1, ]                 # beta = 1 or sigma = 1
  for (met in names(cols)) {
    lines(s$n, s[[paste0("PCS_", met)]], col = cols[met])
    points(s$n, s[[paste0("PCS_", met)]], col = cols[met], pch = pchs[met], cex = 0.6)
  }
  if (fam == "Weibull") {
    legend("bottomright", legend = leg_labels, col = cols, pch = pchs, lty = 1, bty = "n", cex = 0.85)
  } else {
    legend("bottomright", legend = c(leg_labels, expression(1 - alpha == 0.95)),
           col = c(cols, "black"), pch = c(pchs, NA), lty = c(1, 1, 1, 2), bty = "n", cex = 0.85)
  }
}
dev.off()


# =============================================================================
# 5. SECTION 6: APPLICATIONS
# =============================================================================
aircon   <- c(23, 261, 87, 7, 120, 14, 62, 47, 225, 71, 246, 21, 42, 20, 5, 12, 120,
              11, 3, 14, 71, 11, 14, 11, 16, 90, 1, 16, 52, 95)
bearings <- c(17.88, 28.92, 33.00, 41.52, 42.12, 45.60, 48.40, 51.84, 51.96, 54.12,
              55.56, 67.80, 68.64, 68.64, 68.88, 84.12, 93.12, 98.64, 105.12,
              105.84, 127.92, 128.04, 173.40)
# Tensile strength (GPa) of 69 carbon fibres tested at a gauge length of
# 20 mm (Bader and Priest, 1982), as transcribed by Kundu and Gupta (2006).
# We use the original values, without the subtraction of 0.75 made by those
# authors.
fibres <- c(1.312, 1.314, 1.479, 1.552, 1.700, 1.803, 1.861, 1.865, 1.944, 1.958,
            1.966, 1.997, 2.006, 2.021, 2.027, 2.055, 2.063, 2.098, 2.140, 2.179,
            2.224, 2.240, 2.253, 2.270, 2.272, 2.274, 2.301, 2.301, 2.359, 2.382,
            2.382, 2.426, 2.434, 2.435, 2.478, 2.490, 2.511, 2.514, 2.535, 2.554,
            2.566, 2.570, 2.586, 2.629, 2.633, 2.642, 2.648, 2.684, 2.697, 2.726,
            2.770, 2.773, 2.800, 2.809, 2.818, 2.821, 2.848, 2.880, 2.954, 3.012,
            3.067, 3.084, 3.090, 3.096, 3.128, 3.233, 3.433, 3.585, 3.585)
stopifnot(length(aircon) == 30, length(bearings) == 23, length(fibres) == 69)

# Weibull MLE through the profile equation for the shape (uniroot), and
# log-normal MLE
fit_weibull <- function(x) {
  y <- log(x)
  h <- function(k) { w <- exp(k * (y - max(y))); sum(w * y) / sum(w) - 1 / k - mean(y) }
  k   <- uniroot(h, c(1e-3, 100), tol = 1e-12)$root
  lam <- mean(x^k)^(1 / k)
  c(shape = k, scale = lam, loglik = sum(dweibull(x, shape = k, scale = lam, log = TRUE)))
}
fit_lognormal <- function(x) {
  y <- log(x); mu <- mean(y); s <- sqrt(mean((y - mu)^2))
  c(mu = mu, sigma = s, loglik = sum(dlnorm(x, mu, s, log = TRUE)))
}

analyse <- function(x, name) {
  n  <- length(x); y <- log(x)
  k2 <- mean((y - mean(y))^2)
  k3 <- mean((y - mean(y))^3)
  Y  <- matrix(y, nrow = 1)
  T3 <- stat_T3(Y)
  t0 <- null_tabs[[as.character(n)]]      # null tabulation with the seeds of Table 1
  fw <- fit_weibull(x); fl <- fit_lognormal(x)
  r  <- kappenman_ratio(Y)
  data.frame(
    data_set = name, n = n, kappa2 = k2, kappa3 = k3,
    kappa3_corrected = correction_factor(n) * k3, gamma1 = k3 / k2^1.5, T3 = T3,
    p_normal = pnorm(T3), p_monte_carlo = mean(t0 <= T3),
    q005_monte_carlo = quantile(t0, ALPHA, names = FALSE),
    decision_T3 = ifelse(T3 <= Z_ALPHA, "Weibull", "log-normal"),
    loglik_W = unname(fw["loglik"]), loglik_LN = unname(fl["loglik"]),
    shape_W = unname(fw["shape"]), scale_W = unname(fw["scale"]),
    mu_LN = unname(fl["mu"]), sigma_LN = unname(fl["sigma"]),
    decision_RML = ifelse(unname(fw["loglik"]) > unname(fl["loglik"]), "Weibull", "log-normal"),
    r_Kappenman = r,
    decision_Kappenman = ifelse(r > KAPPENMAN_CUTOFF, "log-normal", "Weibull"))
}
tab_appl <- rbind(analyse(aircon, "air conditioning"), analyse(bearings, "ball bearings"),
                  analyse(fibres, "carbon fibres"))
print(t(tab_appl))

# Figure 3: the three data sets in the log-cumulant diagram. Each point is the
# pair (corrected kappa3, kappa2) that enters T3, labelled with the value of
# T3. The curves are those of Figure 1, drawn in black, with the Weibull solid
# and the log-normal dashed. Panel (a) shows the full diagram and marks, with
# a dotted rectangle, the region near the origin that is enlarged in panel
# (b), where data sets 2 and 3 lie.
col_data   <- c("#D62728", "#2CA02C", "#1F77B4")   # data sets 1, 2 and 3: red, green and blue
pch_data   <- c(19, 17, 15)
data_names <- c("Set 1: air conditioning", "Set 2: ball bearings", "Set 3: carbon fibres")
px <- tab_appl$kappa3_corrected; py <- tab_appl$kappa2
label_T3 <- function(v)                              # "T3 = -0.735" with a true minus sign
  if (v < 0) bquote(tilde(T)[3] == -.(sprintf("%.3f", -v))) else bquote(tilde(T)[3] == .(sprintf("%.3f", v)))
limits <- list(a = list(x = c(-3.3, 0.25), y = c(0, 2)),
               b = list(x = c(-0.08, 0.045), y = c(0, 0.3)))  # region enlarged in (b)
pdf("fig_data_diagram.pdf", width = 9.2, height = 4.6)
par(mfrow = c(1, 2), mar = c(4.3, 4.6, 2.2, 1))
for (panel in c("a", "b")) {
  lim <- limits[[panel]]
  plot(NA, xlim = lim$x, ylim = lim$y, las = 1, main = paste0("(", panel, ")"),
       xaxt = "n", xlab = expression(kappa[3] ~~ "(third log-cumulant)"),
       ylab = expression(kappa[2] ~~ "(second log-cumulant)"))
  axis(1, at = if (panel == "a") -3:0 else seq(-0.08, 0.04, by = 0.04))
  k2 <- seq(0, 1.1 * lim$y[2], length.out = 1000)
  lines(gamma1_W * k2^1.5, k2, lwd = 1.5)                   # Weibull: kappa3 = gamma1_W kappa2^(3/2)
  segments(0, 0, 0, 1.1 * lim$y[2], lwd = 1.5, lty = 2)    # log-normal: kappa3 = 0
  inside <- which(px >= lim$x[1] & px <= lim$x[2] & py >= lim$y[1] & py <= lim$y[2])
  points(px[inside], py[inside], pch = pch_data[inside], col = col_data[inside],
         cex = if (panel == "a") 1.3 else 1.1)
  if (panel == "a") {
    with(limits$b, rect(x[1], y[1], x[2], y[2], border = "grey40", lty = 3))
    text(limits$b$x[2], limits$b$y[2], "(b)", pos = 4, col = "grey40", cex = 0.85)
    text(px[1], py[1], label_T3(tab_appl$T3[1]), pos = 2, col = col_data[1], cex = 1)
    legend("bottomleft", c("Weibull", "Log-normal", data_names), lty = c(1, 2, NA, NA, NA),
           lwd = c(1.5, 1.5, NA, NA, NA), pch = c(NA, NA, pch_data),
           col = c("black", "black", col_data), bty = "n", cex = 0.9)
  } else {
    # data set 2 to the left of its point, data set 3 to the right of the log-normal half-line
    text(px[2], py[2], label_T3(tab_appl$T3[2]), pos = 2, col = col_data[2], cex = 1)
    text(0.004, py[3], label_T3(tab_appl$T3[3]), adj = c(0, 0.5), col = col_data[3], cex = 1)
  }
}
dev.off()


# =============================================================================
# 6. OUTPUTS
# Workbook with bold headers, frozen first row and automatic column widths
# (CSV files if the openxlsx package is not installed), and rows of Tables 1
# and 2 ready for LaTeX.
# =============================================================================
save_results <- function(tables, file) {
  if (requireNamespace("openxlsx", quietly = TRUE)) {
    wb   <- openxlsx::createWorkbook()
    bold <- openxlsx::createStyle(textDecoration = "bold")
    for (nm in names(tables)) {
      openxlsx::addWorksheet(wb, nm)
      openxlsx::writeData(wb, nm, tables[[nm]], headerStyle = bold)
      openxlsx::freezePane(wb, nm, firstRow = TRUE)
      openxlsx::setColWidths(wb, nm, cols = seq_len(ncol(tables[[nm]])), widths = "auto")
    }
    openxlsx::saveWorkbook(wb, file, overwrite = TRUE)
  } else {
    for (nm in names(tables))
      write.csv(tables[[nm]], paste0(sub("\\.xlsx$", "", file), "_", nm, ".csv"), row.names = FALSE)
  }
}
save_results(list(constants = constants, correction = correction, calibration = tab1,
                  min_sample_size = min_n, power_scan = power_scan,
                  pcs = tab2, applications = tab_appl),
             "article_results.xlsx")

f3 <- function(v) sprintf("%.3f", v)

# Table 2 of the article: one row per n, with T3, RML and Kappenman (r) under
# each family. Each cell shows the smallest and the largest value across the
# parameters of the family (with the same seeds they coincide, and the cell
# shows a single value). In each row, the largest PCS of each family is set
# in bold.
cell <- function(v) {
  a <- f3(min(v)); b <- f3(max(v))
  if (a == b) a else paste0(a, "--", b)
}
rows_tab2 <- sapply(ns_tab2, function(n) {
  blocks <- lapply(c("Weibull", "Log-normal"), function(fam) {
    s   <- tab2[tab2$family == fam & tab2$n == n, ]
    v   <- list(s$PCS_T3, s$PCS_RML, s$PCS_Kappenman)
    txt <- sapply(v, cell)
    top <- sapply(v, function(x) round(max(x), 3))
    ifelse(top == max(top), paste0("\\textbf{", txt, "}"), txt)
  })
  sprintf("%3d & %s \\\\", n, paste(unlist(blocks), collapse = " & "))
})
# in Table 1 the quantile goes in math mode, so that LaTeX prints a true minus sign
out_lines <- c("% Table 1",
               sprintf("%3d & $%s$ & %s & %s \\\\", tab1$n, f3(tab1$q005), f3(tab1$size), f3(tab1$power)),
               "", "% Table 2 (n; Weibull: T3, RML, r; log-normal: T3, RML, r)",
               rows_tab2)
writeLines(out_lines, "latex_tables.tex")
