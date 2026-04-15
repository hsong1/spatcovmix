remove(list = ls())
set.seed(1312)
# example
n0 = 25
coords <- expand.grid(x1 = runif(n0), x2 = runif(n0))
n = nrow(coords)
n_obs = round(0.9 * n)
n_new = n - n_obs

# set up mu = X(s) %*% beta
X <- cbind(1, coords$x1, coords$x2, coords$x1^2, coords$x2^2, coords$x1 * coords$x2)
colnames(X) <- c("intercept", "s1", "s2", "s1^2", "s2^2", "s1*s2")
beta0 <- rnorm(ncol(X)) # true coefficients
mu0 <- as.vector(X %*% beta0)

# X <- matrix(rep(1, n), nrow = n, ncol = 1)
# beta0 = 3
# mu0 <- X * beta0

#  covariance C(h) = w_nug * I(h == 0) + sum_s w_s  * exp(-(h / tau_s)^2)
distMat <- as.matrix(dist(coords))
s_h <- max(distMat)

# support points: nugget, short-range, long-range
true_support <- c(0, 0.1 * s_h, 0.4 * s_h)
w0 <- c(0.5, 1.5, 0.5) # weights: nugget, short, long

Sigma_true <- makeSigma_mLS_r(
  beta = w0,
  inds = 1:length(w0),
  allSigma_i = lapply(1:length(w0), function(i) {
    makeSigma_i(distMat / s_h, true_support[i] / s_h)
  })
)

cat("true support:", true_support, ", true weights:", w0, "\n")

# Y ~ N(mu,Sigma_true)
library(mvnfast)
Y <- as.vector(rmvn(n = 1, mu = mu0, sigma = Sigma_true))

###### train & test split ######
coords_obs <- coords[1:n_obs, ]
X_obs <- X[1:n_obs, , drop = F]
Y_obs <- Y[1:n_obs]

coords_new <- coords[(n_obs + 1):n, ]
X_new <- X[(n_obs + 1):n, , drop = F]
Y_new <- Y[(n_obs + 1):n]

###### fitting ######
fit_reml <- cov_SR(coords = coords_obs, Y = Y_obs, X = X_obs, method = "REML", lenAlphaGrid = 30, print = TRUE)
fit_ml <- cov_SR(coords = coords_obs, Y = Y_obs, X = X_obs, method = "ML", lenAlphaGrid = 30, print = TRUE)

fit_LS <- lm(Y_obs ~ X_obs - 1)
r_obs = resid(fit_LS) # regress out fixed effects
fit_mz <- cov_SR(coords = coords_obs, Y = r_obs, method = "MeanZero", lenAlphaGrid = 30, print = TRUE)
# test S input:
fit_mz_S <- cov_SR(coords = coords_obs, S = tcrossprod(r_obs), method = "MeanZero", lenAlphaGrid = 30, print = TRUE)


###### compare results ####
h_vals <- seq(0, s_h, length.out = 200)

# True C(h): nugget shows only at h = 0
C_true <- evalC_mLS(h_vals, weights = w0, support = true_support)

# evalC_mLS computes sum_j w_j * exp(-(h / alpha_j)^2)
C_reml <- evalC_mLS(dist = h_vals, weights = fit_reml$mu$weights, support = fit_reml$mu$support)
C_ml <- evalC_mLS(dist = h_vals, weights = fit_ml$mu$weights, support = fit_ml$mu$support)
C_mz <- evalC_mLS(dist = h_vals, weights = fit_mz$mu$weights, support = fit_mz$mu$support)

result <- data.frame(h = h_vals, True = C_true, REML = C_reml, ML = C_ml, MeanZero = C_mz)

rmse <- function(est, tru) sqrt(mean((est - tru)^2))

cat("\nRMSE of estimated C(h):\n")
rmse_out <- c(rmse(C_reml, C_true), rmse(C_ml, C_true), rmse(C_mz, C_true))
names(rmse_out) = c("REML", "ML", "MZ")
print(rmse_out)

library(ggplot2)
library(tidyr)

result_long <- pivot_longer(result, cols = -h, names_to = "Method", values_to = "C")

p <- ggplot(result_long, aes(h, C, color = Method, linetype = Method)) +
  geom_line(linewidth = 1) +
  scale_linetype_manual(
    values = c(True = "solid", REML = "dashed", ML = "dotdash", MeanZero = "dotted")
  ) +
  labs(title = "cov_SR: Estimated vs. True Covariance Function", x = "Distance  h", y = "C(h)") +
  theme_bw(base_size = 13) +
  theme(legend.position = "right")

print(p)


Ypred_reml <- krigingPred_mLS(
  weights = fit_reml$mu$weights,
  support = fit_reml$mu$support,
  method = "REML",
  coords_obs = coords_obs,
  Y_obs = Y_obs,
  X_obs = X_obs,
  coords_new = coords_new,
  X_new = X_new
)

Ypred_ml <- krigingPred_mLS(
  weights = fit_ml$mu$weights,
  support = fit_ml$mu$support,
  method = "ML",
  coords_obs = coords_obs,
  Y_obs = Y_obs,
  X_obs = X_obs,
  coords_new = coords_new,
  X_new = X_new
)

Ypred_mz <- krigingPred_mLS(
  weights = fit_mz$mu$weights,
  support = fit_mz$mu$support,
  method = "MeanZero",
  coords_obs = coords_obs,
  Y_obs = r_obs,
  X_obs = X_obs,
  coords_new = coords_new,
  X_new = X_new
)
Ypred_mz$pred <- Ypred_mz$pred + X_new %*% coef(fit_LS)

mspe_out <- c(mean((Ypred_reml$pred - Y_new)^2), mean((Ypred_ml$pred - Y_new)^2), mean((Ypred_mz$pred - Y_new)^2))
names(mspe_out) = c("REML", "ML", "MZ")

cat("\nMSPE of kriging predictions at new locations:\n")
print(mspe_out)
