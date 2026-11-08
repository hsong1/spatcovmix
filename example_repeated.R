remove(list = ls())
set.seed(1312)
# example
n0 = 25
coords <- expand.grid(x1 = runif(n0), x2 = runif(n0))
n = nrow(coords)

# assume each location has 2 repeats; j=1,2
n_repeats = 2

# set up mu_j = X(s)_j %*% beta
location_id <- rep(1:n, each=n_repeats)
X_s <- cbind(1, coords$x1, coords$x2, coords$x1^2, coords$x2^2, coords$x1 * coords$x2)
X_a <- rnorm(n*n_repeats,mean = 1,sd = 1)
X <- cbind(X_s[location_id,,drop=F],X_a)
colnames(X) <- c("intercept", "s1", "s2", "s1^2", "s2^2", "s1*s2","f")
beta0 <- rnorm(ncol(X)) # true coefficients
mu0 <- as.vector(X %*% beta0)

#  covariance C_Z(h) = w_0 * I(h == 0) + sum_s w_s  * exp(-(h / tau_s)^2)
coords_aug <- coords[location_id,,drop=F]
distMat <- as.matrix(dist(coords_aug))
s_h <- max(distMat)

# support points for spatial component: zero, short-range, long-range
true_support <- c(0, 0.1 * s_h, 0.4 * s_h)
true_support_iter <- c(true_support, -1) # for adding measurement error variance
spatial_weights <- c(0.5, 1.5, 0.5) # spatial weights: theta_0, short, long
sigma_e2 = 1 # measurement error variance
w0 = c(spatial_weights, sigma_e2)

Sigma_true <- makeSigma_mLS_r(
  beta = w0,
  inds = 1:length(w0),
  allSigma_i = lapply(1:length(w0), function(i) {
    makeSigma_i(distMat / s_h, true_support_iter[i] / s_h)
  })
)

cat("spatial component: (true support:", true_support, ", true weights:", spatial_weights, ")\n")
cat("measurement error variance:", sigma_e2, "\n")

# Y ~ N(mu,Sigma_true)
library(mvnfast)
Y <- as.vector(rmvn(n = 1, mu = mu0, sigma = Sigma_true))

###### train & test split ######
uloc_id <- unique(location_id)
train_loc_id <- sample(uloc_id ,round(0.8*n))
test_loc_id <- uloc_id[!uloc_id %in% train_loc_id]

coords_obs <- coords_aug[location_id %in% train_loc_id, ]
X_obs <- X[location_id %in% train_loc_id, , drop = F]
Y_obs <- Y[location_id %in% train_loc_id]

coords_new <- coords_aug[location_id %in% test_loc_id, ]
X_new <- X[location_id %in% test_loc_id, , drop = F]
Y_new <- Y[location_id %in% test_loc_id]

###### fitting ######
fit_reml <- cov_SR(coords = coords_obs, Y = Y_obs, X = X_obs, method = "REML", lenAlphaGrid = 30, print = TRUE)
fit_ml <- cov_SR(coords = coords_obs, Y = Y_obs, X = X_obs, method = "ML", lenAlphaGrid = 30, print = TRUE)

fit_LS <- lm(Y_obs ~ X_obs - 1)
r_obs = resid(fit_LS) # regress out fixed effects
fit_mz <- cov_SR(coords = coords_obs, Y = r_obs, method = "MeanZero", lenAlphaGrid = 30, print = TRUE)
# test S input:
fit_mz_S <- cov_SR(coords = coords_obs, S = tcrossprod(r_obs), method = "MeanZero", lenAlphaGrid = 30, print = TRUE)

###### estimated fixed effects ######
# example for fit_ml
data.frame(
  true = beta0,
  coef = fit_ml$regression_fit$regression_coef,
  se = sqrt(diag(fit_ml$regression_fit$regression_vcov))) |> 
  dplyr::mutate(z = coef/se)

###### compare results ####
h_vals <- seq(0, s_h, length.out = 200)

# True C_Z(h) for spatial component Z
C_true <- evalC_mLS(h_vals, weights = spatial_weights, support = true_support)

# evalC_mLS computes sum_j w_j * exp(-(h / alpha_j)^2)
C_reml <- evalC_mLS(dist = h_vals, weights = fit_reml$mu$weights, support = fit_reml$mu$support)
C_ml <- evalC_mLS(dist = h_vals, weights = fit_ml$mu$weights, support = fit_ml$mu$support)
C_mz <- evalC_mLS(dist = h_vals, weights = fit_mz$mu$weights, support = fit_mz$mu$support)

result <- data.frame(h = h_vals, True = C_true, REML = C_reml, ML = C_ml, MeanZero = C_mz)

rmse <- function(est, tru) sqrt(mean((est - tru)^2))

cat("\nRMSE of estimated C_Z(h):\n")
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
  labs(title = "cov_SR: Estimated vs. True Covariance Function for Spatial Component", x = "Distance  h", y = "C(h)") +
  theme_bw(base_size = 13) +
  theme(legend.position = "right")

print(p)


Ypred_reml <- krigingPred_mLS(
  weights = fit_reml$mu$weights,
  support = fit_reml$mu$support,
  sigma_e2 = fit_reml$sigma_e2,
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
  sigma_e2 = fit_ml$sigma_e2,
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
  sigma_e2 = fit_mz$sigma_e2,
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
