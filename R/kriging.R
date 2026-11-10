
krigingPred_mLS <- function(weights, support, sigma_e2 = NULL, method, coords_obs, Y_obs, X_obs, coords_new, X_new) {
   
  # input:
  # (weights, support, sigma_e2): mLS fit from (coords_obs, Y_obs, X_obs)
  # method: "MeanZero", "REML", or "ML"
  # (coords_new, X_new): prediction location

  M_loc <- nrow(unique(as.data.frame(coords_obs))) # number of unique locations in coords_obs
  repeated <- (M_loc < nrow(coords_obs)) # repeated measurements in the observed data
  if(repeated & is.null(sigma_e2)) stop("when repeated, sigma_e2 needs to be provided")
  s_h <- max(as.matrix(dist(coords_obs)))

  # number of observations
  n_obs <- nrow(coords_obs)
  n_new <- nrow(coords_new)

  # distances (n_new+n_obs) x (n_new+n_obs)
  dist_ext <- as.matrix(dist(rbind(coords_new, coords_obs)))
  dist_obs <- dist_ext[(n_new + 1):(n_new + n_obs), (n_new + 1):(n_new + n_obs)]
  #all.equal(dist_ext[1:n_new,1:n_new], as.matrix(dist(coords_new)))
  #all.equal(dist_obs, as.matrix(dist(coords_obs)), check.attributes = FALSE)

  # check overlap between coords_obs and coords_new
  dist_new_obs <- dist_ext[(1:n_new),(n_new + 1):(n_new + n_obs)]  
  if(any(dist_new_obs ==0) & !repeated) stop("without repeated measurements, predictions are only allowed for new locations")

  # covariance at observed locations and its inverse
  Sigma_obs <- evalC_mLS(dist_obs, weights, support) 
  if(repeated) Sigma_obs = Sigma_obs + sigma_e2*diag(n_obs)
  invSigma_obs <- chol2inv(chol(Sigma_obs))

  # cross-covariance: n_new x n_obs
  C_ext <- evalC_mLS(dist_ext, weights, support)
  C_new_obs <- C_ext[1:n_new, (n_new + 1):(n_new + n_obs)]

  # kriging prediction
  # Y_pred = X_new%*%beta_hat + C_new_obs %*% Sigma_obs^{-1} %*% (Y_obs-X_obs %*% beta_hat)
  if (method == "MeanZero") {
    beta_hat = 0
    muhat_obs = 0
    muhat_new = 0
    vcov_beta = NULL                                   
  } else {
    invSigma_X <- invSigma_obs %*% X_obs
    Xt_invSigma_X <- crossprod(X_obs, invSigma_X)
    beta_hat <- solve(Xt_invSigma_X, crossprod(invSigma_X, Y_obs))
    muhat_obs = X_obs %*% beta_hat
    muhat_new = X_new %*% beta_hat
    vcov_beta = solve(Xt_invSigma_X)                   # Var(beta_hat) = (X' Sigma^{-1} X)^{-1}
  }

  pred <- as.numeric(muhat_new + C_new_obs %*% invSigma_obs %*% (Y_obs - muhat_obs))

  #   var_pred = c_Z(0) - C_new_obs Sigma^{-1} C_new_obs'            (kriging variance)
  #            + b0' (X' Sigma^{-1} X)^{-1} b0, b0 = x_new - X' Sigma^{-1} c0   (uncertainty in beta_hat)
  # with c_Z(0) = sum(weights) = theta_1 + sum_{i>=2} theta_i. Only defined with repeats
  if (repeated) {
    A <- C_new_obs %*% invSigma_obs                    # kriging weights, n_new x n_obs
    var_pred <- sum(weights) - rowSums(A * C_new_obs)
    if (method != "MeanZero") {
      B <- X_new - A %*% X_obs                         # rows = b0
      var_pred <- var_pred + rowSums((B %*% vcov_beta) * B)
    }
    var_pred <- as.numeric(pmax(var_pred, 0))
  } else {
    var_pred <- rep(NA_real_, n_new)
  }

  return(list(beta_hat = beta_hat, muhat_new = muhat_new, pred = pred,
              vcov_beta = vcov_beta, var_pred = var_pred))
}
