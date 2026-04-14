evalC_mLS <- function(dist, weights, support) {
  gauss_K <- function(h, weight, support) {
    if (support == 0) {
      weight * 1 * (h == 0)
    } else {
      weight * exp(-(h / support)^2)
    }
  }

  if (is.null(dim(dist))) {
    # when dist is vector
    # return C(h; (weights,support))
    out <- matrix(sapply(1:length(weights), function(j) gauss_K(dist, weights[j], support[j])), ncol = length(weights))

    out <- apply(out, 1, sum)
  } else {
    #when dist is matrix
    # returh Sigma(weights,support)
    out <- makeSigma_mLS_r(
      beta = weights,
      inds = 1:length(weights),
      allSigma_i = lapply(1:length(weights), function(i) makeSigma_i(dist, support[i]))
    )
  }

  return(out)
}

krigingPred_mLS <- function(weights, support, method, coords_obs, Y_obs, X_obs, coords_new, X_new) {
  s_h <- max(as.matrix(dist(coords_obs)))

  # number of observations
  n_obs <- nrow(coords_obs)
  n_new <- nrow(coords_new)

  # distances (n_new+n_obs) x (n_new+n_obs)
  dist_ext <- as.matrix(dist(rbind(coords_new, coords_obs)))
  dist_obs <- dist_ext[(n_new + 1):(n_new + n_obs), (n_new + 1):(n_new + n_obs)]
  #all.equal(dist_ext[1:n_new,1:n_new], as.matrix(dist(coords_new)))
  #all.equal(dist_obs, as.matrix(dist(coords_obs)), check.attributes = FALSE)

  # covariance at observed locations and its inverse
  Sigma_obs <- evalC_mLS(dist_obs, weights, support)
  invSigma <- chol2inv(chol(Sigma_obs))

  # cross-covariance: n_new x n_obs
  C_ext <- evalC_mLS(dist_ext, weights, support)
  C_new_obs <- C_ext[1:n_new, (n_new + 1):(n_new + n_obs)]

  # kriging prediction
  # Y_pred = X_new%*%beta_hat + C_new_obs %*% Sigma^{-1} %*% (Y_obs-X_obs %*% beta_hat)
  if (method == "MeanZero") {
    beta_hat = 0
    muhat_obs = 0
    muhat_new = 0
  } else {
    invSigma_X <- invSigma %*% X_obs
    Xt_invSigma_X <- crossprod(X_obs, invSigma_X)
    beta_hat <- solve(Xt_invSigma_X, crossprod(invSigma_X, Y_obs))
    muhat_obs = X_obs %*% beta_hat
    muhat_new = X_new %*% beta_hat
  }

  pred <- as.numeric(muhat_new + C_new_obs %*% invSigma %*% (Y_obs - muhat_obs))
  return(list(beta_hat = beta_hat, muhat_new = muhat_new, pred = pred))
}
