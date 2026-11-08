# Evaluate the spatial covariance from mixing measure
# With repeats, this returns the c_Z (spatial) part of the covariance 
# Without repeats, this returns the covariance of the observed Y
evalC_mLS <- function(dist, weights, support) {
  
  if (length(weights) == 0) return(dist * 0)
  
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
    # when dist is matrix
    # return Sigma(weights,support)
    out <- makeSigma_mLS_r(
      beta = weights,
      inds = 1:length(weights),
      allSigma_i = lapply(1:length(weights), function(i) makeSigma_i(dist, support[i]))
    )
  }

  return(out)
}

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
  if(any(dist_new_obs ==0)) stop("currently only allow predictions for new locations")

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
  } else {
    invSigma_X <- invSigma_obs %*% X_obs
    Xt_invSigma_X <- crossprod(X_obs, invSigma_X)
    beta_hat <- solve(Xt_invSigma_X, crossprod(invSigma_X, Y_obs))
    muhat_obs = X_obs %*% beta_hat
    muhat_new = X_new %*% beta_hat
  }

  pred <- as.numeric(muhat_new + C_new_obs %*% invSigma_obs %*% (Y_obs - muhat_obs))
  return(list(beta_hat = beta_hat, muhat_new = muhat_new, pred = pred))
}
