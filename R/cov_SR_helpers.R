makeAlphaGrid <- function(minEdist, maxEdist, s_h, len, gridScale = c("log", "equidist"), effThreshold = 1 / 20) {
  transScale <- function(edist) {
    edist / sqrt(-log(effThreshold))
  }
  gridScale <- match.arg(gridScale, c("log", "equidist"))
  if (gridScale == "log") {
    u <- exp(seq(log(transScale(minEdist)), log(transScale(maxEdist)), length.out = len - 1))
  } else {
    u <- seq(transScale(minEdist), transScale(maxEdist), length.out = len - 1)
  }
  return(c(0, u) * s_h)
}

makeSigma_i <- function(distMat, alpha) {
  # (Sigma_i)[l,m] = exp(-h_lm^2/alpha^2)

  if (alpha == 0) {
    out <- 1 * (distMat == 0)
  } else if (alpha >0) {
    out <- exp(-(distMat / alpha)^2)
  } else if (alpha <0){
    # this is for the repeated measure variance
    out <- diag(nrow = dim(distMat)[1]) # I_M = identity matrix with M = nrow(dist_mat)

  }else{
    stop("should not be executed")
  }
  return(out)
}

makeSigma_mLS_r <- function(beta, inds, allSigma_i) {
  # Sigma(theta) = sum_i Sigma(alpha_i) * theta_i
  # beta = theta[theta>0]; inds = which(theta>0)
  stopifnot(length(inds) > 0)
  res <- lapply(1:length(inds), function(j) beta[j] * allSigma_i[[inds[j]]])
  res <- Reduce("+", res)

  return(res)
}

computeMY <- function(method, invSigma_np, X, Y) {
  # compute MY
  # Sigma^{-1}%*% Y for MZ, P_theta %*% Y for ML and REML
  if (method == "MeanZero") {
    return(list(P_theta = NULL, MY = invSigma_np %*% Y, reml_penalty = 0))
  } else {
    # compute P_theta
    invSigma_X <- invSigma_np %*% X
    Xt_invSigma_X <- crossprod(X, invSigma_X)
    P_theta <- invSigma_np - invSigma_X %*% solve(Xt_invSigma_X, t(invSigma_X))

    reml_penalty <-
      if (method == "REML") {
        2 * sum(log(diag(chol(Xt_invSigma_X))))
      } else {
        0
      }

    return(list(P_theta = P_theta, MY = P_theta %*% Y, reml_penalty = reml_penalty))
  }
}

computeTraceMats <- function(method, invSigma_np, X, allSigma_i, P_theta = NULL) {
  # compute all_trace_mats
  # Sigma^{-1}%*%Sigma_i for MeanZero; P_theta%*%Sigma_i for ML and REML (I_R Hessian)
  s = length(allSigma_i)
  if (method == "MeanZero") {
    all_trace_mats <- lapply(1:s, function(i) {
      invSigma_np %*% allSigma_i[[i]]
    })
  } else {
    if (is.null(P_theta)) {
      # compute P_theta
      invSigma_X <- invSigma_np %*% X
      Xt_invSigma_X <- crossprod(X, invSigma_X)
      P_theta <- invSigma_np - invSigma_X %*% solve(Xt_invSigma_X, t(invSigma_X))
    }
    all_trace_mats <- lapply(1:s, function(i) P_theta %*% allSigma_i[[i]])
  }

  return(all_trace_mats)
}


objFun = function(
  method,
  X,
  Y = NULL,
  S = NULL,
  beta = NULL,
  inds = NULL,
  coords = NULL,
  alphaGrid = NULL,
  L = NULL # cholesky decomposition of Sigma_np
) {
  if (is.null(L)) {
    stopifnot(!is.null(coords), !is.null(alphaGrid), !is.null(beta), !is.null(inds))
    distances <- as.matrix(dist(coords))
    s_h <- max(distances)
    allSigma_i <- lapply(alphaGrid, function(a) {
      makeSigma_i(distances / s_h, a / s_h)
    })
    Sigma_np <- makeSigma_mLS_r(beta, inds, allSigma_i)
    L <- chol(Sigma_np)
  }

  invSigma_np <- chol2inv(L)
  logdet <- 2 * sum(log(diag(L)))

  if (!is.null(S)) {
    return(0.5 * (logdet + sum(invSigma_np * S)))
  }

  out <- computeMY(method, invSigma_np, X, Y)
  return(0.5 * (logdet + as.numeric(crossprod(Y, out$MY)) + out$reml_penalty))
}


gradFun <- function(allSigma_i, trace_M, MY = NULL, S = NULL, invSigma_np = NULL) {
  # gradient of negative log-likelihood
  # grad[i] = 0.5 * (tr(trace_M%*%Sigma_i) - (MY)' Sigma_i (MY))
  s <- length(allSigma_i)
  grad <- numeric(s)

  if (!is.null(S)) {
    K <- invSigma_np %*% S %*% invSigma_np
    for (i in 1:s) {
      trace_term <- sum(trace_M * allSigma_i[[i]])
      quad_term <- sum(K * allSigma_i[[i]])
      grad[i] <- 0.5 * (trace_term - quad_term)
    }
  } else {
    for (i in 1:s) {
      trace_term <- sum(trace_M * allSigma_i[[i]])
      quad_term <- as.numeric(crossprod(MY, allSigma_i[[i]] %*% MY))
      grad[i] <- 0.5 * (trace_term - quad_term)
    }
  }
  return(grad)
}


#minimize grad(beta_0)^{\top}(beta-\beta0)+(1/2)(beta-beta0)^\top H(beta-beta0)
# over beta, subject to beta_i\geq 0, for each i
supportReduction_subproblem <- function(grad, XtXFun, beta0, inds, gradTrace = FALSE, tol = 1e-6, maxit = 1000) {
  ### support reduction algorithm ###

  ### Input:
  #grad = vector containing unscaled gradient of negative loglikelihood, evaluated at beta0
  #corresponds to

  #XtXFun = function to compute elements of positive semidefinite H as needed --- H_{ij}=XtXFun(i,j)

  #beta0 = initial point (where grad is computed); contains only the positive (nonzero) entries
  #inds = the location in the alphaGrid vector corresponding to the positive weights in beta

  #note alphaGrid is not passed in explicitly: the length is inferred via length(grad) and the grid details should be dealt with by grad vector computation and the Hessian function XtXFun

  #fill in the needed Hessian matrix elements, standardize Hessian and gradient
  s = length(grad)

  #fill in diagonal of Hessian and elements containing a member of inds
  s_alpha = rep(0, s)

  #rename
  beta = beta0

  XtX = matrix(NA, s, s)
  for (i in 1:s) {
    s_alpha[i] = sqrt(XtXFun(i, i))
    grad[i] = grad[i] / s_alpha[i]
    XtX[i, i] = 1
  }

  for (i in 1:length(inds)) {
    beta[i] = beta[i] * s_alpha[inds[i]]

    for (j in 1:s) {
      if (is.na(XtX[inds[i], j])) {
        XtX[inds[i], j] = XtXFun(inds[i], j) / (s_alpha[inds[i]] * s_alpha[j])
        XtX[j, inds[i]] = XtX[inds[i], j]
      }
    }
  }

  #compute gradient at 0: f(beta)=grad(beta_t)*(beta-beta_t)+(1/2)(beta-beta_t)H(beta-beta_t)
  #f'(0)=grad(beta_t)-H(0-beta_t)=grad(beta_t)-H*beta_t
  #Xtr=negative gradient at beta=0:
  Xtr = -(grad - XtX[, inds, drop = FALSE] %*% beta)

  if (gradTrace) {
    gradientExt = list(gradMin = NULL, gradMax = NULL, gradMinScaled = NULL, gradMaxScaled = NULL)
  } else {
    gradientExt = NULL
  }

  #### Step2 - 4: Loop ####
  count = 0 #iteration count
  convergence = FALSE
  ll1 = NA

  while (count < maxit) {
    count = count + 1

    ### sequential unrestricted minimization and support reduction ###
    SR_conv = FALSE
    if (any(beta < 0)) {
      stop("weights should be always non-negative")
    }
    current = beta #current feasible solution

    while (!SR_conv) {
      ##### unrestricted minimization ####
      if (length(inds) == 0) {
        betaNew = numeric(0) # inds = empty; then solution = null measure
      } else {
        cholesky = FALSE

        #update XtX

        if (!cholesky) {
          H <- XtX[inds, inds, drop = FALSE]
          betaNew <- tryCatch(
            solve(H, Xtr[inds]),
            error = function(e) {
              warning("Singular Fisher information matrix, solve H+tol*I: ", conditionMessage(e))
              solve(H + 1e-6 * diag(nrow(H)), Xtr[inds])
            }
          )
          #betaNew = solve(XtX[inds, inds, drop = FALSE], Xtr[inds])
        } else {
          R = chol(XtX[inds, inds, drop = FALSE])
          Qty = solve(t(R), Xtr[inds])
          betaNew = solve(R, Qty)
        }
      }

      if (all(betaNew > 0)) {
        # if all betaNew > 0, we are good; go to step 2
        SR_conv = TRUE
        beta = betaNew
      } else {
        # current = (c1,...,cK)
        # betaNew = (d1,...,dK)
        # Mass(t) = (1-t)(currentMass) + t newMass ( t=0: old; t=1: new )
        # largest 0<=t<=1 to make Mass(t)>=0 for all K?

        t = -current / (betaNew - current)
        t = t[-length(t)]
        # t[t<0] = -Inf; t[t>1] = -Inf
        # t = t[which.max(t)]
        t[t < 0] = Inf
        t[t > 1] = Inf
        t_star = unique(t[which.min(t)])
        current = (1 - t_star) * current + t_star * betaNew
        inds = inds[-which.min(t)]
        current = current[-which.min(t)]
      }
    }

    ### Step 2: find new support point
    ### 2-1: compute a gradient vector ###
    XtXw = XtX[, inds, drop = FALSE] %*% beta
    gradient = XtXw - Xtr
    if (gradTrace) {
      gradientExt = updateGradient(gradientExt, gradient, s_alpha)
    }

    ### 2-2 decide whether to add support point ###
    gradMin = min(gradient)
    # cat("gradMin at count",count,":",gradMin,"\n")

    #### if all gradients > -tol at current beta, we are done
    if (gradMin > -tol) {
      #changed to: if (gradMin > (-tol*sqrt(norm2_r))){...} ?
      convergence = TRUE
      break
    }

    ### Step 3: extended current support set ###
    ##### otherwise, add the support point where the gradient is largest
    newIndex = which.min(gradient)
    inds = c(inds, which.min(gradient))

    #update Hessian
    for (j in 1:s) {
      if (is.na(XtX[newIndex, j])) {
        XtX[newIndex, j] = XtXFun(newIndex, j) /
          (s_alpha[newIndex] * s_alpha[j])
        XtX[j, newIndex] = XtX[newIndex, j]
      }
    }

    beta = c(beta, 0)
  }

  if (!convergence) {
    warning("SR algorithm did not converge after iter=", maxit, " iterations")
  }

  #rescale beta
  beta = beta / s_alpha[inds]

  return(list(
    beta = beta,
    active_inds = inds,
    iter = count,
    convergence = convergence,
    gradient = gradient,
    gradientExt = gradientExt
  ))
}


updateGradient <- function(gradientExt, grad, scale) {
  gradientExt$gradMin = c(gradientExt$gradMin, min(-grad))
  gradientExt$gradMax = c(gradientExt$gradMax, max(-grad))
  gradientExt$gradMinScaled = c(gradientExt$gradMinScaled, min(-grad * scale))
  gradientExt$gradMaxScaled = c(gradientExt$gradMaxScaled, max(-grad * scale))
  return(gradientExt)
}
