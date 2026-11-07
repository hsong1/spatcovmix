#'@export
cov_SR <- function(
  coords,
  Y = NULL,
  S = NULL,
  X = NULL,
  method = c("REML", "ML", "MeanZero"),
  alphaGrid = NULL,
  minEdist = NULL, # minimum effective distance
  maxEdist = 0.9, # maximum effective distance
  gridScale = c("log", "equidist"),
  effThreshold = 1 / 20,
  lenAlphaGrid = 100,
  tol = 1e-4,
  maxit = 100,
  print = FALSE,
  init = NULL
) {
  method <- match.arg(method, choices = c("REML", "ML", "MeanZero"))
  if (print) {
    cat("Method:", method, "| tol =", tol, "\n\n")
  }

  # input validation
  if (!is.null(S)) {
    # when S is provided, we ignore Y and X and only use S for estimation
    if (!is.null(Y)) {
      stop("Provide either Y or S, not both")
    }
    if (method != "MeanZero") {
      stop("S can only be used with method = 'MeanZero'")
    }
    S <- as.matrix(S)
    if (nrow(S) != ncol(S)) {
      stop("S must be a square matrix")
    }
    M <- nrow(S)
  } else {
    M <- length(Y)
  }

  if (nrow(coords) != M) {
    stop("nrow(coords) needs to be the same as length(Y)")
  }
  M_loc <- nrow(unique(as.data.frame(coords))) # number of unique locations
  repeated <- (M_loc < M)  

  if (is.null(X)) {
    X <- matrix(1, nrow = M, ncol = 1)
  }
  if (is.null(minEdist)) {
    minEdist <- 1 / M_loc
  }

  # check X contains an intercept column (all 1s)
  if (method != "MeanZero") {
    has_intercept <- any(apply(X, 2, function(col) all(col == 1)))
    if (!has_intercept) stop("X must contain a column of 1s (intercept).")
  }

  obj_trace <- numeric(0)

  # pairwise distances and s_h
  distances = as.matrix(dist(coords))
  s_h = max(distances)
  if (s_h <= 0) {
    stop("coords must contain at least two distinct locations")
  }

  # alphaGrid setup
  if (is.null(alphaGrid)) {
    # alpha_min = first positive-scale kernel drops to effThreshold at h=minEdist*s_h
    # alpha_max = last kernel drops to effThreshold at h=maxEdist*s_h
    # [0, alpha_min = minEdist*s_h/sqrt(-log(effThreshold)),...,alpha_max= maxEdist*s_h/sqrt(-log(effThreshold)) ]
    stopifnot(lenAlphaGrid>=3)
    alphaGrid <- makeAlphaGrid(
      minEdist,
      maxEdist,
      s_h = s_h,
      lenAlphaGrid,
      gridScale,
      effThreshold = effThreshold
    )
  } else {
    # a few checks on alphaGrid
    if (length(unique(alphaGrid)) < length(alphaGrid)) {
      stop("alphaGrid cannot contain duplicates")
    }
    if (!(any(alphaGrid == 0))) {
      stop("alphaGrid must contain 0")
    }

    if (any(alphaGrid < 0)) {
      stop("alphaGrid must contain nonnegative scales")
    }
  }

  # precompute Sigma_i
  # when repeated=T, sigma_e^2 can be estimated separately; last matrix is I_M

  # n_params = length(alphaGrid) when not repeated, 
  # length(alphaGrid)+1 (for sigma_e2) when repeated
  n_params <- if(repeated) length(alphaGrid)+1L else length(alphaGrid)

  # compute allSigma_i
  D_scaled <- distances / s_h
  alpha_scaled <- alphaGrid / s_h
  ascaled_iter = if(repeated) c(alpha_scaled,-1) else alpha_scaled
  allSigma_i <- lapply(ascaled_iter, function(a) makeSigma_i(D_scaled, a))
  sigma_e2_index <- if(repeated) length(allSigma_i) else NULL 

  # weight initialization
  if (is.null(init)) {
    # if not provided, use I_M with weight var(Y)
    inds = if (repeated) sigma_e2_index else which(alphaGrid == 0)
    beta <- if (method == "MeanZero") {
      if (!is.null(S)) sum(diag(S)) / M else sum(Y^2) / M
    } else {
      var(Y)
    }
  } else {
    stopifnot(all(init >= 0))
    stopifnot(length(init) == n_params)
    inds = which(init > 0)
    beta = init[inds]
  }

  # main loop
  iteration = 1
  for (iteration in 1:maxit) {
    if (print) {
      cat("Iteration: ", iteration, "\n")
    }

    # compute Sigma(beta) = sum_i (beta_i * Sigma_i) at current beta
    Sigma_np = makeSigma_mLS_r(beta, inds, allSigma_i)
    # compute Sigma^{-1}(beta)
    L <- chol(Sigma_np)
    invSigma_np <- chol2inv(L)
    logdet <- 2 * sum(log(diag(L)))

    # function value: 0.5*(logdet(Sigma) + reml_penalty + Y'MY)
    # MY = invSigma_np Y for MeanZero, P_theta Y for ML/REML.

    if (!is.null(S)) {
      out <- list(P_theta = NULL, MY = NULL, reml_penalty = 0)
      quad_term <- sum(invSigma_np * S)
    } else {
      out <- computeMY(method, invSigma_np, X, Y)
      quad_term <- as.numeric(crossprod(Y, out$MY))
    }
    ll0 <- 0.5 * (logdet + out$reml_penalty + quad_term)

    # gradient: 0.5*(tr(trace_M %*% Sigma_i) - (MY)'Sigma_i MY) for i = 1,...,n_params
    # trace_M = invSigma_np for MeanZero and ML, P_theta for REML
    trace_M <- if (method == "REML") out$P_theta else invSigma_np
    grad <- gradFun(allSigma_i = allSigma_i, trace_M = trace_M, MY = out$MY, S = S, invSigma_np = invSigma_np) # current gradient

    # hessian
    all_trace_mats <- computeTraceMats(
      method = method,
      invSigma_np = invSigma_np,
      X = X,
      allSigma_i = allSigma_i,
      P_theta = out$P_theta
    ) # for computing I(theta)

    XtXFun <- function(i, j) {
      0.5 * sum(all_trace_mats[[i]] * t(all_trace_mats[[j]]))
    } # Hessian functions

    ##### call support reduction subproblem
    res = supportReduction_subproblem(
      grad = grad,
      XtXFun = XtXFun,
      beta0 = beta,
      inds = inds
    )

    # line search
    w = rep(0, n_params)
    wNew = rep(0, n_params)
    w[inds] = beta # current beta
    wNew[res$active_inds] = res$beta # new beta after supportReduction

    # Armijo backtracking setup
    step = wNew - w # direction
    alpha = 1 # initial step size
    c = 0.5 # sufficient decrease constant
    m = sum(step * grad) # directional derivative at w

    # initialization
    conv = FALSE
    flag = FALSE
    lsCount = 0

    while (!flag) {
      wTemp = w + alpha * step # candidate step

      # compute objective value at wTemp
      Sigma_wTemp <- makeSigma_mLS_r(beta = wTemp[which(wTemp > 0)], inds = which(wTemp > 0), allSigma_i = allSigma_i)
      # if chol(Sigma_wTemp) fails, we set objective to Inf and reduce step size
      L_wTemp <- tryCatch(chol(Sigma_wTemp), error = function(e) NULL)
      if (is.null(L_wTemp)) {
        ll1 <- Inf
      } else {
        ll1 <- objFun(method = method, X = X, Y = Y, S = S, L = L_wTemp)
      }

      # Armijo condition: sufficient decrease from ll0
      if ((ll1 - ll0) <= alpha * c * m) {
        if (lsCount > 0) {
          if (print) {
            cat("Rejected line search steps: ", lsCount, "\n")
          }
        }
        # when line search successful
        flag <- TRUE
        w <- wTemp
        obj_trace <- c(obj_trace, ll1)

        # update inds, beta
        inds <- which(w > 0)
        beta <- w[inds]

        # projGrad
        invSigma_wTemp <- chol2inv(L_wTemp)
        out <- list(P_theta = NULL, MY = NULL, reml_penalty = 0)
        if (is.null(S)) {
          out <- computeMY(method = method, invSigma_np = invSigma_wTemp, X = X, Y = Y)
        }
        trace_M <- if (method == "REML") out$P_theta else invSigma_wTemp
        grad <- gradFun(allSigma_i = allSigma_i, trace_M = trace_M, MY = out$MY, S = S, invSigma_np = invSigma_wTemp)
        projGrad <- grad * (w > 0) + grad * (grad < 0) * (w == 0)

        if (print) {
          # cat(abs((ll1 - ll0) / M), "\n\n")
          cat(
            "obj decrease:",
            abs((ll1 - ll0) / M),
            "| maxGrad:",
            max(abs(projGrad)) / M,
            "\n"
          )
        }

        # convergence: maximum absolute projected gradient per observation below tolerance
        if (max(abs(projGrad) / M) <= tol) {
          conv = TRUE
          mu_Z_weights <- w[1:length(alphaGrid)]
          mu <- list(
            support = alphaGrid[mu_Z_weights > 0], 
            weights = mu_Z_weights[mu_Z_weights > 0])

          # return output
          return(list(
            mu = mu,
            inds = inds,
            alphaGrid = alphaGrid,
            alphaGridScaled = alphaGrid / s_h,
            params = w,
            theta0 = mu_Z_weights[alphaGrid == 0],
            sigma_e2 = if (repeated) w[sigma_e2_index] else NA,
            repeated = repeated,
            objective = ll1 / M,
            obj_trace = obj_trace / M,
            projGrad = projGrad / M,
            converged = conv
          ))
        }
      }

      # step rejected
      alpha = alpha / 2 # reduce stepsize
      lsCount = lsCount + 1 # increment lineSearch Count
    }
  }

  warning("maxit reached without convergence")
  inds <- which(w > 0)
  mu_Z_weights <- w[1:length(alphaGrid)]
  mu <- list(support = alphaGrid[mu_Z_weights > 0], 
    weights = mu_Z_weights[mu_Z_weights > 0])

  return(list(
    mu = mu,
    inds = inds,
    alphaGrid = alphaGrid,
    alphaGridScaled = alphaGrid / s_h,
    params = w,
    theta0 = mu_Z_weights[alphaGrid == 0],
    sigma_e2 = if (repeated) w[sigma_e2_index] else NA,
    repeated = repeated,
    objective = ll1 / M,
    obj_trace = obj_trace / M,
    projGrad = projGrad / M,
    converged = conv
  ))
}
