cov_SR <- function(
  coords,
  Y,
  X = matrix(1, nrow = length(Y), ncol = 1),
  method = c("REML", "ML", "MeanZero"),
  alphaGrid = NULL,
  minEdist = 1 / length(Y), # minimum effective distance
  maxEdist = 0.9, # maximum effective distance
  gridScale = c("log", "equidist"),
  effThreshold = 1 / 20,
  lenAlphaGrid = 100,
  tol = 1e-5,
  maxit = 100,
  print = FALSE,
  init = NULL
) {
  method <- match.arg(method, choices = c("REML", "ML", "MeanZero"))
  if (print) {
    cat("Method:", method, "| tol =", tol, "\n\n")
  }

  # check X contains an intercept column (all 1s)
  if (!identical(method, "MeanZero")) {
    has_intercept <- any(apply(X, 2, function(col) all(col == 1)))
    if (!has_intercept) stop("X must contain a column of 1s (intercept).")
  }

  M = length(Y) # nobs
  obj_trace <- numeric(0)

  distances = as.matrix(dist(coords))
  s_h = max(distances)

  # alphaGrid setup
  if (is.null(alphaGrid)) {
    # [0, minEdist*s_h/sqrt(-log(effThreshold)),...,maxEdist*s_h/sqrt(-log(effThreshold)) ]
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
  }

  s = length(alphaGrid) # grid size

  # weight initialization
  if (is.null(init)) {
    # if not provided, pointmass at 0 with weight var(Y)
    inds = which(alphaGrid == 0)
    beta = var(Y)
  } else {
    stopifnot(all(init >= 0))
    stopifnot(length(init) == length(alphaGrid))
    inds = which(init > 0)
    beta = init[inds]
  }

  # precompute Sigma_i for i = 1,...,s
  allSigma_i <- lapply(1:length(alphaGrid), function(i) {
    makeSigma_i(distances / s_h, alphaGrid[i] / s_h)
  })

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

    out <- computeMY(method, invSigma_np, X, Y)
    ll0 <- 0.5 * (logdet + out$reml_penalty + as.numeric(crossprod(Y, out$MY)))
    all_trace_mats <- computeTraceMats(
      method = method,
      invSigma_np = invSigma_np,
      X = X,
      allSigma_i = allSigma_i,
      P_theta = out$P_theta
    ) # for computing I(theta)
    trace_M <- if (method == "REML") out$P_theta else invSigma_np
    grad <- gradFun(MY = out$MY, allSigma_i = allSigma_i, trace_M = trace_M) # current gradient
    XtXFun <- function(i, j) {
      0.5 * sum(all_trace_mats[[i]] * t(all_trace_mats[[j]]))
    } # Hessian functions

    ##### call support reduction subproblem
    res = supportReduction_subproblem(
      grad = grad,
      XtXFun = XtXFun,
      beta = beta,
      inds = inds
    )

    # line search
    w = rep(0, s)
    wNew = rep(0, s)
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
      Sigma_wTemp <- makeSigma_mLS_r(
        wTemp[which(wTemp > 0)],
        which(wTemp > 0),
        allSigma_i
      )
      L_wTemp <- chol(Sigma_wTemp)
      ll1 <- objFun(method = method, X = X, Y = Y, L = L_wTemp)

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
        out <- computeMY(
          method = method,
          invSigma_np = invSigma_wTemp,
          X = X,
          Y = Y
        )
        trace_M <- if (method == "REML") out$P_theta else invSigma_wTemp
        grad <- gradFun(MY = out$MY, allSigma_i = allSigma_i, trace_M = trace_M)
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

        # convergence: relative decrease in objective below tolerance
        #if (abs(-(ll1 - ll0) / M) <= tol) {
        if (max(abs(projGrad) / M) <= tol) {
          conv = TRUE
          mu <- list(support = alphaGrid[w > 0], weights = w[w > 0])

          # return output
          return(list(
            mu = mu,
            inds = inds,
            alphaGrid = alphaGrid,
            alphaGridScaled = alphaGrid / s_h,
            params = w,
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
  mu <- list(support = alphaGrid[w > 0], weights = w[w > 0])

  return(list(
    mu = mu,
    inds = inds,
    alphaGrid = alphaGrid,
    alphaGridScaled = alphaGrid / s_h,
    params = w,
    objective = ll1 / M,
    obj_trace = obj_trace / M,
    projGrad = projGrad / M,
    converged = conv
  ))
}
