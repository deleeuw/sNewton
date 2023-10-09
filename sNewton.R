source("utils.R")
source("basis.R")
source("exampleRun.R")
source("derivatives.R")

library(MASS)
library(microbenchmark)

mSmacof <-
  function(delta,
           # dissimilarities, class dist
           w = oneDist(attr(delta, "Size")),
           # weights, class dist
           p = 2L,
           # dimensionality
           xold = NULL,
           # initial configuration, matrix, or NULL
           basis = "B",
           # choice of basis
           itmax = 100000L,
           relax = FALSE,
           strategy = 1L,
           eps1 = 15L,
           eps2 = 10L,
           eps3 = 15L,
           verbose = FALSE) {
    n <- attr(delta, "Size")
    w <- w / sum(w)
    v <- mmatrix(w)
    vinv <- solve(w + (1 / n)) - (1 / n)
    delta <- delta / sqrt(sum(w * (delta ^ 2)))
    if (is.null(xold)) {
      xold <- torgerson(delta, p)
    }
    xold <- apply(xold, 2, function(x)
      x - mean(x))
    y <- makeBasis(p, basis, w)
    told <-
      lapply(1:p, function(s)
        drop(crossprod(y[[s]], xold[, s])))
    xold <- sapply(1:p, function(s)
      y[[s]] %*% told[[s]])
    dold <- dist(xold)
    lbd  <- sum(w * delta * dold) / sum(w * (dold ^ 2))
    dold <- dold * lbd
    xold <- xold * lbd
    bold <- mmatrix(w * (delta / dold))
    sold <- sum(w * (delta - dold) ^ 2) / 2
    rold <- Inf
    itel <- 1L
    eopt <- 0.0
    repeat {
      # smacof step
      tnew <-
        lapply(1:p, function(s)
          drop(crossprod(y[[s]], bold %*% xold[, s]))) # Guttman transform at told
      grat <-
        lapply(1:p, function(s)
          told[[s]] - tnew[[s]]) # gradient at told as list
      rgrd <-
        sqrt(sum(sapply(grat, function(x)
          sum(x ^ 2)))) # norm of change in theta
      if (relax) {
        eopt <- min(1, eopt)
        tnew <- lapply(1:p, function(s)
          (1 + eopt) * tnew[[s]] - eopt * told[[s]])
      } # compute relaxed update of theta
      xnew <- sapply(1:p, function(s)
        drop(y[[s]] %*% tnew[[s]])) # compute new x after smacof step
      dnew <-
        dist(xnew) # compute new d after smacof step
      bnew <- mmatrix(w * delta / dnew)
      snew <-
        sum(w * (delta - dnew) ^ 2) / 2 # compute new stress after smacof step
      type <- "SMACOF"
      # newton step
      if ((strategy > 1) || ((strategy == 1) && (rgrd < (10 ^ -eps3)))) {
        hest <-
          hessiant(xold, delta, w, y) # hessian at told
        hist <- ginv(par2mat(hest))
        tnwt <-
          unlist(told) - drop(hist %*% unlist(grat)) # new theta with newton step
        tnwt <- vec2List(tnwt, y)
        xnwt <-
          sapply(1:p, function(s)
            y[[s]] %*% tnwt[[s]]) # compute new x after newton step
        dnwt <- dist(xnwt) # new d after newton step
        bnwt <- mmatrix(w * delta / dnwt) # bmat at tnwt
        snwt <-
          sum(w * (delta - dnwt) ^ 2) / 2 # new stress after newton step
      }
      if (((strategy == 3) && (snwt < snew)) ||
          ((strategy == 2) && (snwt < sold)) ||
          ((strategy == 1) && (rgrd < (10 ^ -eps3)))) {
        snew <- snwt
        tnew <- tnwt
        xnew <- xnwt
        dnew <- dnwt
        bnew <- bnwt
        type <- "NEWTON"
      }
      chng <-
        lapply(1:p, function(s)
          told[[s]] - tnew[[s]])  # change in theta
      rnew <-
        sqrt(sum(sapply(chng, function(x)
          sum(x ^ 2)))) # norm of change in theta
      erat <- rnew / rold # empirical rate
      eopt <- min(1, erat / (2 - erat))
      sdif <- sold - snew
      if (verbose) {
        cat(
          type,
          "itel =",
          formatC(itel, digits = 2, format = "d"),
          "snew =",
          formatC(snew, digits = min(15, eps1), format = "f"),
          "sdif =",
          formatC(sdif, digits = min(15, eps1), format = "f"),
          "chng =",
          formatC(rnew, digits = min(10, eps2), format = "f"),
          "rate =",
          formatC(erat, digits = min(10, eps2), format = "f"),
          "\n"
        )
      }
      if ((itel == itmax) ||
          ((abs(sold - snew)) < (10 ^ -eps1)) &&
          (rnew < (10 ^ -eps2))) {
        break()
      }
      itel <- itel + 1
      told <- tnew
      xold <- xnew
      dold <- dnew
      bold <- bnew
      sold <- snew
      rold <- rnew
    }
    hest <-
      hessiant(xnew, delta, w, y) # hessian after convergence
    hesx <-
      hessianx(xnew, delta, w) # hessian after convergence
    grax <- lapply(1:p, function(s)
      drop(xnew[, s] - vinv %*% bnew %*% xnew[, s])) # gradient at xnew as list
    evlt <- 
      eigen(par2mat(hest))$values # eigenvalues of hessian at t
    evlx <-
      eigen(par2mat(hesx))$values # eigenvalues of hessian at x
    evlb <-
      eigen(ginv(mmatrix(w)) %*% bnew)$values # eigenvalues of V^+B matrix
    if (basis == "A") {
      trat <- rev(1 - evlt)[1] # theoretical convergence rate
    } else {
      trat <- rev(1 - evlt)[2] # theoretical convergence rate
    }
    grat <-
      lapply(1:p, function(s)
        drop(crossprod(y[[s]], grax[[s]]))) # gradient wrt t as list
    grax <- matrix(unlist(grax), n, p) # gradient wrt x as matrix
    if (itel > itmax) {
      warning("Maximum number of iterations reached")
    }
    return(
      list(
        wmat = w,
        delt = delta,
        ybas = y,
        xmat = xnew,
        bmat = bnew,
        dmat = dnew,
        itel = itel,
        loss = snew,
        thet = tnew,
        grax = grax,
        grat = grat,
        hesx = hesx,
        hest = hest,
        evlt = evlt,
        evlx = evlx,
        evlb = evlb,
        erat = erat,
        trat = trat
      )
    )
  }
