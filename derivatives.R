# returns the gradient of stress at X as n x p matrix

gradientx <- function(x, b, v) {
  return((v - b) %*% x)
}

# returns the gradient of stress at theta as list of vectors

gradientt <- function(x, b, v, y) {
  gradx <- gradientx(x, b, v)
  return(lapply(1:length(y), function(s) drop(crossprod(y[[s]], gradx[, s]))))
}

# returns the Guttman transform at x as matrix

guttmanx <- function (x, b, vinv) {
  return(vinv %*% b %*% x)
}

# returns the Guttman transform at theta as list of vectors

guttmant <- function(x, b, vinv, y) {
  guttx <- guttmanx(x, b, vinv)
  return(lapply(1:length(y), function(s) drop(crossprod(y[[s]], guttx[, s]))))
}

# returns the hessian of stress at x as partioned matrix

hessianx <- function(x, delta, w) {
  p <- ncol(x)
  n <- nrow(x)
  d <- dist(x)
  aux1 <- w * delta / d
  aux2 <- w * delta / (d ^ 3)
  hesx <- list()
  for (s in 1:p) {
    hessrows <- list()
    di <- as.dist(outer(x[, s], x[, s], "-"))
    for (t in 1:p) {
      dj <- as.dist(outer(x[, t], x[, t], "-"))
      hesspart <- mmatrix(aux2 * di * dj)
      if (s == t) {
        hesspart <- mmatrix(w) - (mmatrix(aux1) - hesspart)
      }
      hessrows <- c(hessrows, list(hesspart))
    }
    hesx <- c(hesx, list(hessrows))
  }
  return(hesx)
}

# returns the hessian of stress at theta as partioned matrix

hessiant <- function(x, delta, w, y) {
  hesx <- hessianx(x, delta, w)
  hest <- hesx
  p <- length(y)
  for (s in 1:p) {
    for (t in 1:p) {
      hest[[s]][[t]] <- crossprod(y[[s]], (hesx[[s]][[t]] %*% y[[t]]))
    }
  }
  return(hest)
}


# computes first, second, and third derivatives (in x space for now).

thirdx <- function(delta, x, w = 1 - diag(nrow(x))) {
  n <- nrow(x)
  p <- ncol(x)
  np <- n * p
  w <- 2 * w / sum(w)
  delta <- delta / sqrt(sum(w * delta ^ 2) / 2)
  d <- as.matrix(dist(x))
  y <- as.vector(x)
  func <- 0.0
  fij <- array(0, np)
  ftot <- array(0, np)
  sij <- array(0, c(np, np))
  stot <- array(0, c(np, np))
  tij <- array(0, c(np, np, np))
  ttot <- array(0, c(np, np, np))
  for (j in 1:(n - 1)) {
    ej <- ifelse(j == 1:n, 1, 0)
    for (i in (j + 1):n) {
      dij <- d[i, j]
      wij <- w[i, j]
      gij <- delta[i, j]
      func <- func + wij * (dij - gij) ^ 2
      ei <- ifelse(i == 1:n, 1, 0)
      aij <- outer(ei - ej, ei - ej)
      a <- directSum(rep(list(aij), p))
      u <- drop(a %*% y)
      for (s in 1:np) {
        fij[s] <- u[s] / dij
        for (t in 1:np) {
          sij[s, t] <- a[s, t] / dij - (u[s] * u[t]) / (dij ^ 3)
          for (r in 1:np) {
            tij[s, t, r] <- (u[s] * u[t] * u[r]) / (dij ^ 5)
            tij[s, t, r] <-
              tij[s, t, r] - ((u[s] * a[t, r]) + (u[t] + a[s, r]) + (u[r] * a[s, t])) / (dij ^ 3)
          }
        }
      }
      ttot <- ttot - w[i, j] * delta[i, j] * tij
      stot <- stot - w[i, j] * delta[i, j] * sij
      ftot <- ftot - w[i, j] * delta[i, j] * fij
      for (s in 1:np) {
        ftot[s] <- ftot[s] + w[i, j] * u[s]
        for (t in 1:np) {
          stot[s, t] <- stot[s, t] + w[i, j] * a[s, t]
        }
      }
    }
  }
  return(list(
    func = func / 2.0,
    ftot = ftot,
    stot = stot,
    ttot = ttot
  ))
}
