

aTriBas <- function(n, p) {
  np <- n - (p - 1)
  x <- diag(np)
  x[, 1] <- 1
  x <- rbind(matrix(0, p - 1, np), x)
  return(qr.Q(qr(x))[, -1, drop = FALSE])
}

aBasis <- function(n, p) {
  return(lapply(1:p, function(s)
    aTriBas(n, s)))
}

bFulBas <- function(n) {
  x <- diag(n)
  x[, 1] <- 1
  return(qr.Q(qr(x))[, -1])
}

bBasis <- function(n, p) {
  x <- bFulBas(n)
  return(lapply(1:p, function(s)
    x))
}

uBasis <- function(n, p) {
  x <- bBasis(n, p)
  for (i in 1:n) {
    ei <- ifelse(i == 1:n, 1, 0)
    ei <- as.matrix(ei - mean(ei))
    x <- c(x, list(ei))
  }
  return(x)
}

makeBasis <- function(p, basis, w) {
  v <- mmatrix(w)
  n <- nrow(v)
  y <-
    switch(basis,
           A = aBasis(n, p),
           B = bBasis(n, p),
           C = cBasis(n, p))
  ev <-
    lapply(1:p, function(s)
      eigen(crossprod(y[[s]], v %*% y[[s]])))
  y <-
    lapply(1:p, function(s)
      y[[s]] %*% ev[[s]]$vectors %*% diag(1 / sqrt(ev[[s]]$values)))
  return(y)
}
