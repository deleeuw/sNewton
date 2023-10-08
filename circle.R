library(MASS)

regpol <- function(n) {
  u <- seq(0, 2 * pi, length = (n + 1))
  return(cbind(sin(u), cos(u))[-(n + 1), ])
}


pSmacof <-
  function(delta,
           x,
           y,
           w = 1 - diag(nrow(delta)),
           eps = 1e-10,
           itmax = 1000,
           verbose = TRUE) {
    nx <- nrow(x)
    ny <- nrow(y)
    nt <- nx + ny
    p <- ncol(x)
    v <- -w
    diag(v) <- -rowSums(v)
    vv <- ginv(v)
    x <- apply(x, 2, function(x)
      x - mean(x))
    y <- apply(y, 2, function(x)
      x - mean(x))
    delta <- delta / sqrt(sum(delta ^ 2))
    x <- x / sqrt(sum(x ^ 2))
    y <- y / sqrt(sum(y ^ 2))
    xold <- rbind(x, y / 2)
    dold <- as.matrix(dist(xold))
    sold <- sum(w * (delta - dold) ^ 2)
    itel <- 1
    repeat {
      bold <- -delta / (dold + diag(nt))
      diag(bold) <- -rowSums(bold)
      xgut <- vv %*% bold %*% xold
      dgut <- as.matrix(dist(xgut))
      sgut <- sum(w * (delta - dgut) ^ 2)
      sgat <- 1 - 2 * sum(xgut * (v %*% xgut))
      smit <- 1 + 2 * sum((xgut - xold) * (v %*% (xgut - xold))) - 
        2 * sum(xgut * (v %*% xgut))
      u <- sum(x * xgut[1:nx, ])
      s <- svd(crossprod(y, xgut[nx + 1:ny, ]))
      k <- tcrossprod(s$u, s$v)
      v <- sum(s$d)
      xnew <- rbind(u * x, v * y %*% k)
      smat <- 1 + nt * sum((xgut - xnew) ^ 2) - nt * sum(xgut ^ 2)
      dnew <- as.matrix(dist(xnew))
      snew <- sum((delta - dnew) ^ 2)
      loss <- c(sgut, sgat, smit, smat, sold, snew)
      names(loss) <- c("sgut","sgat","smit","smat","sold","snew")
      print(loss)
      if (itel == itmax) {
        break
      }
      itel <- itel + 1
      xold <- xnew
      dold <- dnew
      sold <- snew
    }
    return(list(
      x = xnew,
      d = dnew,
      s = snew,
      itel = itel
    ))
  }
