# returns unit vector i of length n

ei <- function(i, n) {
  return(ifelse(i == 1:n, 1, 0))
}

# return distancing matrix (i,j) of order n

aij <- function(i, j, n) {
  df <- ei(i, n) - ei(j, n)
  return(outer(df, df))
}

# y is a list of matrices
# returns a matrix with the direct sum of the y[[i]]

directSum <- function(y) {
  mr <- sum(sapply(y, nrow))
  mc <- sum(sapply(y, ncol))
  nr <- c(0, cumsum(sapply(y, nrow)))
  nc <- c(0, cumsum(sapply(y, ncol)))
  p <- length(y)
  z <- matrix(0, mr, mc)
  for (i in 1:p) {
    ind <- (nr[i] + 1):nr[i + 1]
    jnd <- (nc[i] + 1):nc[i + 1]
    z[ind, jnd] <- y[[i]]
  }
  return(z)
}

# y is a matrix
# returns a matrix with the direct sum of p copies of y

directExpand <- function(y, p) {
  return(directSum(rep(list(y), p)))
}

# double centers a symmetric matrix

doubleCenter <- function(x) {
  rs <- apply(x, 1, mean)
  ss <- mean(x)
  return(x - outer(rs, rs, "+") + ss)
}

# mPrint() formats a matrix (or vector, or scalar) of numbers
# for printing 

mPrint <- function(x,
                   digits = 6,
                   width = 8,
                   format = "f",
                   flag = "+") {
  print(noquote(
    formatC(
      x,
      digits = digits,
      width = width,
      format = format,
      flag = flag
    )
  ))
}

# classical MDS

torgerson <- function(delta, p = 2) {
  e <- eigen(-.5 * doubleCenter(as.matrix(delta) ^ 2))
  l <- sqrt(pmax(0, e$values[1:p]))
  if (p == 1) {
    return(as.matrix(e$vectors[, 1] * l))
  } else {
    return(e$vectors[, 1:p] %*% diag(l))
  }
}

# converts dist object w to M-matrix

mmatrix <- function(w) {
  v <- -as.matrix(w)
  diag(v) <- -rowSums(v)
  return(v)
}

# create dist with all elements equal to one

oneDist <- function(n) {
  return(as.dist(matrix(1, n, n)))
}

## partitioned matrix to single matrix

par2mat <- function(x) {
  nc <- length(x[[1]])
  nr <- length(x)
  h <- NULL
  for (i in 1:nr) {
    hi <- NULL
    for (j in 1:nc) {
      hi <- cbind(hi, x[[i]][[j]])
    }
    h <- rbind(h, hi)
  }
  return(h)
}

# vector to partitioned list, inverse of unlist()

vec2List <- function(x, y) {
  z <- list()
  h <- sapply(y, ncol)
  p <- length(h)
  ind <- cumsum(c(1, h))[-(p + 1)]
  jnd <- ind + (h - 1)
  for (k in 1:p) {
    z <- c(z, list(x[ind[k]:jnd[k]]))
  }
  return(z)
}

# tmForm of a matrix X has x_{ij}=0 if j geq i
# while the distances satisfy D(tlForm(X)) = D(X)

tmForm <- function(x) {
  x <- x - matrix(x[1,], nrow(x), ncol(x), byrow = TRUE)
  x[-1,] <- x[-1,] %*% qr.Q(qr(t(x[-1,])))
  return(x)
}

# dlForm of a matrix X has e'X = 0 and x_{ij}=0 if j > i
# while the distances satisfy D(dlForm(X)) = D(X)

dlForm <- function(x) {
  x <- apply(x, 2, function(x)
    x - mean(x))
  return(x %*% qr.Q(qr(t(x))))
}

# returns a configuration scaled to minimize stress

scalConf <- function(delta, x, w = oneDist(attr(delta, "Size"))) {
  d <- dist(x)
  lbd <- sum(w * delta * d) / sum(w * (d ^ 2))
  return(lbd * x)
}

# just returns stress

stressComp <- function(delta, x, w = oneDist(attr(delta, "Size"))) {
  return(sum(w * (delta - dist(x)) ^ 2) / 2)
}

# returns a weighted sum of matrices in a list x with weights in a

listSum <- function(x, a = rep(1, length(x))) {
  n <- length(x)
  y <- array(0, dim(x[[1]]))
  for (i in 1:n) {
    y <- y + a[i] * x[[i]]
  }
  return(y)
}

# returns the Guttman transform of a matrix

guttman <- function(delta, x, w = oneDist(nrow(x))) {
  v <- mmatrix(w)
  b <- mmatrix(w * delta / dist(x))
  return(solve(v + 1, b %*% x))
}

# returns the weighted root-mean-square of a scalar, vector, or matrix

rms <- function(x, w = rep(1, length(x))) {
  return(sqrt(sum(w * (x ^ 2)) / sum(w)))
}

# X is a list of matrices, V a matrix
# return the matrix with elements tr X[[i]]'VX[[j]]

listCrossprod <- function(x, v) {
  m <- length(x)
  h <- matrix(0, m, m)
  for (i in 1:m) {
    for (j in 1:m) {
      h[i, j] <- sum(x[[i]] * (v %*% x[[j]]))
    }
  }
  return(h)
}

last <- function(x) {
  return(x[length(x)])
}

butLast <- function(x) {
  return(x[-length(x)])
}
