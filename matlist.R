# y is a list of matrices, x is a list of vectors
# returns the products y[[i]]x[[i]] in a list of vectors

mLvL <- function(y, x) {
  return(lapply(1:length(y), function(i)
    drop(y[[i]] %*% x[[i]])))
}

# y is a list of matrices, x is a list of vectors
# returns the products crossprod(y[[i]],x[[i]]) in a list of vectors

mLtvL <- function(y, x) {
  return(lapply(1:length(y), function(i)
    drop(crossprod(y[[i]], x[[i]]))))
}

# y is a list of matrices, x is a matrix
# returns the products crossprod(y[[i]],x[,i]) as a list of vectors

mLtm <- function(y, x) {
  return(lapply(1:length(y), function(i)
    drop(crossprod(y[[i]], x[, i]))))
}

# y is a list of matrices, a is a matrix
# returns the products crossprod(y[[i]], a %*% y[[i]])
# in a list of matrices

mLtdmL <- function(y, a) {
  return(lapply(1:length(y), function(i)
    crossprod(y[[i]], a %*% y[[i]])))
}

# y is a list of matrices, x is is a list of matrices
# returns the products x[[i %M*% y[[i]] in a list of matrices

mLmL <- function(x, y) {
  return(lapply(1:length(x), function(i)
    x[[i]] %*% y[[i]]))
}

# x is matrix
# returns the columns of x as list of vectors

m2vL <- function(x) {
  return(lapply(1:ncol(x), function(i)
    x[, i]))
}

# x is a list of vectors of the same length
# returns a matrix with the elements of x as columns

vL2m <- function(x) {
  p <- length(x)
  return(sapply(1:p, function(i)
    x[[i]]))
}
