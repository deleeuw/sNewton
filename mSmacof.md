---
title: "A Souped-Up Version of **smacof**"
author: 
- Jan de Leeuw - University of California Los Angeles
date: 'Started September 02 2023, Version of September 10, 2023'
output:
  bookdown::html_document2:
    keep_md: yes
    css: preamble.css
    toc: true
    toc_depth: 3
    number_sections: yes
  bookdown::pdf_document2:
    latex_engine: lualatex 
    includes:
      in_header: preamble.tex
    keep_tex: yes
    toc: true
    toc_depth: 3
    number_sections: yes
graphics: yes
mainfont: Times New Roman
fontsize: 12pt
bibliography: ["mypubs.bib","total.bib"]
abstract: We discuss theory and give R code for a version of **smacof** that (1) allows
  for separate and possibly different linear constraints on different dimensions and that    (2) uses second derivatives to accelerate convergence.
---





**Note:** This is a working paper which will be expanded/updated frequently. All suggestions for improvement are welcome. 

# Introduction

The **smacof** algorithm (@deleeuw_mair_A_09c, @mair_groenen_deleeuw_A_22) computes a (local) minimum and a corresponding (local) minimizer of the (metric, Euclidean) multidimensional scaling (MDS) loss function *stress*, defined as
\begin{equation}
\sigma(X):=\frac12\jis w_{ij}(\delta_{ij}-d_{ij}(X))^2.
(\#eq:sdef)
\end{equation}
(the symbol $:=$ is used for definitions). Of course we hope that the local minimum we find is actually a global minimum, but, except for some special cases, there is no way to guarantee that **smacof** will find the global minimizer.

In equation \@ref(eq:sdef) $W$ and $\Delta$ are known square symmetric and non-negative matrices of order $n$, containing *weights* and
*dissimilarities*, which are the *data* of the MDS data analysis problem. Both $W$ are $\Delta$ are *hollow*, i.e. have zero diagonals. 
The matrix-valued function $D(X)$ also takes values in the symmetric, non-negative, and
hollow matrices of order $n$. Its elements $d_{ij}(X)$ are the *Euclidean distances*
between points $i$ and $j$. 

We minimize stress over the matrix $X$, generally called the *configuration*.
Matrix $X$ has $n$ rows, corresponding with $n$ *objects* that will be represented by
$n$ elements of $\mathbb{R}^p$, i.e. by $n$ *points* in $p$-dimensional space (where 
$p<n$). Thus $X\in\mathbb{R}^{n\times p}$, the linear space of all $n\times p$ matrices. Because distances are invart\iant over translations we will actually only  consider $X$ in the subspace $\overline{\mathbb{R}}^{n\times p}$ of matrices that are *column-centered*. 

# Smacof Notation

Convenient notation for **smacof** problems was introduced by @deleeuw_C_77 and @deleeuw_heiser_C_77. It is used in the vignettes @deleeuw_mair_A_09c and @mair_groenen_deleeuw_A_22, in the review article @groenen_vandevelden_16, as well as in the comprehensive textbook of @borg_groenen_05, and in the textbook in statu nascendi of @deleeuw_B_21. We refer to these publications
for more details about this standard **smacof** notation. 


Next we define, as in traditional **smacof**,
\begin{equation}
A_{ij}:=(e_i-e_j)(e_i-e_j)',
(\#eq:adef)
\end{equation}
as well as
\begin{align}
V&:=\jis  w_{ij}A_{ij}(\#eq:vdef),\\
B(X)&:=\jis w_{ij}\frac{\delta_{ij}}{d_{ij}(X)}(X)A_{ij}.(\#eq:bdef).
\end{align}
Equation \@ref(eq:bdef) assume that $d_{ij}(X)\not= 0$. If $d_{ij}(X)=0$
we replace $\delta_{ij}/d_{ij}(X)$ by zero.

Assume, without loss of generality, that the dissimilarities are normalized
by
\begin{equation}
\frac12\jis w_{ij}\delta_{ij}^2=1.
(\#eq:disnorm)
\end{equation}
Now stress becomes
\begin{equation}
\sigma(X)=1-\text{tr}\ X'B(X)X+\frac12\text{tr}\ X'VX.
(\#eq:sigqf)
\end{equation}

## Using Bases

Suppose $Y_s$ with $s=1,\cdots,p$ are bases for subspaces of $\overline{\mathbb{R}}^n$, the space of all centered vectors with $n$ elements. We suppose, without loss of generality, that $Y_s'VY_s=I$ for all $s$. In the rest of this paper we impose linear constraints on the columns of $X$ of the form $x_s=Y_s\theta_s$. Observe that this  eliminates translational invariance from the original formulation of the MDS problem, because we work entirely in $\overline{\mathbb{R}}^n$.

In the terminology of @takane_kiers_deleeuw_A_95 we allow for DCDD (Different sets of Constraints on Different Dimensions). Although it is not within the scope of this paper, DCDD makes it possible to incorporate various linear constraints on the configuration that are discussed in detail in @deleeuw_heiser_C_80, which covers simplexes, circumplexes, and unique variances. 

t is convenient to write stress as a function of the concatenated vector $\theta=(\theta_1,\cdots,\theta_s)$. This gives

$$
\sigma(\theta):=1-2\jis w_{ij}\delta_{ij}d_{ij}(\theta)+\frac12\theta'\theta=1-\theta'B(\theta)\theta+\frac12\theta'\theta,
$$
where
$$
d_{ij}(\theta):=\sqrt{\sum_{s=1}^p \theta_s'Y_s'A_{ij}Y_s\theta_s},
$$

and, correspondingly,
$$
B(\theta):=\jis w_{ij}\frac{\delta_{ij}}{d_{ij}(\theta)}A_{ij}.
$$
Now our task is to minimize stress over $\theta$.

Of course it is possible to choose the $Y_s$ in such a way their span is all of
$\overline{\mathbb{R}}^n$, in which case there are no restrictions on $x_s$.
It is also possible to eliminate rotational indeterminacy, by choosing $Y_s$
as an $n\times (n-p)$ matrix whose first $p-1$ rows are zero.

The iterations of the smacof algorithm are
$$
\theta^{(k+1)}=B(\theta^{(k)})\theta^{(k)}
$$
and the usual qualitative and quantitative convergence analysis (@deleeuw_heiser_C_80, @deleeuw_A_88b) applies.







# Smacof Algorithm

In the vec-torized version of **smacof** we define the **Guttman Transform**
of a configuration $\theta$ as
\begin{equation}
\Gamma(\theta)=V^+B(\theta)\theta,
(\#eq:guttran)
\end{equation}
and the **smacof** majorization (MM) algorithm is simply,
\begin{equation}
x^{(\nu+1)}=\Gamma(x^{(\nu)})
(\#eq:smaalg)
\end{equation}
with $\nu$ indexing iterations. As shown in @deleeuw_C_77, or in more detail in @deleeuw_A_88b, algorithm \@ref(eq:smaalg) converges in the sense that $\sigma(x^{(\nu)})$
is a decreasing sequence, converging to, say, $\sigma_\infty$. Moreover all 
accumulation points of the sequence $x^{(\nu)}$ are stationary points of $\sigma$,
and have the function value $\sigma_\infty$. And finally 
\begin{equation}
\|x^{(\nu+1)}-x^{(\nu)}\|\rightarrow 0,
(\#eq:asyreg)
\end{equation}
which implies that either $x^{(\nu)}$ converges to, say, $x_\infty$ or that the accumulation points of $x^{(\nu)}$ form a continuum (all with the same function value).
For almost all starting points accumulation points will be local minima, and not saddle points or local maxima (in fact, stress has no local maxima, except for one at the origin
$x=0$).


# Theory

$$
d_k(\xi)=\sqrt{\xi'Y'A_kY\xi}
$$
$$
\mathcal{D}d_k^2(\xi)=2Y'A_kY\xi
$$
$$
\mathcal{D}^2d_k^2(\xi)=2Y'A_kY
$$
$$
\mathcal{D}d_k(\xi)=\frac{1}{d_k(\xi)}Y'A_kY\xi
$$
$$
\mathcal{D}^2d_k(\xi)=\frac{1}{d_k(\xi)}Y'\left\{A_k-\frac{A_kY\xi\xi'Y'A_k}{d_k^2(\xi)}\right\}Y
$$

$$
\sigma(\xi)=1-\sum_{k=1}^Kw_k\delta_kd_k(\xi)+\frac12\sum_{k=1}^Kw_kd_k^2(\xi)
$$
$$
\mathcal{D}\sigma(\xi)=-\sum_{k=1}^Kw_k\frac{\delta_k}{d_k(\xi)}Y'A_kY\xi+\sum_{k=1}^Kw_kY'A_kY\xi=Y'(V-B(\xi))Y\xi
$$
$$
\mathcal{D}^2\sigma(\xi)=Y'VY-Y'B(\xi)Y+\sum_{k=1}^Kw_k\frac{\delta_k}{d_k^3(\xi)}Y'A_kY\xi\xi'Y'A_kY=Y'(V-B(\xi)+H(\xi))Y
$$
$Y\xi=x$ $Y'VY=I$


Simplification if
$$
Y=Y_1\oplus\cdots\oplus Y_p
$$
$$
V=\underbrace{V\oplus\cdots\oplus V}_{\text{p times}}
$$
$$
A_k=\underbrace{A_{ij}\oplus\cdots\oplus A_{ij}}_{\text{p times}}
$$

$$
A_{ij}=(e_{i(k)}-e_{j(k)})(e_{i(k)}-e_{j(k)})'
$$
$$
\{H(\xi)\}_{st}=\sum_{k=1}^Kw_k\frac{\delta_k}{d_k^3(\xi)}(x_{i(k)s}-x_{j(k)s})(x_{i(k)t}-x_{j(k)t})
$$
# Appendix: Code

## utils.R


```r
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

# x is a list of matrices
# returns a matrix with the direct sum of the x[[i]]

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

# double centers a symmetric matrix

doubleCenter <- function(x) {
  rs <- apply(x, 1, mean)
  ss <- mean(x)
  return(x - outer(rs, rs, "+") + ss)
}

# mPrint() formats a matrix (or vector, or scalar) of numbers 
# for printing to the console

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

torgerson <- function(delta, p) {
  e <- eigen(-.5 * doubleCenter(delta ^ 2))
  l <- sqrt(pmax(0, e$values[1:p]))
  if (p == 1) {
    return(as.matrix(e$vectors[, 1] * l))
  } else {
    return(e$vectors[, 1:p] %*% diag(l))
  }
}
```

## basis.R


```r
aTriBas <- function(n, p) {
  np <- n - (p - 1)
  x <- matrix(rnorm(np * np), np, np)
  x[, 1] <- 1
  x <- rbind(matrix(0, p - 1, np), x)
  return(qr.Q(qr(x))[, -1, drop = FALSE])
}

aBasis <- function(n, p) {
  return(lapply(1:p, function(s)
    aTriBas(n, s)))
}

bFulBas <- function(n) {
  x <- matrix(rnorm(n * n), n, n)
  x[, 1] <- 1
  return(qr.Q(qr(x))[, -1])
}

bBasis <- function(n, p) {
  x <- bFulBas(n)
  return(lapply(1:p, function(s)
    x))
}
```

## msmacof.R


```r
library(numDeriv)
library(MASS)

source("basis.R")
source("utils.R")

hessian <- function(delta, w, x, y) {
  p <- ncol(x)
  n <- nrow(x)
  d <- as.matrix(dist(x))
  v <- -w
  diag(v) <- -rowSums(v)
  aux1 <- w * delta / (d + diag(n))
  b <- -aux1
  diag(b) <- -rowSums(b)
  aux2 <- w * delta / ((d ^ 3) + diag(n))
  hess <- NULL
  for (i in 1:p) {
    hessrow <- NULL
    for (j in 1:p) {
      hesspart <-
        -aux2 * outer(x[, i], x[, i], "-") * outer(x[, j], x[, j], "-")
      diag(hesspart) <- -rowSums(hesspart)
      if (i == j) {
        hesspart <- v - (b - hesspart)
      }
      hessrow <-
        cbind(hessrow, crossprod(y[[i]], hesspart %*% y[[j]]))
    }
    hess <- rbind(hess, hessrow)
  }
  return(hess)
}

hFunc <- function(theta, delta, w, y, n, p) {
  xmat <- matrix(drop(directSum(y) %*% theta), n, p)
  dmat <- as.matrix(dist(xmat))
  return(sum(w * (delta - dmat) ^ 2) / 2)
}

numHess <- function(delta, w, theta, y) {
  p <- length(y)
  n <- nrow(y[[1]])
  return(numDeriv::hessian(
    hFunc,
    theta,
    n = n,
    p = p,
    delta = delta,
    w = w,
    y = y
  ))
}

mSmacof <-
  function(delta,
           w,
           p = 2,
           xold = torgerson(delta, p),
           basis = "B",
           itmax = 1000,
           burn = itmax,
           eps1 = 15,
           eps2 = 6,
           verbose = TRUE) {
    n <- nrow(xold)
    p <- ncol(xold)
    v <- -w
    diag(v) <- -rowSums(v)
    if (basis == "A") {
      y <- aBasis(n, p)
    } else {
      y <- bBasis(n, p)
    }
    told <- NULL
    for (i in 1:p) {
      yi <- y[[i]]
      vv <- crossprod(yi, v %*% yi)
      ev <- eigen(vv)
      y[[i]] <- yi %*% ev$vectors %*% diag(1 / sqrt(ev$values))
      yi <- y[[i]]
      ti <- drop(crossprod(yi, v %*% xold[, i]))
      told <- c(told, ti)
    }
    xold <- matrix(drop(directSum(y) %*% told), n, p)
    dold <- as.matrix(dist(xold))
    sold <- sum(w * (delta - dold) ^ 2) / 4
    itel <- 1
    repeat {
      bmat <- -w * delta / (dold + diag(n))
      diag(bmat) <- -rowSums(bmat)
      xnew <- xold
      tnew <- drop(directSum(mLtdmL(y, bmat)) %*% told)
      grad <- told - tnew
      grms <- sqrt(sum(grad ^ 2) / p)
      type <- "SMACOF"
      if (itel > burn) {
        hess <- hessian(delta, w, xold, y)
        ness <- numHess(delta, w, told, y)
        hiss <- ginv(hess)
        tnew <- told - drop(hiss %*% grad)
        type <- "NEWTON"
      }
      xnew <- matrix(drop(directSum(y) %*% tnew), n, p)
      grms <- sqrt(grms / p)
      dnew <- as.matrix(dist(xnew))
      snew <- sum(w * (delta - dnew) ^ 2) / 4
      if (verbose) {
        cat(
          type,
          "itel =",
          formatC(itel, digits = 4, format = "d"),
          "sold =",
          formatC(sold, digits = eps1, format = "f"),
          "snew =",
          formatC(snew, digits = eps1, format = "f"),
          "grms =",
          formatC(grms, digits = eps2, format = "f"),
          "\n"
        )
      }
      if ((itel == itmax) ||
          ((abs(sold - snew)) < (10 ^ -eps1)) && (grms < (10 ^ -eps2))) {
        break()
      }
      itel <- itel + 1
      told <- tnew
      xold <- xnew
      dold <- dnew
      sold <- snew
    }
    hess <- hessian(delta, w, xnew, y)
    evlh <- eigen(hess)$values
    evlb <- eigen(bmat)$values
    rate <- rev(1 - evlh)[2]
    return(list(
      xmat = xnew,
      bmat = bmat,
      dmat = dnew,
      itel = itel,
      loss = snew,
      thet = tnew,
      grad = grad,
      hess = hess,
      evlh = evlh,
      evlb = evlb,
      rate = rate
    ))
  }
```

# References
