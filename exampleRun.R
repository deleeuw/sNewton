exampleRun <- function(delta) {
  args <- list(
    list(delta, itmax = 1000000),
    list(delta, relax = TRUE, itmax = 1000000),
    list(delta, basis = "A", itmax = 1000000),
    list(delta, basis = "A", relax = TRUE, itmax = 1000000),
    list(delta, eps3 = -1, itmax = 250),
    list(delta, basis = "A", eps3 = -1, itmax = 250),
    list(delta, strategy = 2),
    list(delta, strategy = 2, basis = "A"),
    list(delta, strategy = 3),
    list(delta, strategy = 3, basis = "A"),
    list(delta, eps3 = 2),
    list(delta, eps3 = 4),
    list(delta, eps3 = 6)
  )
  n <- length(args)
  run <- as.list(1:n)
  trun <- as.list(1:n)
  for (i in 1:n) {
    run[[i]] <- do.call("mSmacof", args[[i]])
    trun[[i]] <-
      suppressWarnings(fivenum(microbenchmark(do.call(
        "mSmacof", args[[i]]
      ), times = 11L)$time / (10 ^ 6)))
  }
  results <- matrix(0, n, 5)
  row.names(results) <-
    c(
      "smacof      ",
      "smacofrelax ",
      "smacofa     ",
      "smacofarelax",
      "newton      ",
      "newtona     ",
      "strategy2   ",
      "strategy2a  ",
      "strategy3   ",
      "strategy3a  ",
      "strategy1-2 ",
      "strategy1-4 ",
      "strategy1-6 "
    )
  colnames(results) <-
    c("itel", "stress", "emprate", "minhess", "time")
  for (i in 1:n) {
    results[i,] <-
      c(run[[i]]$itel,
        run[[i]]$loss,
        run[[i]]$erat,
        min(run[[i]]$evlt),
        trun[[i]][[3]])
  }
  return(list(results = results, run = run))
}

