
jeffreyPrior <- function(x) {
  return( sqrt( (1/(x^2)) + (1/(2 - (exp(x) - exp(-x)))) ) )
}

partOfJeffreyPrior <-function(x) {
  return(2 - exp(x) + exp(-x))
}

plotJeffreyPrior <- function() {
  x = seq(0.1,20, 0.1)
  y = jeffreyPrior(x)
  plot(x,y, type="l")
}

plotPartOfJeffreyPrior <- function() {
  x = seq(0, 20, 0.1)
  y = partOfJeffreyPrior(x)
  plot(x, y, type="l")
}

plotJeffreyPrior()
plotPartOfJeffreyPrior()