hat.plot <- function(fit){
  n = length(coefficients(fit))
  p = length(fitted(fit))
  plot(hatvalues(fit), main = "index plot of hat value")
  abline(h = c(2,3) * n/p, col = "red", lyt = 2)
  identify(1:p, hatvalues(fit), names(hatvalues(fit)))
}
