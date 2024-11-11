subsTable <- function(leaps, scale = "adjr2"){
  x <- summary(leaps)
  m <- cbind(round(x[[scale]], 3), x$which[, -1])
  colnames(m)[1] <- scale
  m[order(m[,1], decreasing = T),]
}