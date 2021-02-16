mse_split <- function(data, model) 
  as.data.frame(t(sapply(data, function(obj) metamisc:::mse(p = obj[, model], y = obj$Y))))

rema <- function(object) {
  m <- metafor::rma.uni(unlist(object$estimate), unlist(object$variances), method = "REML", test="knha")
  p <- predict(m)
  c("estimate" = m$b,
    "ci.lb" = m$ci.lb,
    "ci.ub" = m$ci.ub,
    "se" = m$se,
    "tau2" = m$tau2,
    "pi.lb" = p$cr.lb,
    "pi.ub" = p$cr.ub)
}

