begin <- proc.time()
library(metamisc)
library(metafor)
source("sampling.R")
source("validation.R")


nsim <- 1000
set.seed(202101)

coef_fw <- coef_bw <- list()
seed <- matrix(ncol = length(.Random.seed), nrow = nsim)

rema12 <- function(object, ...) {
  metamisc:::rema(object, method = "DL", lambda = 1/2)
}


perf_fw <- perf_bw <-  data.frame(matrix(ncol = 7, nrow = nsim))
colnames(perf_bw) <- colnames(perf_fw) <- c("sim", "model", "measure", 
                                            "estimate", "ci.lb", "ci.ub", "se")

perf_fw_rema <- perf_bw_rema <- data.frame(matrix(ncol = 10, nrow = nsim))
colnames(perf_fw_rema) <- colnames(perf_bw_rema) <- c("sim", "model", "measure", 
                                                      "estimate", "ci.lb", "ci.ub", "se",
                                                      "tau2", "pi.lb", "pi.ub")

f <- Y ~ 1
s <- Y ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9
# f <- dvt ~ sex + malign + surg + ddimdich +
#   notraum  + vein + calfdif3 + I(log(durat)) + 
#   age25 + I(age25^2)  + surg:notraum + age25:malign + ddimdich:sex
# d <- sample_clusters( args below)
# str(d)
# this translates into the following (for names see sim_parameter$meta_coefs)
s <- Y ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + I(X8^2)  + X3:X4 + X8:X2 + X7:X1


for (sim in seq_len(nsim)) {
  seed[sim, ] <- .Random.seed
  save(seed, file = "output/nonlinear/seed.RData")

  d <- sample_clusters(n = sim_parameters$n,
                       k = length(sim_parameters$n),
                       a_mean = sim_parameters$meta_coefs["(Intercept)"],
                       a_sd = sim_parameters$meta_tau["(Intercept)"],
                       b_bin_mean = sim_parameters$meta_coefs[c("sex", "malign", "surg", "notraum", "vein", "calfdif3", "ddimdich")],
                       b_norm_mean = sim_parameters$meta_coefs[c("age", "durat")],
                       b_bin_sd = sim_parameters$meta_tau[c("sex", "malign", "surg", "notraum", "vein", "calfdif3", "ddimdich")],
                       b_norm_sd = sim_parameters$meta_tau[c("age", "durat")],
                       covariance = sim_parameters$optim_cor,
                       p_x_bin = sim_parameters$props_studyspecific)
  
  m_fw <- metapred(data = d, strata = "cluster", formula = f, scope = s, estFUN = "glm",
                   center = TRUE, meta.method = "DL", genFUN = "rema12")
  m_bw <- update(m_fw, formula = s, scope = f)
  
  coef_fw[[length(coef_fw) + 1]] <- coef(m_fw)
  coef_bw[[length(coef_bw) + 1]] <- coef(m_bw)
  
  save(coef_fw, file = "output/nonlinear/coef_fw.RData")
  save(coef_bw, file = "output/nonlinear/coef_bw.RData")
  
  v <- sample_clusters(n = sim_parameters$n,
                       k = length(sim_parameters$n),
                       a_mean = sim_parameters$meta_coefs["(Intercept)"],
                       a_sd = sim_parameters$meta_tau["(Intercept)"],
                       b_bin_mean = sim_parameters$meta_coefs[c("sex", "malign", "surg", "notraum", "vein", "calfdif3", "ddimdich")],
                       b_norm_mean = sim_parameters$meta_coefs[c("age", "durat")],
                       b_bin_sd = sim_parameters$meta_tau[c("sex", "malign", "surg", "notraum", "vein", "calfdif3", "ddimdich")],
                       b_norm_sd = sim_parameters$meta_tau[c("age", "durat")],
                       covariance = sim_parameters$optim_cor,
                       p_x_bin = sim_parameters$props_studyspecific)
  
  p_fw <- predict(m_fw, newdata = v)
  p_bw <- predict(m_bw, newdata = v)
  
  ## Average performance
  # Forward
  mse_fw <- metamisc:::mse(p = p_fw, y = v$Y)
  mse_fw_cilb <- mse_fw$estimate + mse_fw$se * qnorm(.025)
  mse_fw_ciub <- mse_fw$estimate + mse_fw$se * qnorm(.975)
  
  perf_fw[sim, ] <- c(sim, "fw", "mse", mse_fw$estimate, mse_fw_cilb, mse_fw_ciub, mse_fw$se)
  
  # Backward
  mse_bw <- metamisc:::mse(p = p_bw, y = v$Y)
  mse_bw_cilb <- mse_bw$estimate + mse_bw$se * qnorm(.025)
  mse_bw_ciub <- mse_bw$estimate + mse_bw$se * qnorm(.975)
  
  perf_bw[sim, ] <- c(sim, "bw", "mse", mse_bw$estimate, mse_bw_cilb, mse_bw_ciub, mse_bw$se)
  
  save(perf_fw, file = "output/nonlinear/perf_fw.RData")
  save(perf_bw, file = "output/nonlinear/perf_bw.RData")
  
  
  ## Stratified and meta-analyzed performance
  v_with_p <- cbind(v, fw = p_fw, bw = p_bw)
  v_split <- split(v_with_p, v_with_p$cluster)
  
  mse_fw_split <- mse_split(v_split, "fw")
  mse_fw_rema <- rema(mse_fw_split)
  perf_fw_rema[sim, ] <- c(sim, "fw", "mse", mse_fw_rema)
  
  mse_bw_split <- mse_split(v_split, "bw")
  mse_bw_rema <- rema(mse_bw_split)
  perf_bw_rema[sim, ] <- c(sim, "bw", "mse", mse_bw_rema)
  
  save(perf_fw_rema, file = "output/nonlinear/perf_fw_rema.RData")
  save(perf_bw_rema, file = "output/nonlinear/perf_bw_rema.RData")
  
  print(sim)
}

perf <- rbind(perf_fw, perf_bw)
save(perf, file = "output/nonlinear/perf.RData")

perf_rema <- rbind(perf_fw_rema, perf_bw_rema)
save(perf_rema, file = "output/nonlinear/perf_rema.RData")

end <- proc.time()
print(end - begin)


# 216 out of 1000 iterations the models were exactly the same
bw_selected <- sapply(coef_bw, function(x) paste0(sort(names(x)), collapse = ""))
fw_selected <- sapply(coef_fw, function(x) paste0(sort(names(x)), collapse = ""))
sum(bw_selected == fw_selected)

# BW selects 8 predictors on average
# FW 7
length(unlist(coef_bw))/1000
length(unlist(coef_fw))/1000

# , the out-of-sample performance was almost exactly the same
mean(as.numeric(perf_fw_rema$estimate) - as.numeric(perf_bw_rema$estimate))
t.test(as.numeric(perf_fw_rema$estimate), as.numeric(perf_bw_rema$estimate), paired = TRUE)

# and generalizability
mean(as.numeric(perf_fw_rema$tau2) - as.numeric(perf_bw_rema$tau2))
t.test(as.numeric(perf_fw_rema$tau2), as.numeric(perf_bw_rema$tau2), paired = TRUE)
