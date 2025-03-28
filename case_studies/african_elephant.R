## Ivory coast elephant case study

## packages
library(mvtnorm)
# devtools::install_github("janoleko/LaMa") # development version
library(LaMa)
# install.packages("scales") # for muted colors
library(scales)


## data
data = read.csv("./data/elephant_data.csv")
data$tod = 2 * data$tod
head(data)
nrow(data)


## colors
color = c("orange", "deepskyblue")



# Initialising IC table ---------------------------------------------------

IC_table <- data.frame(model = 1:5,
                       llk = rep(NA, 5),
                       AIC = rep(NA, 5),
                       BIC = rep(NA, 5))


# Fitting homogeneous and parametric HMMs ---------------------------------

nll_hom = function(par){
  getAll(par, dat)
  ## transforming parameters
  mu = exp(logmu)
  sigma = exp(logsigma)
  kappa = exp(logkappa)
  Gamma = tpm(eta)
  delta = stationary(Gamma)
  ## reporting quantities of interest later
  REPORT(mu)
  REPORT(sigma)
  REPORT(kappa)
  ## computing state-dependent densities
  allprobs = matrix(1, nrow = length(step), ncol = N)
  ind = which(!is.na(step) & !is.na(angle))
  for(j in 1:N){
    allprobs[ind,j] = dgamma2(step[ind],mu[j],sigma[j]) *
      dvm(angle[ind], 0, kappa[j])
  }
  ## forward algorithm
  -forward(delta, Gamma, allprobs)
}

N = 2
par = list(logmu = log(c(0.2, 2)),
           logsigma = log(c(0.2, 2)),
           logkappa = log(c(0.2, 1)),
           eta = rep(-2, 2))

dat = list(step = data$step, angle = data$angle, N = 2)

obj_hom = MakeADFun(nll_hom, par, silent = TRUE)
opt_hom = nlminb(obj_hom$par, obj_hom$fn, obj_hom$gr)

mod_hom = obj_hom$report()

mu = mod_hom$mu
sigma = mod_hom$sigma
kappa = mod_hom$kappa
delta = mod_hom$delta

mod_hom$states = viterbi(mod_hom$delta, mod_hom$Gamma, mod_hom$allprobs)

mod_hom$AIC = 2 * obj_hom$fn() + 2 * length(obj_hom$par)
mod_hom$BIC = 2 * obj_hom$fn() + log(nrow(data)) * length(obj_hom$par)

IC_table[1, 2:4] = c(-obj_hom$fn(), mod_hom$AIC, mod_hom$BIC)

# pdf("./case_studies/figs/elephant_marginal.pdf", width = 8, height = 4)

par(mfrow = c(1,2))
hist(data$step, breaks = 100, prob = T, bor = "white", xlim = c(0,5), main = "", xlab = "step length", ylab = "density")
curve(delta[1]*dgamma2(x, mu[1], sigma[1]), add = T, lwd = 2, col = color[1], n = 500)
curve(delta[2]*dgamma2(x, mu[2], sigma[2]), add = T, lwd = 2, col = color[2], n = 500)
curve(delta[1]*dgamma2(x, mu[1], sigma[1])+delta[2]*dgamma2(x, mu[2], sigma[2]),
      add = T, lwd = 2, lty = 2, n = 500)
legend("topright", legend = c("encamped", "exploratory", "marginal"), col = c(color[1], color[2], "black"),
       lty = c(1,1,2), bty = "n")

hist(data$angle, breaks = 20, prob = T, bor = "white", main = "", xlab = "turning angle", ylab = "density")
curve(delta[1]*dvm(x, 0, kappa[1]), add = T, lwd = 2, col = color[1], n = 500)
curve(delta[2]*dvm(x, 0, kappa[2]), add = T, lwd = 2, col = color[2], n = 500)
curve(delta[1]*dvm(x, 0, kappa[1])+delta[2]*dvm(x, 0, kappa[2]),
      add = T, lwd = 2, lty = 2, n = 500)

# dev.off()




# Univariate smooth models ------------------------------------------------

## penalised likelihood function for either tod or doy
pnll_uni = function(par) {
  getAll(par, dat)

  mu = exp(logmu); REPORT(mu)
  sigma = exp(logsigma); REPORT(sigma)
  kappa = exp(logkappa); REPORT(kappa)

  Gamma = tpm_g(Z, cbind(beta0, beta_spline))
  delta = stationary_p(Gamma[,,1:24], t = 1)

  allprobs = matrix(1, nrow = length(step), ncol = N)
  ind = which(!is.na(step) & !is.na(angle))
  for(j in 1:N){
    allprobs[ind,j] = dgamma2(step[ind], mu[j], sigma[j]) *
      dvm(angle[ind], 0, kappa[j])
  }

  -forward_g(delta, Gamma, allprobs) +
    penalty(beta_spline, S, lambda) # computes quadratic penalty
}



# Tday model --------------------------------------------------------------

k = 12
modmat_tday = make_matrices(~ s(tod, bs = "cc", k = k),
                            data = data,
                            knots = list(tod = c(0, 24))) # telling mgcv where to wrap the basis

Z = modmat_tday$Z
S = modmat_tday$S

# initial parameter values
N = 2
par = list(logmu = log(c(0.3, 1.1)),      # state-dependent step mean
           logsigma = log(c(0.25, 0.75)),  # state-dependent step sd
           logkappa = log(c(0.2, 0.7)),    # state-dependent angle concentration
           beta0 = c(-2,-2),               # state process intercept
           beta_spline = matrix(0, N*(N-1), modmat_tday$pardim$`s(tod)`)) # state process tod smooth

dat = list(step = data$step, angle = data$angle, N = 2,
           Z = Z, S = S,
           lambda = rep(1e4, 2))

system.time(
  mod_tday <- qreml(pnll_uni, par, dat, random = "beta_spline", silent = 0)
)

IC_table[2, 2:4] = c(mod_tday$llk, AIC(mod_tday), BIC(mod_tday))



# Julian model ------------------------------------------------------------

k = 12
modmat_julian <- make_matrices(~ s(doy, bs = "cc", k = k),
                               data = data,
                               knots = list(tod = c(0, 366))) # telling mgcv where to wrap the basis

Z = modmat_julian$Z
S = modmat_julian$S

# initial parameter values
N = 2
par = list(logmu = log(c(0.3, 1.1)),      # state-dependent step mean
           logsigma = log(c(0.25, 0.75)),  # state-dependent step sd
           logkappa = log(c(0.2, 0.7)),    # state-dependent angle concentration
           beta0 = c(-2,-2),               # state process intercept
           beta_spline = matrix(0, N*(N-1), modmat_julian$pardim$`s(doy)`)) # state process tod smooth

dat = list(step = data$step, angle = data$angle, N = 2,
           Z = Z, S = S,
           lambda = rep(1e4, 2))

system.time(
  mod_julian <- qreml(pnll_uni, par, dat, random = "beta_spline", silent = 0)
)

IC_table[3, 2:4] = c(mod_julian$llk, AIC(mod_julian), BIC(mod_julian))



# Additive model ----------------------------------------------------------

# likelihood function with additive effects only
pnll_add = function(par) {
  getAll(par, dat)

  mu = exp(logmu); REPORT(mu)
  sigma = exp(logsigma); REPORT(sigma)
  kappa = exp(logkappa); REPORT(kappa)

  Gamma = tpm_g(Z, cbind(beta0, beta_tod, beta_doy))
  delta = stationary_p(Gamma[,,1:24], t = 1)

  allprobs = matrix(1, nrow = length(step), ncol = N)
  ind = which(!is.na(step) & !is.na(angle))
  for(j in 1:N){
    allprobs[ind,j] = dgamma2(step[ind], mu[j], sigma[j]) *
      dvm(angle[ind], 0, kappa[j])
  }

  -forward_g(delta, Gamma, allprobs) +
    penalty(list(beta_tod, beta_doy), S, lambda) # computes quadratic penalty
}


k = 12
modmat_add <- make_matrices(~ s(tod, bs = "cc", k = k) + s(doy, bs = "cc", k = k),
                               data = data,
                               knots = list(tod = c(0, 24) ,tod = c(0, 366))) # telling mgcv where to wrap the basis

Z = modmat_add$Z
S = modmat_add$S

# initial parameter values
N = 2
par = list(logmu = log(c(0.3, 1.1)),      # state-dependent step mean
           logsigma = log(c(0.25, 0.75)),  # state-dependent step sd
           logkappa = log(c(0.2, 0.7)),    # state-dependent angle concentration
           beta0 = c(-2,-2),               # state process intercept
           beta_tod = matrix(0, N*(N-1), modmat_add$pardim$`s(tod)`),     # state process tod smooth
           beta_doy = matrix(0, N*(N-1), modmat_add$pardim$`s(doy)`))     # state process doy smooth

dat = list(step = data$step, angle = data$angle, N = 2,
           Z = Z, S = S,
           lambda = rep(1e4, 4))

system.time(
  mod_add <- qreml(pnll_add, par, dat,
                   random = c("beta_tod", "beta_doy"), silent = 0)
)

IC_table[4, 2:4] = c(mod_add$llk, AIC(mod_add), BIC(mod_add))


## pseudo res
pres_step <- pseudo_res(data$step,
                        dist = "gamma2",
                        par = list(mean = mod_add$mu, sd = mod_add$sigma),
                        mod = mod_add)
pres_step[is.infinite(pres_step)] <- NA

pres_angle <- pseudo_res(data$angle,
                         dist = "vm",
                         par = list(kappa = mod_add$kappa),
                         mod = mod_add)
pres_angle[is.infinite(pres_angle)] <- NA

par(mfrow = c(2,3))
hist(pres_step, bor = "white", xlim = c(-3,3), prob = TRUE)
curve(dnorm(x), lty = 2, add = TRUE)
qqnorm(na.omit(pres_step), pch = 20, col = "#00000030")
abline(a = 0, b = 1, col = "red")
acf(na.omit(pres_step))

hist(pres_angle, bor = "white", xlim = c(-3,3), prob = TRUE)
curve(dnorm(x), lty = 2, add = TRUE)
qqnorm(na.omit(pres_angle), pch = 20, col = "#00000030")
abline(a = 0, b = 1, col = "red")
acf(na.omit(pres_angle))



# Tensor-product model ----------------------------------------------------

## penalised likelihood function
pnll_ti = function(par) {
  getAll(par, dat)

  mu = exp(logmu); REPORT(mu)
  sigma = exp(logsigma); REPORT(sigma)
  kappa = exp(logkappa); REPORT(kappa)

  Gamma = tpm_g(Z, cbind(beta0, beta_tod, beta_doy, beta_ti))
  delta = stationary_p(Gamma[,,1:24], t = 1)

  allprobs = matrix(1, nrow = length(step), ncol = N)
  ind = which(!is.na(step) & !is.na(angle))
  for(j in 1:N){
    allprobs[ind,j] = dgamma2(step[ind], mu[j], sigma[j]) *
      dvm(angle[ind], 0, kappa[j])
  }

  -forward_g(delta, Gamma, allprobs) +
    penalty2(list(beta_tod, beta_doy, beta_ti), S, lambda) # computes quadratic penalty
}

k = 12
modmat_ti = make_matrices(~ s(tod, bs = "cc", k = k) +
                         s(doy, bs = "cc", k = k) +
                         ti(tod, doy, bs = "cc", k = k),
                         data = data,
                         knots = list(tod = c(0, 24),
                                      doy = c(0, 366))) # telling mgcv where to wrap the basis

Z = modmat_ti$Z
S = modmat_ti$S

# initial parameter values
N = 2
par = list(logmu = log(c(0.3, 1.1)),      # state-dependent step mean
           logsigma = log(c(0.25, 0.75)),  # state-dependent step sd
           logkappa = log(c(0.2, 0.7)),    # state-dependent angle concentration
           beta0 = c(-2,-2),               # state process intercept
           beta_tod = matrix(0, N*(N-1), modmat_ti$pardim$`s(tod)`),     # state process tod smooth
           beta_doy = matrix(0, N*(N-1), modmat_ti$pardim$`s(doy)`),     # state process doy smooth
           beta_ti = matrix(0, N*(N-1), modmat_ti$pardim$`ti(tod,doy)`)) # state process interaction smooth

dat = list(step = data$step, angle = data$angle, N = 2,
           Z = Z, S = S,
           lambda = c(rep(1e4, 4), rep(1e5, 4)))

system.time(
  mod_ti <- qreml(pnll_ti, par, dat,
                  random = c("beta_tod", "beta_doy", "beta_ti"),
                  silent = 0)
)

IC_table[5, 2:4] = c(mod_ti$llk, AIC(mod_ti), BIC(mod_ti))

print(IC_table)

diff_IC_table <- IC_table
diff_IC_table$AIC = IC_table$AIC - IC_table$AIC[5]
diff_IC_table$BIC = IC_table$BIC - IC_table$BIC[4]

round(diff_IC_table, 2)

summary(mod)


## pseudo res
pres_step <- pseudo_res(data$step,
                        dist = "gamma2",
                        par = list(mean = mod_ti$mu, sd = mod_ti$sigma),
                        mod = mod_ti)
pres_step[is.infinite(pres_step)] <- NA

pres_angle <- pseudo_res(data$angle,
                         dist = "vm",
                         par = list(kappa = mod_ti$kappa),
                         mod = mod_ti)
pres_angle[is.infinite(pres_angle)] <- NA

par(mfrow = c(2,3))
hist(pres_step, bor = "white", xlim = c(-3,3), prob = TRUE)
curve(dnorm(x), lty = 2, add = TRUE)
qqnorm(na.omit(pres_step), pch = 20, col = "#00000030")
abline(a = 0, b = 1, col = "red")
acf(pres_step, na.action = na.pass)

hist(pres_angle, bor = "white", xlim = c(-3,3), prob = TRUE)
curve(dnorm(x), lty = 2, add = TRUE)
qqnorm(na.omit(pres_angle), pch = 20, col = "#00000030")
abline(a = 0, b = 1, col = "red")
acf(pres_angle, na.action = na.pass)



sdr <- sdreport_outer(mod, invert = TRUE)
sdr$report

AIC(mod)
BIC(mod)

beta <- mod$beta

## Visualising the fitted model
n_plot = 200
tod_seq = seq(0,24, length = n_plot)
tod_fine = matrix(NA, n_plot, 12)
for(i in 1:n_plot){
  tod_fine[i,] = (tod_seq[i] + 1:12*2) %% 24
}
tod_fine = c(t(tod_fine))
Z_predict = pred_matrix(modmat,
                        newdata = data.frame(tod = rep(tod_fine, 366),
                                             doy = rep(1:366, each = length(tod_fine))))

Deltas = array(dim = c(n_plot, 2, 366))

Z_day_list <- lapply(seq_len(366), function(i) Z_predict[((i-1)*n_plot*12 + 1):(i*n_plot*12), , drop = FALSE])
for(day in 1:366){
  print(day)
  Z_day = Z_day_list[[day]]
  for(i in 1:n_plot){
    thisZ = Z_day[((i-1)*12+1):(i*12),]
    Gamma = tpm_g(thisZ, beta)
    Deltas[i,,day] = stationary_p(Gamma, t = 1)
  }
}

# par(mfrow = c(1,1))
# for(day in 1:366){
#   plot(tod_seq, Deltas[,2,day], type = "l", lwd = 2, col = color[2],
#        xlab = "time of day", ylab = "Pr(state 2)",
#        ylim = c(0,1), bty = "n",
#        main = paste("doy =", day))
#   Sys.sleep(0.1)
# }

# par(mfrow = c(2,3))
# source("./utils/gen_sun_colors.R")
# sun_cycle_colors <- gen_sun_cycle_colors()
# daynames = paste(c("Jan", "March", "May", "Jul", "September", "November"), 1)
# days = c(1, 61, 121, 182, 243, 304)
# for(day in days){
#   plot(NA, xlim = c(0, 24), ylim = c(0,1), bty = "n",
#        xlab = "time of day", ylab = "Pr(state 2)",
#        main = daynames[which(days == day)])
#   polygon(x = c(0, 6.5, 6.5, 0), y = c(0, 0, 1, 1), col = "gray95", border = "white")
#   polygon(x = c(18.5, 24, 24, 18.5), y = c(0, 0, 1, 1), col = "gray95", border = "white")
#   for(t in 0:47){
#     polygon(x = c(t/2, (t+1)/2, (t+1)/2, t/2), y = c(-0.04, -0.04, -0.01, -0.01), col = sun_cycle_colors[t+1], border = sun_cycle_colors[t+1])
#   }
#   lines(tod_seq, Deltas[,2,day], lwd = 2)
# }



# With confidence intervals -----------------------------------------------

## sample parametes from posterior
n_samples <- 2000
set.seed(123)
pars <- rmvnorm(n_samples, mod$par_vec, solve(mod$Hessian_conditional))
pars <- apply(pars, 1, mod$relist_par)
betas <- lapply(pars, function(p) cbind(p$beta0, p$beta_tod, p$beta_doy, p$beta_ti))

allDeltas = array(dim = c(n_plot, 2, 366, n_samples))
for(day in 1:366){
  print(day)
  Z_day = Z_day_list[[day]]
  for(i in 1:n_plot){
    thisZ = Z_day[((i-1)*12+1):(i*12),]
    for(s in 1:n_samples){
      Gamma = tpm_g(thisZ, betas[[s]])
      allDeltas[i,,day,s] = stationary_p(Gamma, t = 1)
    }
  }
}
saveRDS(allDeltas, file = "./case_studies/objects/allDeltas_elephant.rds")

Delta_q <- apply(allDeltas, 1:3, quantile, probs = c(0.025, 0.975))


pdf("./case_studies/figs/elephant_stationary.pdf", width = 7, height = 5)
par(mfrow = c(2,3))
source("./utils/gen_sun_colors.R")
sun_cycle_colors <- gen_sun_cycle_colors()
daynames = paste(c("Jan", "March", "May", "Jul", "September", "November"), 1)
days = c(1, 61, 121, 182, 243, 304)
for(day in days){
  plot(NA, xlim = c(0, 24), ylim = c(0,1), bty = "n",
       xlab = "time of day", ylab = "Pr(state 2)",
       main = daynames[which(days == day)])
  polygon(x = c(0, 6.5, 6.5, 0), y = c(0, 0, 1, 1), col = "gray95", border = "white")
  polygon(x = c(18.5, 24, 24, 18.5), y = c(0, 0, 1, 1), col = "gray95", border = "white")
  for(t in 0:47){
    polygon(x = c(t/2, (t+1)/2, (t+1)/2, t/2), y = c(-0.04, -0.04, -0.01, -0.01), col = sun_cycle_colors[t+1], border = sun_cycle_colors[t+1])
  }
  lines(tod_seq, Deltas[,2,day], lwd = 2)
  polygon(c(tod_seq, rev(tod_seq)),
          c(Delta_q[1,,2,day], rev(Delta_q[2,,2,day])), col = "#00000020", border = NA)
}
dev.off()








### predict tpm over the year
Z_predict = pred_matrix(modmat,
                        newdata = data.frame(tod = rep(tod_seq, 366),
                                             doy = rep(1:366, each = length(tod_seq))))
Z_day_list <- lapply(seq_len(366), function(i) Z_predict[((i-1)*n_plot + 1):(i*n_plot), , drop = FALSE])
Gammas = array(dim = c(2,2,366,n_plot))
for(day in 1:366){
  Gammas[,,day,] = tpm_g(Z_day_list[[day]], beta)
}

par(mfrow = c(1,2))
for(day in 1:366){
  for(i in 1:2){
    for(j in 1:2){
      if(i!=j){
        plot(tod_seq, Gammas[i,j,day,], type = "l", lwd = 2, ylim = c(0,1),
             ylab = paste0("Pr(", i, " -> ", j, ")"), xlab = "time of day", main = day, bty = "n")
      }
    }
  }
  Sys.sleep(0.1)
}


### look at surface
pred_gamma = function(tod, doy, i, j){
  Z_predict = pred_matrix(modmat, newdata = data.frame(tod = tod, doy = doy))
  tpm_g(Z_predict, beta)[i,j,]
}
pred_eta = function(tod, doy, i){
  Z_predict = pred_matrix(modmat, newdata = data.frame(tod = tod, doy = doy))
  Z_predict %*% beta[i,]
}

tod = seq(0, 24, length = 70)
doy = seq(0, 366, length = 70)

gamma_12 = outer(tod, doy, pred_gamma, i = 1, j = 2)
gamma_21 = outer(tod, doy, pred_gamma, i = 2, j = 1)

eta1 = outer(tod, doy, pred_eta, i = 1)
eta2 = outer(tod, doy, pred_eta, i = 2)

library(plot3D)
par(mfrow = c(1,2))
persp3D(x = tod, y = doy, z = gamma_12, theta = -35, phi = 30, main = "Pr(1->2)",
      xlab = "time of day", ylab = "day of year", zlab = "Pr(state 1)")
persp3D(x = tod, y = doy, z = gamma_21, theta = -35, phi = 30, main = "Pr(1->2)",
        xlab = "time of day", ylab = "day of year", zlab = "Pr(state 1)")



library(plotly)
# Create interactive 3D plot using plotly
plot_ly(x = tod, y = doy, z = gamma_12, type = "surface") %>%
  layout(
    title = "Pr(1->2)",
    scene = list(
      xaxis = list(title = "Time of Day"),
      yaxis = list(title = "Day of Year"),
      zaxis = list(title = "Pr(state 1)")
    )
  )

plot_ly(x = tod, y = doy, z = gamma_21, type = "surface") %>%
  layout(
    title = "Pr(1->2)",
    scene = list(
      xaxis = list(title = "Time of Day"),
      yaxis = list(title = "Day of Year"),
      zaxis = list(title = "Pr(state 1)")
    )
  )

# Create interactive 3D plot using plotly
plot_ly(x = tod, y = doy, z = eta1, type = "surface") %>%
  layout(
    title = "Pr(1->2)",
    scene = list(
      xaxis = list(title = "Time of Day"),
      yaxis = list(title = "Day of Year"),
      zaxis = list(title = "Pr(state 1)")
    )
  )

plot_ly(x = tod, y = doy, z = eta2, type = "surface") %>%
  layout(
    title = "Pr(1->2)",
    scene = list(
      xaxis = list(title = "Time of Day"),
      yaxis = list(title = "Day of Year"),
      zaxis = list(title = "Pr(state 1)")
    )
  )

# for(day in 1:365){
#   cat(day, "\n")
#   Z_predict = pred_matrix(modmat,
#                           newdata = data.frame(tod = 1:12*2, doy = day))
#   Gamma = tpm_g(Z_predict, beta)
#   Deltas[,,day] = stationary_p(Gamma)
# }

par(mfrow = c(1,1))
for(day in 1:365){
  plot(1:12*2, Deltas[,2,day], type = "b", lwd = 2, col = color[2],
       xlab = "time of day", ylab = "Pr(state 1)",
       ylim = c(0,1), bty = "n",
       main = paste("doy =", day))
  Sys.sleep(0.1)
}

par(mfrow = c(3,4))
for(day in c(1:12*30)){
  plot(1:12*2, Deltas[,2,day], type = "b", lwd = 2, col = color[2],
       xlab = "time of day", ylab = "Pr(state 1)",
       ylim = c(0,1), bty = "n",
       main = paste("doy =", day))
}




# Space-time interaction --------------------------------------------------

## penalised likelihood function
pnll3 = function(par) {
  getAll(par, dat)

  # fixpar
  sigma = exp(logsigma); REPORT(sigma)
  kappa = exp(logkappa); REPORT(kappa)

  # state process
  Gamma = tpm(beta0)
  delta = stationary(Gamma)

  # state-dependent process
  beta_mu = cbind(beta0_mu, beta_mu_xy, beta_mu_t, beta_mu_ti); REPORT(beta_mu)
  Mu = exp(Z %*% t(beta_mu))

  allprobs = matrix(1, nrow = length(step), ncol = N)
  ind = which(!is.na(step) & !is.na(angle))
  for(j in 1:N){
    allprobs[ind,j] = dgamma2(step[ind], Mu[ind,j], sigma[j]) *
      dvm(angle[ind], 0, kappa[j])
  }

  -forward(delta, Gamma, allprobs) +
    penalty2(par[random], S, lambda) # computes quadratic penalty
}

# modmat_hid = make_matrices(~ s(tod, bs = "cc"),
#                            data = data,
#                            knots = list(tod = c(0, 24),
#                                         doy = c(1, 366))) # telling mgcv where to wrap the basis
# Z_hid = modmat_hid$Z
# S_hid = modmat_hid$S

data$t = 1:nrow(data)
modmat = make_matrices(~ s(x, y, k = 30) +
                         s(t, bs = "cr", k = 10)+
                         ti(x, y, t, bs = c("tp", "cr"), k = c(30, 10), d = c(2,1)),
                       data = data)
Z = modmat$Z
S = modmat$S
pardim = modmat$pardim


# initial parameter values
N = 2
par = list(
  beta0 = c(-2,-2),               # state process intercept
  beta0_mu = log(c(0.35, 1.1)),   # state-dependent step mean intercept
  logsigma = log(c(0.25, 0.75)),  # state-dependent step sd
  logkappa = log(c(0.2, 0.7)),    # state-dependent angle concentration
  beta_mu_xy = matrix(0, N*(N-1), pardim$`s(x,y)`),
  beta_mu_t = matrix(0, N*(N-1), pardim$`s(t)`),
  beta_mu_ti = matrix(0, N*(N-1), pardim$`ti(x,y,t)`)
  )

random = c("beta_mu_xy", "beta_mu_t", "beta_mu_ti")

dat = list(step = data$step, angle = data$angle, N = 2,
           Z = Z, S = S, random = random,
           lambda = c(rep(1e4, 4), rep(1e6, 4)))

map = list(lambda = c(1,1,
                      2,2,
                      3,4,3,4))

system.time(
  mod3 <- qreml2(pnll3, par, dat,
                 random = random ,
                 silent = 0,
                 map = map)
)

par = mod3$par
beta_mu_xy = matrix(0, 2, ncol(mod3$beta_mu))
beta_mu_xy[,1+1:pardim$`s(x,y)`] = par$beta_mu_xy

# plot the estimated field (spatial component)
# Define grid range
x_seq <- seq(min(data$x), max(data$x), length.out = 50)  # 100 grid points in x
y_seq <- seq(min(data$y), max(data$y), length.out = 50)  # 100 grid points in y
grid_data <- expand.grid(x = x_seq, y = y_seq, t = 1)
Z_pred = pred_matrix(modmat, newdata = grid_data)
Mu_pred = exp(Z_pred %*% t(beta_mu_xy))

# Convert predictions to matrices for plotting
Z1_matrix <- matrix(Mu_pred[, 1], nrow = length(x_seq), ncol = length(y_seq))
Z2_matrix <- matrix(Mu_pred[, 2], nrow = length(x_seq), ncol = length(y_seq))

par(mfrow = c(1,2))
image(x_seq, y_seq, Z1_matrix, col = hcl.colors(100),
      xlab = "X", ylab = "Y")
image(x_seq, y_seq, Z2_matrix, col = hcl.colors(100),
      xlab = "X", ylab = "Y")


# plot the estimated field
# Define grid range
x_seq <- seq(min(data$x), max(data$x), length.out = 50)  # 100 grid points in x
y_seq <- seq(min(data$y), max(data$y), length.out = 50)  # 100 grid points in y
t_seq <- seq(min(data$numID), max(data$numID), length.out = 50)  # 50 time points

# Create a data frame with all combinations
grid_data <- expand.grid(x = x_seq, y = y_seq, numID = t_seq)

Z_pred = pred_matrix(modmat_obs, newdata = grid_data)
Mu_pred = exp(Z_pred %*% t(cbind(par$beta0_mu, par$beta_mu)))

times = unique(grid_data$numID)

par(mfrow = c(1,2))
for (i in seq_along(times)) {
  # Extract predicted values for each field
  Mu_pred1 <- Mu_pred[grid_data$numID == times[i], 1]
  Mu_pred2 <- Mu_pred[grid_data$numID == times[i], 2]

  # Convert predictions to matrices for plotting
  Z1_matrix <- matrix(Mu_pred1, nrow = length(x_seq), ncol = length(y_seq))
  Z2_matrix <- matrix(Mu_pred2, nrow = length(x_seq), ncol = length(y_seq))

  image(x_seq, y_seq, Z1_matrix, col = hcl.colors(100),
        xlab = "X", ylab = "Y", main = paste("Heatmap of Field 1 at Time =", times[i]))
  image(x_seq, y_seq, Z2_matrix, col = hcl.colors(100),
        xlab = "X", ylab = "Y", main = paste("Heatmap of Field 1 at Time =", times[i]))

  Sys.sleep(0.2)  # Pause to create animation effect
}





par(mfrow = c(1,2))
for(i in 1:length(times)){
  # Extract predicted values for each field
  Mu_pred1 <- Mu_pred[grid_data$numID == times[i], 1]
  Mu_pred2 <- Mu_pred[grid_data$numID == times[i], 2]

  # Convert predictions to matrices for plotting
  Z1_matrix <- matrix(Mu_pred1, nrow = length(x_seq), ncol = length(y_seq))
  Z2_matrix <- matrix(Mu_pred2, nrow = length(x_seq), ncol = length(y_seq))

  # Plot field 1
  contour(x_seq, y_seq, Z1_matrix, col = terrain.colors(10),
          xlab = "X", ylab = "Y", main = paste("Estimated Field 1 at Time =", time_fixed))

  # Plot field 2
  contour(x_seq, y_seq, Z2_matrix, col = terrain.colors(10),
          xlab = "X", ylab = "Y", main = paste("Estimated Field 2 at Time =", time_fixed))

  Sys.sleep(0.25)
}

library(plot3D)
par(mfrow = c(1,1))
for (i in seq_along(times)) {
  # Extract predicted values for each field
  Mu_pred1 <- Mu_pred[grid_data$numID == times[i], 1]
  Mu_pred2 <- Mu_pred[grid_data$numID == times[i], 2]

  # Convert predictions to matrices for plotting
  Z1_matrix <- matrix(Mu_pred1, nrow = length(x_seq), ncol = length(y_seq))
  Z2_matrix <- matrix(Mu_pred2, nrow = length(x_seq), ncol = length(y_seq))

  persp3D(x = x_seq, y = y_seq, z = Z1_matrix,
          colvar = Z1_matrix, col = terrain.colors(100),
          theta = 135, phi = 30, shade = 0.5,
          xlab = "X", ylab = "Y", zlab = "Field 1 Prediction",
          main = paste("Time =", times[i]))

  Sys.sleep(0.25)  # Pause to create animation effect
}




