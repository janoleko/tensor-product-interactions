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

# number of observations
nrow(data)


## color vector for plotting later
color = c("orange", "deepskyblue")



# Initialising IC table ---------------------------------------------------

# information criteria
IC_table <- data.frame(model = 1:5,
                       llk = rep(NA, 5),
                       AIC = rep(NA, 5),
                       BIC = rep(NA, 5))


# Fitting homogeneous and parametric HMMs ---------------------------------

# likelihood function for homogeneous model
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

N = 2 # number of states

# initial parameter list
par = list(logmu = log(c(0.2, 2)),
           logsigma = log(c(0.2, 2)),
           logkappa = log(c(0.2, 1)),
           eta = rep(-2, 2))

# data and hyperparameter list
dat = list(step = data$step, angle = data$angle, N = 2)

# creating AD function
obj_hom = MakeADFun(nll_hom, par, silent = TRUE)

# fitting the model
opt_hom = nlminb(obj_hom$par, obj_hom$fn, obj_hom$gr)

# reporting for easy access
mod_hom = obj_hom$report()

mu = mod_hom$mu # state-dependent step means
sigma = mod_hom$sigma # state-dependent step sds
kappa = mod_hom$kappa # state-dependent angle concentration
delta = mod_hom$delta # stationary distribution

# decoding states
mod_hom$states = viterbi(mod_hom$delta, mod_hom$Gamma, mod_hom$allprobs)

# computing AIC and BIC
mod_hom$AIC = 2 * obj_hom$fn() + 2 * length(obj_hom$par)
mod_hom$BIC = 2 * obj_hom$fn() + log(nrow(data)) * length(obj_hom$par)

IC_table[1, 2:4] = c(-obj_hom$fn(), mod_hom$AIC, mod_hom$BIC) # saving


## plotting state-dependent distributions

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
# (one for both)
pnll_uni = function(par) {
  getAll(par, dat)

  # state-dependent parameters
  mu = exp(logmu); REPORT(mu)
  sigma = exp(logsigma); REPORT(sigma)
  kappa = exp(logkappa); REPORT(kappa)

  # state process
  Gamma = tpm_g(Z, cbind(beta0, beta_spline)) # tpm GAM
  delta = stationary_p(Gamma[,,1:L], t = 1) # periodically stationary initial dist

  # state-dependent probabilities
  allprobs = matrix(1, nrow = length(step), ncol = N)
  ind = which(!is.na(step) & !is.na(angle))
  for(j in 1:N){
    allprobs[ind,j] = dgamma2(step[ind], mu[j], sigma[j]) *
      dvm(angle[ind], 0, kappa[j])
  }

  # forward algorithm
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

dat = list(step = data$step, angle = data$angle, N = 2, L = 24,
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

dat = list(step = data$step, angle = data$angle, N = 2, L = 366,
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
  delta = stationary_p(Gamma[,,1:24], t = 1) # approx initial p-stationary

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
par = list(logmu = log(c(0.3, 1.1)),       # state-dependent step mean
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

summary(mod_ti)


## pseudo-residuals
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



sdr <- sdreport_outer(mod_ti)
sdr$report

AIC(mod_ti)
BIC(mod_ti)

beta <- mod_ti$beta

## Visualising the fitted model
n_plot = 200
tod_seq = seq(0,24, length = n_plot)
tod_fine = matrix(NA, n_plot, 12)
for(i in 1:n_plot){
  tod_fine[i,] = (tod_seq[i] + 1:12*2) %% 24
}
tod_fine = c(t(tod_fine))
Z_predict = pred_matrix(modmat_ti,
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

# ## sample parameters from posterior
# n_samples <- 2000
# set.seed(123)
# pars <- rmvnorm(n_samples, mod_ti$par_vec, solve(mod_ti$Hessian_conditional))
# pars <- apply(pars, 1, mod_ti$relist_par)
# betas <- lapply(pars, function(p) cbind(p$beta0, p$beta_tod, p$beta_doy, p$beta_ti))
#
# allDeltas = array(dim = c(n_plot, 2, 366, n_samples))
# for(day in 1:366){
#   print(day)
#   Z_day = Z_day_list[[day]]
#   for(i in 1:n_plot){
#     thisZ = Z_day[((i-1)*12+1):(i*12),]
#     for(s in 1:n_samples){
#       Gamma = tpm_g(thisZ, betas[[s]])
#       allDeltas[i,,day,s] = stationary_p(Gamma, t = 1)
#     }
#   }
# }
# Delta_q <- apply(allDeltas, 1:3, quantile, probs = c(0.025, 0.975))

# saveRDS(Delta_q, file = "./case_studies/objects/elephant_delta_quantiles.rds")

Delta_q = readRDS("./case_studies/objects/elephant_delta_quantiles.rds")

## plotting the estimated state distribution over the year

# pdf("./case_studies/figs/elephant_stationary.pdf", width = 7, height = 5)

par(mfrow = c(2,3))
source("./utils/gen_sun_colors.R")
sun_cycle_colors <- gen_sun_cycle_colors()
daynames = paste(c("Jan", "March", "May", "Jul", "September", "November"), 1)
days = c(1, 61, 121, 182, 243, 304)
for(day in days){
  plot(NA, xlim = c(0, 24), ylim = c(0,1), bty = "n",
       xlab = "time of day", ylab = "Pr(state exploratory)",
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

# dev.off()








### predict tpm over the year
Z_predict = predict(modmat_ti,
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
  Z_predict = predict(modmat_ti, newdata = data.frame(tod = tod, doy = doy))
  tpm_g(Z_predict, beta)[i,j,]
}
pred_eta = function(tod, doy, i){
  Z_predict = predict(modmat_ti, newdata = data.frame(tod = tod, doy = doy))
  Z_predict %*% beta[i,]
}

tod = seq(0, 24, length = 70)
doy = seq(0, 366, length = 70)

gamma_12 = outer(tod, doy, pred_gamma, i = 1, j = 2)
gamma_21 = outer(tod, doy, pred_gamma, i = 2, j = 1)

eta1 = outer(tod, doy, pred_eta, i = 1)
eta2 = outer(tod, doy, pred_eta, i = 2)


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
