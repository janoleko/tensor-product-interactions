## African lion central Kalahari Botswana
# https://www.movebank.org/cms/webapp?gwt_fragment=page=studies,path=study3791354435

## packages
library(mvtnorm)
# devtools::install_github("janoleko/LaMa") # development version
library(LaMa)
# install.packages("scales") # for muted colors
library(scales)


## colors
color = c("orange", "deepskyblue", "seagreen2")


## data
data = read.csv("./data/lion1.csv")
head(data)
nrow(data)

plot(data$x, data$y)



# First simple HMM to understand patterns ---------------------------------

nll = function(par){
  getAll(par, dat)
  ## transforming parameters
  mu = exp(logmu); REPORT(mu)
  sigma = exp(logsigma); REPORT(sigma)
  kappa = exp(logkappa); REPORT(kappa)
  Gamma = tpm(beta0)
  delta = stationary(Gamma)
  ## computing state-dependent densities
  allprobs = matrix(1, nrow = length(step), ncol = N)
  ind = which(!is.na(step) & !is.na(angle))
  for(j in 1:N){
    allprobs[ind,j] = dgamma2(step[ind], mu[j], sigma[j]) *
      dvm(angle[ind], c(pi, 0)[j], kappa[j])
  }
  ## forward algorithm
  -forward(delta, Gamma, allprobs)
}


hist(data$step, breaks = 200, xlim = c(0,4), ylim = c(0, 2), prob = TRUE)
hist(data$angle, breaks = 50, prob = TRUE)

par = list(logmu = log(c(0.01, 1)),
           logsigma = log(c(0.01, 0.5)),
           logkappa = log(c(0.5, 5)),
           beta0 = rep(-2, 2))
N = 2

dat = list(step = data$step, angle = data$angle, N = N)

obj = MakeADFun(nll, par)
opt = nlminb(obj$par, obj$fn, obj$gr)

mod = obj$report()
(mu = mod$mu)
(sigma = mod$sigma)
(kappa = mod$kappa)
Gamma = mod$Gamma
round(Gamma, 3)
delta = mod$delta


hist(data$step, xlim = c(0, 4), ylim = c(0, 2), breaks = 200, prob = TRUE)
for(j in 1:N){
  curve(delta[j] * dgamma2(x, mu[j], sigma[j]), add = TRUE, n = 500, lwd = 2, col = color[j])
}
curve(delta[1] * dgamma2(x, mu[1], sigma[1]) +
        delta[2] * dgamma2(x, mu[2], sigma[2]), add = TRUE, n = 500, lwd = 2, lty = 2)

hist(data$angle, breaks = 40, prob = TRUE)
for(j in 1:N){
  curve(delta[j] * dvm(x, c(pi, 0, 0)[j], kappa[j]), add = TRUE, n = 500, lwd = 2, col = color[j])
}
curve(delta[1] * dvm(x, pi, kappa[1]) +
        delta[2] * dvm(x, 0, kappa[2]) +
        delta[3] * dvm(x, 0, kappa[3]), add = TRUE, n = 500, lwd = 2, lty = 2)

states = viterbi(mod = mod)

plot(data$step[1:500], col = color[states[1:500]], type = "h")




# Model with Spatial effect on the mean step length -----------------------

## penalised likelihood function
pnll = function(par) {
  getAll(par, dat)

  # fixpar
  sigma = exp(logsigma); REPORT(sigma)
  kappa = exp(logkappa); REPORT(kappa)

  # state process
  Gamma = tpm(beta0)
  delta = stationary(Gamma)

  # state-dependent process
  beta = c(beta0_field, beta_field); REPORT(beta)
  logmu2 = Z %*% beta
  Mu = exp(cbind(logmu1, logmu2))

  allprobs = matrix(1, nrow = length(step), ncol = N)
  ind = which(!is.na(step) & !is.na(angle))
  for(j in 1:N){
    allprobs[ind,j] = dgamma2(step[ind], Mu[ind,j], sigma[j]) *
      dvm(angle[ind], c(pi, 0)[j], kappa[j])
  }

  -forward(delta, Gamma, allprobs) +
    penalty2(list(beta_field), S, lambda) # computes quadratic penalty
}


modmat = make_matrices(~ s(x, y, k = 50), data = data)
Z = modmat$Z
S = modmat$S
pardim = modmat$pardim


# initial parameter values
N = 2
par = list(
  beta0 = rep(-2, N*(N-1)),                  # state process intercept
  logsigma = log(c(0.01, 0.5)),           # state-dependent step sd
  logkappa = log(c(0.5, 5)),            # state-dependent angle concentration
  beta0_field = log(0.55),
  beta_field = rep(0, pardim$`s(x,y)`),
  logmu1 = log(0.01)
)

dat = list(step = data$step, angle = data$angle, N = N,
           Z = Z, S = S,
           lambda = 1e3)

system.time(
  mod2 <- qreml2(pnll, par, dat,
                 random = "beta_field",
                 silent = 0)
)

states = viterbi(mod = mod2)


beta = mod2$beta
ind = which(states == 2)
mu2 = exp(Z %*% beta)
mu2[-ind] = NA
# Define a color gradient
color_gradient <- col_numeric(palette = "viridis", domain = range(mu2, na.rm = TRUE))
gradient_colors <- color_gradient(mu2[ind])
color = rep(NA, nrow(data))
color[states == 1] = 1
color[states == 2] = gradient_colors

plot(data$x, data$y, col = color, pch = 20)
points(data$x[-ind], data$y[-ind], pch = 20)

plot(data$step, col = color)


# plot the estimated field (spatial component)
# Define grid range
x_seq <- seq(min(data$x), max(data$x), length.out = 50)  # 100 grid points in x
y_seq <- seq(min(data$y), max(data$y), length.out = 50)  # 100 grid points in y
grid_data <- expand.grid(x = x_seq, y = y_seq)
Z_pred = pred_matrix(modmat, newdata = grid_data)
mu2 = exp(Z_pred %*% beta)

# Convert predictions to matrices for plotting
Z_matrix = matrix(mu2, nrow = length(x_seq), ncol = length(y_seq))

par(mfrow = c(1,1))
image(x_seq, y_seq, Z_matrix, col = hcl.colors(100),
      xlab = "long", ylab = "lat", main = paste("state",j))
points(data$x, data$y, col = "#ffffff10", pch = 20)
