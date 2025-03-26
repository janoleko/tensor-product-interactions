## Kruger national park buffalo case study
# https://www.movebank.org/cms/webapp?gwt_fragment=page=studies,path=study1764627

## packages
library(mvtnorm)
# devtools::install_github("janoleko/LaMa") # development version
library(LaMa)
# install.packages("scales") # for muted colors
library(scales)


## colors
color = c("orange", "deepskyblue", "seagreen2")


## data
data = read.csv("./data/african_buffalo.csv")
head(data)
data = subset(data, aniID == "Queen")



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
      dvm(angle[ind], c(pi, 0, 0)[j], kappa[j])
  }
  ## forward algorithm
  -forward(delta, Gamma, allprobs)
}

par = list(logmu = log(c(0.02, 0.3, 1)),
           logsigma = log(c(0.02, 0.3, 1)),
           logkappa = log(c(0.1, 0.5, 1)),
           beta0 = rep(-3, 6))
N = 3

dat = list(step = data$step, angle = data$angle, N = N)

obj = MakeADFun(nll, par)
opt = nlminb(obj$par, obj$fn, obj$gr)

mod = obj$report()
mu = mod$mu
sigma = mod$sigma
kappa = mod$kappa
delta = mod$delta


hist(data$step, xlim = c(0,5), breaks = 200, prob = TRUE)
for(j in 1:N){
  curve(delta[j] * dgamma2(x, mu[j], sigma[j]), add = TRUE, n = 500, lwd = 2, col = color[j])
}
hist(data$angle, breaks = 40, prob = TRUE)
for(j in 1:N){
  curve(delta[j] * dvm(x, c(pi, 0, 0)[j], kappa[j]), add = TRUE, n = 500, lwd = 2, col = color[j])
}

states = viterbi(mod = mod)

plot(data$step[1:200], col = color[states[1:200]], type = "h")



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
  field = Z %*% beta
  Mu = exp(cbind(field, field + logmu[1], field + logmu[2]))

  allprobs = matrix(1, nrow = length(step), ncol = N)
  ind = which(!is.na(step) & !is.na(angle))
  for(j in 1:N){
    allprobs[ind,j] = dgamma2(step[ind], Mu[ind,j], sigma[j]) *
      dvm(angle[ind], c(pi, 0, 0)[j], kappa[j])
  }

  -forward(delta, Gamma, allprobs) +
    penalty2(list(beta_field), S, lambda) # computes quadratic penalty
}


modmat = make_matrices(~ s(x, y, k = 100), data = data)
Z = modmat$Z
S = modmat$S
pardim = modmat$pardim


# initial parameter values
N = 3
par = list(
  beta0 = rep(-2, N*(N-1)),                  # state process intercept
  logsigma = log(c(0.02, 0.3, 1)),           # state-dependent step sd
  logkappa = log(c(0.1, 0.5, 1)),            # state-dependent angle concentration
  beta0_field = log(0.02),
  beta_field = rep(0, pardim$`s(x,y)`),
  logmu = c(2,4)
)

dat = list(step = data$step, angle = data$angle, N = N,
           Z = Z, S = S,
           lambda = 1e4)

system.time(
  mod2 <- qreml2(pnll, par, dat,
                 random = "beta_field",
                 silent = 0)
)


beta = mod2$beta
logmu = mod2$par$logmu

# plot the estimated field (spatial component)
# Define grid range
x_seq <- seq(min(data$x), max(data$x), length.out = 50)  # 100 grid points in x
y_seq <- seq(min(data$y), max(data$y), length.out = 50)  # 100 grid points in y
grid_data <- expand.grid(x = x_seq, y = y_seq)
Z_pred = pred_matrix(modmat, newdata = grid_data)
field = exp(Z_pred %*% beta)

# Convert predictions to matrices for plotting
Z_matrix = matrix(field, nrow = length(x_seq), ncol = length(y_seq))

states = viterbi(mod = mod2)


par(mfrow = c(1,1))
image(x_seq, y_seq, Z_matrix, col = hcl.colors(100),
      xlab = "long", ylab = "lat", main = paste("state",j))
points(data$x, data$y, col = "white")


