## packages
# devtools::install_github("janoleko/LaMa")
library(LaMa)
library(scales)


## data
data <- readRDS("./data/fruitflies_genotype0.rds")
data <- readRDS("./data/fruitflies_genotype1.rds")


## unique animal IDs
aniIDs <- unique(data$aniID)
nAnimals <- length(aniIDs)


## color vector
color <- c("orange", "deepskyblue")


## penalised log-likelihood function
pnll <- function(par) {
  getAll(par, dat) # extract everything from parameter and data list

  # state-dependent process parameters (mean and dispersion)
  mu <- exp(logmu); REPORT(mu)
  phi <- exp(logphi); REPORT(phi)

  # state process coefficient matrix
  beta <- cbind(
    beta0, # intercept LD, DD joint
    beta_riLD, beta_riDD, # random intercepts LD, DD
    beta_timeLD, beta_timeDD, # main effect time of day LD, DD
    beta_tiLD, beta_tiDD # interaction time of day and individual LD, DD
  )

  # computing the transition probability matrix via inverse mlogit
  Gamma <- tpm_g(Z, beta)

  # periodically stationary initial distribution for each track
  nTracks <- length(trackInd)
  Delta <- matrix(0, nTracks, 2)
  for (tr in 1:nTracks) Delta[tr, ] <- stationary_p(Gamma[, , trackInd[tr] + 0:47], t = 1)

  # state-dependent probabilities
  allprobs <- matrix(1, length(activity), 2)
  ind <- which(!is.na(activity))
  size <- 1 / phi # reparametrisation of negbinom
  prob <- size / (size + mu) # reparametrisation of negbinom
  allprobs[ind, ] <- cbind(
    dnbinom(activity[ind], prob[1], size = size[1]),
    dnbinom(activity[ind], prob[2], size = size[2])
  )

  # forward algorithm + penalty()
  -forward_g(Delta, Gamma, allprobs, trackID) +
    penalty2(list(
      beta_riLD, beta_riDD,
      beta_timeLD, beta_timeDD,
      beta_tiLD, beta_tiDD), S, lambda)
}


## build design and penalty matrices for the model
modmat <- make_matrices(
  ~ condition +
    s(aniID, bs = "re", by = condition) +
    s(tod, bs = "cc", by = condition, k = 8) +
    ti(aniID, tod, bs = c("re", "cc"), by = condition, k = c(nAnimals, 8)),
  data = data,
  knots = list(tod = c(0, 24))
)

Z <- modmat$Z # design matrix
S <- modmat$S # list of penalty matrices for s(ID) and s(time) and list of penalty matrices for ti(ID, time)
pardim <- modmat$pardim

# initial parameters
par <- list(
  logmu = log(c(4, 55)),                                        # state-dependent mean
  logphi = log(c(10, 0.5)),                                     # state-dependent dispersion
  beta0 = matrix(-2, 2, 2),                                     # state process intercepts
  beta_riLD = matrix(0, 2, pardim$`s(aniID):conditionLD`),      # state process random intercept DD
  beta_riDD = matrix(0, 2, pardim$`s(aniID):conditionDD`),      # state process random intercept DD
  beta_timeLD = matrix(0, 2, pardim$`s(tod):conditionLD`),      # state process smooth time LD
  beta_timeDD = matrix(0, 2, pardim$`s(tod):conditionDD`),      # state process smooth time DD
  beta_tiLD = matrix(0, 2, pardim$`ti(aniID,tod):conditionLD`), # state process interaction LD
  beta_tiDD = matrix(0, 2, pardim$`ti(aniID,tod):conditionDD`)  # state process interaction DD
)

# data and hyperparameter list
dat <- list(
  activity = data$activity,
  condition = data$condition,
  trackID = data$trackID,
  trackInd = calc_trackInd(as.vector(data$trackID)),
  Z = Z,
  S = S,
  lambda = c(rep(1e3, 8), rep(1e4, 8)) # initial lambda: length equals total number of random effects
)

# mapping lambdas: estimate one for each smooth, but the same for both off-diagonal tpm entries
map <- list(lambda = c(
  1, 1,       # only one smoothness parameter for random intercept LD
  2, 2,       # only one smoothness parameter for random intercept DD
  3, 3,       # only one smoothness parameter for time smooth LD
  4, 4,       # only one smoothness parameter for time smooth DD
  5, 6, 5, 6, # anisotropic smoothing but same of tpm entries ti LD
  7, 8, 7, 8  # anisotropic smoothing but same of tpm entries ti DD
))


## fitting the model with qREML
system.time(
  mod <- qreml2(pnll, par, dat,
    random = c(
      "beta_riLD", "beta_riDD",
      "beta_timeLD", "beta_timeDD",
      "beta_tiLD", "beta_tiDD"
    ),
    map = map,
    silent = 0,
    alpha = 0.4,
    control = list(REPORT = 25),
    maxiter = 50
  )
)
# mod0 2.75 h
# mod1 ?
# save(mod, file = "./case_studies/objects/fruitflies_mod_genotype1.RData")

load(file = "./case_studies/objects/fruitflies_mod_genotype1.RData")

round(sqrt(mod$lambda^-1), 2)

# sdr = sdreport_outer(mod)

# save(mod, file = "./case_studies/models/fruitflies_mod_full.RData")
# load(file = "./case_studies/models/fruitflies_mod_full.RData")

AIC(mod)
BIC(mod)

## extracting parameters
beta <- mod$beta

## plotting periodically stationary distribution
L <- 48 # number of unique times of day
Delta_mean <- array(dim = c(L, 2, 2)) # mean array
Delta <- array(dim = c(L, 2, 2, nAnimals)) # array for all flies
conds <- unique(data$condition) # unique conditions

# estimated stationary distribution for each fly
for (cond in 1:2) {
  Z_pred <- pred_matrix(
    modmat,
    data.frame(
      tod = rep(1:L / 2, nAnimals),
      aniID = rep(aniIDs, each = L),
      condition = conds[cond]
    )
  )
  Z_list <- lapply(seq_len(nAnimals), function(a) Z_pred[((a - 1) * L + 1):(a * L), , drop = FALSE])
  for (a in 1:nAnimals) {
    Delta[, , cond, a] <- stationary_p(tpm_g(Z_list[[a]], beta))
  }
}
# save(Delta, file = "./case_studies/objects/fruitflies_Delta_genotype0.RData")
load(file = "./case_studies/objects/fruitflies_Delta_genotype0.RData")

for (cond in 1:2) {
  Z_pred <- pred_matrix(
    modmat,
    data.frame(
      tod = 1:L / 2,
      aniID = "newID",
      condition = conds[cond]
    )
  )
  Delta_mean[, , cond] <- stationary_p(tpm_g(Z_pred, beta))
}
# save(Delta_mean, file = "./objects/case_studies/fruitflies_Delta_mean_genotype0.RData")
load(file = "./case_studies/objects/fruitflies_Delta_mean_genotype0.RData")

# plotting
conditions <- c("LD", "DD")
state <- 2

color <- c("orange", "deepskyblue")
# pdf("./case_studies/figs/fruitflies_stationary_genotype0.pdf", width = 8, height = 3.5)
par(mfrow = c(1, 4))
for (cond in 1:2) {
  plot(NA,
    xlim = c(0, 24), ylim = c(0, 1), xlab = "time of day", ylab = "Pr(active)",
    bty = "n", main = conditions[cond], xaxt = "n"
  )
  axis(1, at = seq(0, 24, by = 4), labels = seq(0, 24, by = 4))
  # individual flies
  for (a in 1:nAnimals) {
    lines(1:L / 2, Delta[, state, cond, a], col = alpha(a, 0.3), type = "l", pch = 20)
  }
  # mean
  lines(1:L / 2, Delta_mean[, state, cond], type = "l", lwd = 2)
  points(1:L / 2, Delta_mean[, state, cond], pch = 16)

  # indicate light
  polygon(x = c(0, 8, 8, 0), y = c(-0.02, -0.02, 0, 0), col = "black", border = "black")
  if(cond == 1){
    polygon(x = c(8, 20, 20, 8), y = c(-0.02, -0.02, 0, 0), col = "white", border = "black")
  } else{
    polygon(x = c(8, 20, 20, 8), y = c(-0.02, -0.02, 0, 0), col = "#00000098", border = "black")

  }
  polygon(x = c(20, 24, 24, 20), y = c(-0.02, -0.02, 0, 0), col = "black", border = "black")
}


# differences only
DeltaDiff <- Delta
for (a in 1:nAnimals) {
  DeltaDiff[, , , a] <- Delta[, , , a] - Delta_mean
}

# pdf("./case_studies/figs/fruitflies_differences2.pdf", width = 8, height = 5)
# par(mfrow = c(1, 2))
for (cond in 1:2) {
  plot(NA,
    xlim = c(0, 24), ylim = c(-0.7, 0.7), xlab = "time of day", ylab = "difference to main effect Pr(active)",
    bty = "n", main = conditions[cond], xaxt = "n"
  )
  axis(1, at = seq(0, 24, by = 4), labels = seq(0, 24, by = 4))
  # mean
  abline(h = 0, lwd = 2)
  # individual flies
  for (a in 1:nAnimals) {
    lines(1:L / 2, DeltaDiff[, state, cond, a], col = alpha(a, 0.3), type = "l", pch = 20)
  }
  # indicate light
  polygon(x = c(0, 8, 8, 0), y = c(-0.02, -0.02, 0, 0), col = "black", border = "black")
  if(cond == 1){
    polygon(x = c(8, 20, 20, 8), y = c(-0.02, -0.02, 0, 0), col = "white", border = "black")
  } else{
    polygon(x = c(8, 20, 20, 8), y = c(-0.02, -0.02, 0, 0), col = "#00000098", border = "black")

  }
  polygon(x = c(20, 24, 24, 20), y = c(-0.02, -0.02, 0, 0), col = "black", border = "black")
}
#dev.off()



