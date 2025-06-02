## packages
# devtools::install_github("janoleko/LaMa")
library(LaMa)
library(scales)


## loading the data
# wild-type
data0 <- readRDS("./data/fruitflies_genotype0.rds")
# modified genotype
data1 <- readRDS("./data/fruitflies_genotype1.rds")

data = list(data0 = data0, data1 = data1)

## number of observations
nrow(data0)
nrow(data1)


## ranges
range(data0$activity, na.rm = TRUE)
range(data1$activity, na.rm = TRUE)


## unique animal IDs
aniIDs <- list(
  unique(data0$aniID),
  unique(data1$aniID)
  )
(nAnimals <- sapply(aniIDs, length))


## color vector for plotting later
color <- c("orange", "deepskyblue")


## penalised log-likelihood function for the tensor-product model
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

# model matrices for both genotypes
modmat0 <- make_matrices(
  ~ condition +
    s(aniID, bs = "re", by = condition) +
    s(tod, bs = "cc", by = condition, k = 10) +
    ti(aniID, tod, bs = c("re", "cc"), by = condition, k = c(nAnimals[1], 10)),
  data = data0,
  knots = list(tod = c(0, 24))
)
modmat1 <- make_matrices(
  ~ condition +
    s(aniID, bs = "re", by = condition) +
    s(tod, bs = "cc", by = condition, k = 10) +
    ti(aniID, tod, bs = c("re", "cc"), by = condition, k = c(nAnimals[2], 10)),
  data = data1,
  knots = list(tod = c(0, 24))
)
modmat <- list(modmat0, modmat1)


## fitting models
# change this to fit the model to the wild type (1) or modified genotype (2)
genotype = 2

Z <- modmat[[genotype]]$Z # design matrix
S <- modmat[[genotype]]$S # list of penalty matrices for s(ID) and s(time) and list of penalty matrices for ti(ID, time)
pardim <- modmat[[genotype]]$pardim # to easily set up initial parameters

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
thisdata <- data[[genotype]]
dat <- list(
  activity = thisdata$activity, # observations
  condition = thisdata$condition, # LD or DD variable
  trackID = thisdata$trackID, # track ID
  trackInd = calc_trackInd(thisdata$trackID), # index of the first observation of each track (for initial distributions)
  Z = Z, # design matrix for smooth terms
  S = S, # penalty matrices for smooth terms
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


# ## fitting the model with qREML
# system.time(
#   mod <- qreml(pnll, par, dat,
#     random = c(
#       "beta_riLD", "beta_riDD",
#       "beta_timeLD", "beta_timeDD",
#       "beta_tiLD", "beta_tiDD"
#     ),
#     map = map,
#     silent = 0
#   )
# )
# # mod0 2.5 h
# # mod1 2.3 h
# saveRDS(mod0, "./case_studies/objects/fruitflies_mod_genotype0.rds")
# saveRDS(mod1, "./case_studies/objects/fruitflies_mod_genotype1.rds")

########################################################################


## loading saved models
mod0 <- readRDS(file = "./case_studies/objects/fruitflies_mod_genotype0.rds")
mod1 <- readRDS(file = "./case_studies/objects/fruitflies_mod_genotype1.rds")

# # uncertainty quantification in estimated variance paramters
# sdr0 = sdreport_outer(mod0, invert = TRUE)
# saveRDS(sdr0, "./case_studies/objects/fruitflies_sdr_genotype0.rds")

# sdr1 = sdreport_outer(mod1, invert = TRUE)
# saveRDS(sdr1, "./case_studies/objects/fruitflies_sdr_genotype1.rds")

# information criteria
AIC(mod0)
AIC(mod1)
BIC(mod0)
BIC(mod1)


## plotting state-dependent distributions
states0 = viterbi_g(mod = mod0)
delta_hat0 = c(mean(states0 == 1), mean(states0 == 2))
states1 = viterbi_g(mod = mod1)
delta_hat1 = c(mean(states1 == 1), mean(states1 == 2))

# pdf("./case_studies/figs/fruitflies_marginal.pdf", width = 8, height = 4)

par(mfrow = c(1,2))
hist(data[[1]]$activity, breaks = 30, bor = "white",
     xlim = c(0, 150), ylim = c(0, 0.03), prob = TRUE,
     main = "wild type", xlab = "activity count", ylab = "frequency")
size = 1/mod0$phi
mu = mod0$mu
pmf = matrix(NA, 2, 151)
for(j in 1:2){
  x = 0:150
  pmf[j,] = delta_hat0[j] * dnbinom(x, mu = mu[j], size = size[j])
  lines(x, pmf[j,], type = "l", lwd = 2, col = color[j])
}
lines(x, colSums(pmf), type = "l", lwd = 2, col = "black", lty = 2)
legend("topright", legend = c("inactive", "active", "marginal"), col = c(color[1:2], "black"),
       lty = c(1,1,2), bty = "n", lwd = 2)

hist(data[[2]]$activity, breaks = 30, bor = "white",
     xlim = c(0, 150), ylim = c(0, 0.03), prob = TRUE,
     main = "modified genotype", xlab = "activity count", ylab = "frequency")
size = 1/mod1$phi
mu = mod1$mu
for(j in 1:2){
  x = 0:150
  pmf[j,] = delta_hat0[j] * dnbinom(x, mu = mu[j], size = size[j])
  lines(x, pmf[j,], type = "l", lwd = 2, col = color[j])
}
lines(x, colSums(pmf), type = "l", lwd = 2, col = "black", lty = 2)

# dev.off()



## extracting parameters
beta <- list(beta0 = mod0$beta,
             beta1 = mod1$beta)

## plotting periodically stationary distribution
L <- 48 # number of unique times of day
Delta_mean <- array(dim = c(L, 2, 2, 2)) # mean array
Delta <- list(array(dim = c(L, 2, 2, nAnimals[1])),
              array(dim = c(L, 2, 2, nAnimals[2])))
conds <- unique(data0$condition) # unique conditions

# estimated stationary distribution for each fly
for(m in 1:2){
  for (cond in 1:2) {
    Z_pred <- pred_matrix(
      modmat[[m]],
      data.frame(
        tod = rep(1:L / 2, nAnimals[m]),
        aniID = rep(aniIDs[[m]], each = L),
        condition = conds[cond]
      )
    )
    Z_list <- lapply(seq_len(nAnimals[m]), function(a) Z_pred[((a - 1) * L + 1):(a * L), , drop = FALSE])
    for (a in seq_len(nAnimals[m])) {
      Delta[[m]][, , cond, a] <- stationary_p(tpm_g(Z_list[[a]], beta[[m]]))
    }
  }
}

for(m in 1:2){
  for (cond in 1:2) {
    Z_pred <- pred_matrix(
      modmat[[m]],
      data.frame(
        tod = 1:L / 2,
        aniID = "newID",
        condition = conds[cond]
      )
    )
    Delta_mean[, , cond, m] <- stationary_p(tpm_g(Z_pred, beta[[m]]))
  }
}


# plotting
conditions <- c("LD", "DD")
state <- 2

color <- c("orange", "deepskyblue")
dark_color <- adjustcolor(color, red.f = 0.8, green.f = 0.8, blue.f = 0.8)
# darken colors


# pdf("./case_studies/figs/fruitflies_stationary.pdf", width = 8, height = 3.5)
par(mfrow = c(1, 4))

for(m in 1:2){
  for (cond in 1:2) {
    plot(NA,
         xlim = c(0, 24), ylim = c(0, 1), xlab = "time of day", ylab = "Pr(active)",
         bty = "n", main = conditions[cond], xaxt = "n"
    )

    # indicate light
    polygon(x = c(0, 8, 8, 0), y = c(-0.005, -0.005, 0.005, 0.005), col = "black", border = "black")
    if(cond == 1){
      polygon(x = c(8, 20, 20, 8), y = c(-0.005, -0.005, 0.005, 0.005), col = "white", border = "black")
    } else{
      polygon(x = c(8, 20, 20, 8), y = c(-0.005, -0.005, 0.005, 0.005), col = "#00000098", border = "black")

    }
    polygon(x = c(20, 24, 24, 20), y = c(-0.005, -0.005, 0.005, 0.005), col = "black", border = "black")

    axis(1, at = seq(0, 24, by = 4), labels = seq(0, 24, by = 4))
    # individual flies
    for (a in 1:nAnimals[m]) {
      lines(1:L / 2, Delta[[m]][, state, cond, a], col = alpha(color[m], 0.3), type = "l", pch = 20)
    }
    # mean
    lines(1:L / 2, Delta_mean[, state, cond, m], type = "l", lwd = 2, col = "#000000")
    points(1:L / 2, Delta_mean[, state, cond, m], pch = 16, col = "#000000")
  }
}
# dev.off()
