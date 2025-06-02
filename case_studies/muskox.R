## packages

# model fitting
# devtools::install_github("janoleko/LaMa")
library(LaMa)

# spatial stuff
library(sf)
library(sp)
library(rnaturalearth)


## loading the data
data <- readRDS("./data/muskox_winter.rds")


## color palette
color <- c("orange", "deepskyblue", "seagreen2")


## EDA and spatial prep

# Extract the island boundaries
# Create an sf object from UTM data
utm_sf = st_as_sf(data[,c("x", "y")], coords = c("x", "y"),
                  crs = paste0("+proj=utm +zone=", 27, " +datum=WGS84"),
                  na.fail = FALSE)
latlon_sf = st_transform(utm_sf, crs = 4326)

# create a bounding box
bbox = st_bbox(latlon_sf)

# download map data for greenland only
world <- ne_countries(scale = "large", returnclass = "sf")

# filter for greenland
greenland = world[world$name == "Greenland",]
greenland_geom = st_geometry(greenland)

# restrict to bounding box
greenland_bbox = st_intersection(greenland_geom, st_as_sfc(bbox))

# Extract boundary points from the land polygon
boundary_coords <- st_coordinates(greenland_bbox)[, 1:2]

# Convert to a data frame
boundary_df <- data.frame(x = boundary_coords[,1], y = boundary_coords[,2])

# Ensure it's in UTM coordinates
boundary_df <- st_transform(st_as_sf(boundary_df, coords = c("x", "y"), crs = 4326), crs = "+proj=utm +zone=27 +datum=WGS84")
boundary_df <- st_coordinates(boundary_df)
boundary_df <- data.frame(x = boundary_df[,1], y = boundary_df[,2])

# split into two islands
bnd_main = boundary_df[25:142, ]
bnd_island = boundary_df[143:203, ]



###########################################################################
## Fitting models #########################################################
###########################################################################

# Initialising IC table ---------------------------------------------------

IC_table <- data.frame(model = 1:5,
                       llk = rep(NA, 5),
                       AIC = rep(NA, 5),
                       BIC = rep(NA, 5))

# Simple homogeneous model ------------------------------------------------

## negative log-likelihood function
nll <- function(par){
  getAll(par, dat) # make everything available without $
  ## transforming unconstrained parameters to working
  mu = exp(logmu); REPORT(mu)
  sigma = exp(logsigma); REPORT(sigma)
  kappa = exp(logkappa); REPORT(kappa)
  ## state process
  Gamma = tpm(beta) # inverse mlogit
  delta = stationary(Gamma) # stationary distribution
  ## computing state-dependent densities
  allprobs = matrix(1, nrow = length(step), ncol = N)
  ind = which(!is.na(step) & !is.na(angle))
  for(j in 1:N){
    allprobs[ind,j] = dgamma2(step[ind],mu[j],sigma[j]) *
      dvm(angle[ind], c(pi,0,0)[j], kappa[j])
  }
  ## forward algorithm
  -forward(delta, Gamma, allprobs, ID)
}

## initial parameters
par = list(
  logmu = log(c(5, 60, 400)),
  logsigma = log(c(5, 60, 400)),
  logkappa = log(c(0.01, 0.2, 1.2)),
  beta = rep(-2, 6)
)

## data
dat = list(
  step = data$step,
  angle = data$angle,
  N = 3,
  ID = data$ID
)

## creating AD function
obj <- MakeADFun(nll, par)

## model fitting
opt <- nlminb(obj$par, obj$fn, obj$gr)

## extracting results
mod_hom <- obj$report()
mu <- mod_hom$mu
sigma <- mod_hom$sigma
kappa <- mod_hom$kappa
delta <- mod_hom$delta

## information criteria

mod_hom$AIC = 2 * opt$objective + 2 * length(obj$par)
mod_hom$BIC = 2 * opt$objective + log(nrow(data)) * length(obj$par)
IC_table[1, 2:4] = c(- opt$objective, mod_hom$AIC, mod_hom$BIC)

## visualising fitted state-dependent distributions

# pdf("./case_studies/figs/muskox_marginal.pdf", width = 8, height = 4)

par(mfrow = c(1,2))
# step length
hist(data$step, breaks = 200, xlim = c(0, 600), ylim = c(0, 0.01), prob = TRUE, bor = "white",
     main = "", xlab = "step length", ylab = "density")
for(j in 1:3){
  curve(delta[j] * dgamma2(x, mu[j], sigma[j]), add = TRUE,
        col = color[j], lwd = 2, n = 500)
}
curve(delta[1] * dgamma2(x, mu[1], sigma[1]) +
        delta[2] * dgamma2(x, mu[2], sigma[2]) +
        delta[3] * dgamma2(x, mu[3], sigma[3]), add = TRUE,
      lty = 2, lwd = 2, n = 500)
legend("topright", legend = c("resting", "foraging", "travelling", "marginal"), col = c(color[1:3], "black"),
       lty = c(1,1,1,2), bty = "n", lwd = 2)
# turning angle
hist(data$angle, breaks = 30, prob = TRUE, bor = "white",
     main = "", xlab = "step length", ylab = "density")
for(j in 1:3){
  curve(delta[j] * dvm(x, c(pi,0,0)[j], kappa[j]), add = TRUE,
        col = color[j], lwd = 2, n = 500)
}
curve(delta[1] * dvm(x, pi, kappa[1]) +
        delta[2] * dvm(x, 0, kappa[2]) +
        delta[3] * dvm(x, 0, kappa[3]), add = TRUE,
      lty = 2, lwd = 2, n = 500)

# dev.off()

## state decoding
states <- viterbi(mod = mod_hom)

## plotting decoded states
# plot(data$x, data$y, col = scales::alpha(color[states], 0.3), pch = 20)
# plot(data$step, col = color[states], xlim = c(0, 10000), type = "h")



# Model with spatial smooth only ------------------------------------------

## creating model matrices
modmat <- make_matrices(~ s(x, y, k = 50), data = data) # calls mgcv under the hood
Z <- modmat$Z # design matrix
S <- modmat$S # list of penalty matrices (in this case only one)

# more flexible tpm constructor than in LaMa
# Z and beta and be lists, hence different covariates can be included for each
# off-diagonol element
tpm_g2 <- function(Z, beta, byrow = FALSE, report = TRUE) {

  K <- length(beta) # beta is list of parameter vectors
  N <- as.integer(0.5 + sqrt(0.25 + K), 0)
  if (report) RTMB::REPORT(beta)

  nObs <- max(sapply(Z, function(z){
    if(is.null(dim(z))){
      return(length(z))
    } else{
      nrow(z)
    }
  }))

  Eta <- AD(matrix(NaN, nObs, K))
  for(k in seq_len(K)){
    Eta2[,k] <- as.numeric(Z_list[[k]] %*% beta[[k]])
  }
  expEta <- exp(Eta)
  Gamma = AD(array(1, dim = c(N, N, nrow(expEta))))
  col_ind <- 1
  for (i in seq_len(N)) {
    for (j in seq_len(N)) {
      if (j != i) {
        if (byrow) {
          Gamma[i, j, ] <- expEta[, col_ind]
        }
        else {
          Gamma[j, i, ] <- expEta[, col_ind]
        }
        col_ind = col_ind + 1
      }
    }
  }
  for (i in seq_len(N)) {
    Gamma[i, , ] <- t(t(Gamma[i, , ])/rowSums(t(Gamma[i, , ])))
  }
  Gamma
}


## penalised negative log-likelihood function
pnll <- function(par){
  getAll(par, dat, warn = FALSE)
  ## transforming parameters
  mu <- exp(logmu); REPORT(mu)
  sigma <- exp(logsigma); REPORT(sigma)
  kappa <- exp(logkappa); REPORT(kappa)

  # state process
  # beta <- vector("list", 6)
  # beta[-4] <- eta
  # beta[[4]] <- c(beta0, beta_xy)
  # Z_list <- vector("list", 6)
  # Z_list[-4] <- rep(list(matrix(1, nrow(Z), 1)), 5)
  # Z_list[[4]] <- Z
  beta <- matrix(0, 6, ncol(Z))
  beta[-4,1] <- eta
  beta[4,] <- c(beta0, beta_xy)
  Gamma <- tpm_g(Z, beta)
  Delta <- stationary(Gamma[,,trackInd])

  ## computing state-dependent densities
  allprobs <- matrix(1, nrow = length(step), ncol = N)
  ind <- which(!is.na(step) & !is.na(angle))
  for(j in 1:N){
    allprobs[ind,j] <- dgamma2(step[ind],mu[j],sigma[j]) *
      dvm(angle[ind], c(pi,0,0)[j], kappa[j])
  }
  ## forward algorithm
  -forward_g(Delta, Gamma, allprobs, trackID = ID) +
    penalty2(list(beta_xy), S, lambda)
}

par <- list(
  logmu = log(c(5, 50, 300)),
  logsigma = log(c(5, 50, 300)),
  logkappa = log(c(0.01, 0.2, 1.2)),
  eta = rep(-2, 5),
  beta0 = -2,
  beta_xy = rep(0, modmat$pardim$`s(x,y)`)
)

dat <- list(
  step = data$step,
  angle = data$angle,
  ID = data$ID,
  trackInd = calc_trackInd(data$ID),
  N = 3,
  Z = Z,
  S = S,
  lambda = 1e3
)

mod_space = qreml(pnll, par, dat,
                  random = "beta_xy",
                  silent = 0)

saveRDS(mod_space, file = "./case_studies/objects/muskox_mod_space2.rds")

mod_space = readRDS("./case_studies/objects/muskox_mod_space.rds")

summary(mod_space)

## information criteria
IC_table[2, 2:4] = c(mod_space$llk, AIC(mod_space), BIC(mod_space))

eta = mod_space$par$eta
beta = c(mod_space$par$beta0, mod_space$par$beta_xy)

data$Pr23 = mod_space$Gamma[3,2,]

states = viterbi_g(mod = mod_space)
find = which(states == 2) # indices of decoded foraging state

# plot the estimated field
# Define grid range
x_seq <- seq(min(data$x), max(data$x), length.out = 150)  # 100 grid points in x
y_seq <- seq(min(data$y), max(data$y), length.out = 150)  # 100 grid points in y

# Create a data frame with all combinations
grid_data <- expand.grid(x = x_seq, y = y_seq)

hull_idx <- chull(data$x, data$y)  # Get convex hull indices
hull_coords <- cbind(data$x[hull_idx], data$y[hull_idx])  # Hull boundary points

# Check which grid points are inside the convex hull
inside_chull <- point.in.polygon(grid_data$x, grid_data$y, hull_coords[,1], hull_coords[,2])
inside_main = point.in.polygon(grid_data$x, grid_data$y, bnd_main$x, bnd_main$y)
inside_island = point.in.polygon(grid_data$x, grid_data$y, bnd_island$x, bnd_island$y)
inside = (inside_main | inside_island) & inside_chull

Z_pred = pred_matrix(modmat, newdata = grid_data)

beta <- matrix(0, 6, ncol(Z))
beta[-4,1] <- eta
beta[4,] <- c(mod_space$par$beta0, mod_space$par$beta_xy)
Gamma <- tpm_g(Z_pred, beta)

pr = Gamma[3,2,]
grid_data$pr = pr
# grid_data2 = grid_data
grid_data$pr[inside == 0] = NA
# grid_data2$pr[inside == 1] = NA

Pr <- matrix(grid_data$pr, nrow = length(x_seq), ncol = length(y_seq))
# Pr2 <- matrix(grid_data2$pr, nrow = length(x_seq), ncol = length(y_seq))


# image(x_seq, y_seq, Pr, col = hcl.colors(100),
#       xlim = c(465000, 560000),
#       xlab = "UTM easting", ylab = "UTM northing",
#       main = "Pr(travelling → foraging)",
#       bty = "n", asp = 1)
# # add a legend
# legend("top", legend = round(seq(0, 1, length = 6),1),
#        fill = hcl.colors(100)[c(1, 1:5*20)], bty = "n", horiz = TRUE, border = NA)


library(grDevices)  # Needed for color functions
library(graphics)   # Needed for rasterImage


# pdf("./case_studies/figs/muskox_space.pdf", width = 7, height = 7)
# png("./case_studies/figs/muskox_space.png", width = 3000, height = 3000, res = 450)

# Create the main plot
image(x_seq, y_seq, Pr, col = hcl.colors(100),
      xlim = c(465000, 560000), ylim = c(8230000, 8300000),
      xlab = "UTM easting", ylab = "UTM northing",
      main = expression("Pr(travelling" %->% "foraging)"),
      bty = "n", asp = 1)

# Define gradient legend position
legend_x <- c(480000, 545000)  # X range for legend
legend_y <- c(8300000-4000, 8300000-2000)  # Adjust above the plot

# Create a color gradient
legend_colors <- as.raster(matrix(hcl.colors(100), nrow = 1))

# Draw the gradient legend
rasterImage(legend_colors, legend_x[1], legend_y[1], legend_x[2], legend_y[2])

# Add numerical labels
legend_values <- round(seq(0, 1, length.out = 6), 1)
legend_positions <- seq(legend_x[1], legend_x[2], length.out = 6)
text(legend_positions, legend_y[2] + 2500, labels = legend_values, cex = 0.9)

# dev.off()



# Space-time interaction --------------------------------------------------

## creating model matrices
modmat2 = make_matrices(~ s(x, y, bs = "tp", k = 50) + s(tday, bs = "cc", k = 8) +
                         ti(x, y, tday, d = c(2,1), bs = c("tp", "cc"), k = c(50, 8)),
                       knots = list(tday = c(0, 24)),
                       data = data)
Z = modmat2$Z # large design matrix
S = modmat2$S # list of penalty matrices

# penalised negative log-likelihood function
pnll2 <- function(par){
  getAll(par, dat, warn = FALSE)

  ## transforming parameters
  mu = exp(logmu); REPORT(mu)
  sigma = exp(logsigma); REPORT(sigma)
  kappa = exp(logkappa); REPORT(kappa)

  # state process
  # beta = vector("list", 6)
  # beta[-4] = eta
  # beta[[4]] = c(beta0, beta_xy, beta_t, beta_xy_t)
  # Z_list = vector("list", 6)
  # Z_list[-4] = 1
  # Z_list[[4]] = Z
  beta <- matrix(0, 6, ncol(Z))
  beta[-4,1] <- eta
  beta[4,] <- c(beta0, beta_xy, beta_t, beta_xy_t)
  Gamma <- tpm_g(Z, beta)
  Delta <- stationary(Gamma[,,trackInd])

  ## computing state-dependent densities
  allprobs = matrix(1, nrow = length(step), ncol = N)
  ind = which(!is.na(step) & !is.na(angle))
  for(j in 1:N){
    allprobs[ind,j] = dgamma2(step[ind],mu[j],sigma[j]) *
      dvm(angle[ind], c(pi, 0, 0)[j], kappa[j])
  }
  ## forward algorithm
  -forward_g(Delta, Gamma, allprobs, trackID = ID) +
    penalty2(list(beta_xy, beta_t, beta_xy_t), S, lambda)
}


par = list(
  logmu = log(c(5, 50, 300)),
  logsigma = log(c(5, 50, 300)),
  logkappa = log(c(0.001, 0.1, 1)),
  eta = rep(-2, 5),
  beta0 = -2,
  beta_xy = rep(0, modmat2$pardim$`s(x,y)`),
  beta_t = rep(0, modmat2$pardim$`s(tday)`),
  beta_xy_t = rep(0, modmat2$pardim$`ti(x,y,tday)`)
)

dat = list(
  step = data$step,
  angle = data$angle,
  ID = data$ID,
  trackInd = calc_trackInd(data$ID),
  N = 3,
  Z = Z,
  S = S,
  lambda = c(1e3, 1e3, 1e4, 1e4)
)

system.time(
  mod_spacetime <- qreml(pnll2, par, dat,
                         random = c("beta_xy", "beta_t", "beta_xy_t"),
                         silent = 0)
)
# basis functions -> 8.7 h
saveRDS(mod_spacetime, file = "./case_studies/objects/muskox_mod_spacetime2.rds")

mod_spacetime = readRDS("./case_studies/objects/muskox_mod_spacetime.rds")

IC_table[5, 2:4] = c(mod_spacetime$llk, AIC(mod_spacetime), BIC(mod_spacetime))


par = mod_spacetime$par
eta = par$eta
beta = c(par$beta0, par$beta_xy, par$beta_t, par$beta_xy_t)

data$Pr23 = mod_spacetime$Gamma[3,2,]

states = viterbi_g(mod = mod_spacetime)
find = which(states == 2)


# plot the estimated field
# Define grid range
x_seq <- seq(min(data$x), max(data$x), length.out = 70)  # 100 grid points in x
y_seq <- seq(min(data$y), max(data$y), length.out = 70)  # 100 grid points in y
tday_seq <- 1:24 # 24 grid points in y

# Create a data frame with all combinations
grid_data <- expand.grid(x = x_seq, y = y_seq, tday = tday_seq)

hull_idx <- chull(data$x, data$y)  # Get convex hull indices
hull_coords <- cbind(data$x[hull_idx], data$y[hull_idx])  # Hull boundary points

# Check which grid points are inside the convex hull
inside_chull <- point.in.polygon(grid_data$x, grid_data$y, hull_coords[,1], hull_coords[,2])
inside_main = point.in.polygon(grid_data$x, grid_data$y, bnd_main$x, bnd_main$y)
inside_island = point.in.polygon(grid_data$x, grid_data$y, bnd_island$x, bnd_island$y)
inside = (inside_main | inside_island) & inside_chull

Z_pred = predict(modmat2, newdata = grid_data)

beta <- matrix(0, 6, ncol(Z))
beta[-4,1] <- eta
beta[4,] <- c(par$beta0, par$beta_xy, par$beta_t, par$beta_xy_t)
Gamma <- tpm_g(Z_pred, beta)


pr = Gamma[3,2,]
grid_data$pr = pr
# grid_data2 = grid_data
grid_data$pr[inside == 0] = NA
# grid_data2$pr[inside == 1] = NA

par(mfrow = c(1,1))
for(t in 1:24){
  this_grid = grid_data[grid_data$tday == t,]
  thisPr <- matrix(this_grid$pr, nrow = length(x_seq), ncol = length(y_seq))

  image(x_seq, y_seq, thisPr, col = hcl.colors(100),
        xlim = c(465000, 560000), ylim = c(8230000, 8300000),
        xlab = "UTM easting", ylab = "UTM northing",
        main = paste0("Pr(travelling → foraging)", " - time of day: ", t, ":00"),
        bty = "n", asp = 1)

  Sys.sleep(0.25)
}

# plot
par(mfrow = c(1,2))

for(t in c(1, 13)){
  this_grid = grid_data[grid_data$tday == t,]
  thisPr <- matrix(this_grid$pr, nrow = length(x_seq), ncol = length(y_seq))

  image(x_seq, y_seq, thisPr, col = hcl.colors(100),
        xlim = c(465000, 560000), ylim = c(8230000, 8300000),
        xlab = "UTM easting", ylab = "UTM northing",
        main = paste0("time of day: ", t, ":00"),
        bty = "n", asp = 1)

  # Define gradient legend position
  legend_x <- c(480000, 545000)  # X range for legend
  legend_y <- c(8300000-4000, 8300000-2000)  # Adjust above the plot

  # Create a color gradient
  legend_colors <- as.raster(matrix(hcl.colors(100), nrow = 1))

  # Draw the gradient legend
  rasterImage(legend_colors, legend_x[1], legend_y[1], legend_x[2], legend_y[2])

  # Add numerical labels
  legend_values <- round(seq(0, 1, length.out = 6), 1)
  legend_positions <- seq(legend_x[1], legend_x[2], length.out = 6)
  text(legend_positions, legend_y[2] + 2500, labels = legend_values, cex = 0.9)
}




# Additive model ----------------------------------------------------------


## creating model matrices
modmat3 = make_matrices(~ s(x, y, bs = "tp", k = 50) + s(tday, bs = "cc", k = 10),
                        knots = list(tday = c(0, 24)),
                        data = data)
Z = modmat3$Z # large design matrix
S = modmat3$S # list of penalty matrices

# penalised negative log-likelihood function
pnll3 <- function(par){
  getAll(par, dat, warn = FALSE)

  ## transforming parameters
  mu = exp(logmu); REPORT(mu)
  sigma = exp(logsigma); REPORT(sigma)
  kappa = exp(logkappa); REPORT(kappa)

  # state process
  beta <- matrix(0, 6, ncol(Z))
  beta[-4,1] <- eta
  beta[4,] <- c(beta0, beta_xy, beta_t)
  Gamma <- tpm_g(Z, beta)
  Delta <- stationary(Gamma[,,trackInd])

  ## computing state-dependent densities
  allprobs = matrix(1, nrow = length(step), ncol = N)
  ind = which(!is.na(step) & !is.na(angle))
  for(j in 1:N){
    allprobs[ind,j] = dgamma2(step[ind],mu[j],sigma[j]) *
      dvm(angle[ind], c(pi, 0, 0)[j], kappa[j])
  }
  ## forward algorithm
  -forward_g(Delta, Gamma, allprobs, trackID = ID) +
    penalty2(list(beta_xy, beta_t), S, lambda)
}


par = list(
  logmu = log(c(5, 50, 300)),
  logsigma = log(c(5, 50, 300)),
  logkappa = log(c(0.001, 0.1, 1)),
  eta = rep(-2, 5),
  beta0 = -2,
  beta_xy = rep(0, modmat3$pardim$`s(x,y)`),
  beta_t = rep(0, modmat3$pardim$`s(tday)`)
)

dat = list(
  step = data$step,
  angle = data$angle,
  ID = data$ID,
  trackInd = calc_trackInd(data$ID),
  N = 3,
  Z = Z,
  S = S,
  lambda = c(1e3, 1e3)
)

# system.time(
#   mod_add <- qreml(pnll3, par, dat,
#                    random = c("beta_xy", "beta_t"),
#                    silent = 0)
# )
# saveRDS(mod_add, file = "./case_studies/objects/muskox_mod_add.rds")

mod_add = readRDS("./case_studies/objects/muskox_mod_add.rds")
IC_table[4, 2:4] = c(mod_add$llk, AIC(mod_add), BIC(mod_add))




# only tday ---------------------------------------------------------------


## creating model matrices
modmat4 = make_matrices(~ s(tday, bs = "cc", k = 10),
                        knots = list(tday = c(0, 24)),
                        data = data)
Z = modmat4$Z # large design matrix
S = modmat4$S # list of penalty matrices

# penalised negative log-likelihood function
pnll4 <- function(par){
  getAll(par, dat, warn = FALSE)

  ## transforming parameters
  mu = exp(logmu); REPORT(mu)
  sigma = exp(logsigma); REPORT(sigma)
  kappa = exp(logkappa); REPORT(kappa)

  # state process
  beta <- matrix(0, 6, ncol(Z))
  beta[-4,1] <- eta
  beta[4,] <- c(beta0, beta_t)
  Gamma <- tpm_g(Z, beta)
  Delta <- stationary(Gamma[,,trackInd])

  ## computing state-dependent densities
  allprobs = matrix(1, nrow = length(step), ncol = N)
  ind = which(!is.na(step) & !is.na(angle))
  for(j in 1:N){
    allprobs[ind,j] = dgamma2(step[ind],mu[j],sigma[j]) *
      dvm(angle[ind], c(pi, 0, 0)[j], kappa[j])
  }
  ## forward algorithm
  -forward_g(Delta, Gamma, allprobs, trackID = ID) +
    penalty2(list(beta_t), S, lambda)
}


par = list(
  logmu = log(c(5, 50, 300)),
  logsigma = log(c(5, 50, 300)),
  logkappa = log(c(0.001, 0.1, 1)),
  eta = rep(-2, 5),
  beta0 = -2,
  beta_t = rep(0, modmat3$pardim$`s(tday)`)
)

dat = list(
  step = data$step,
  angle = data$angle,
  ID = data$ID,
  trackInd = calc_trackInd(data$ID),
  N = 3,
  Z = Z,
  S = S,
  lambda = 1e3
)

# system.time(
#   mod_time <- qreml(pnll4, par, dat,
#                     random = "beta_t",
#                     silent = 0)
# )
# saveRDS(mod_time, file = "./case_studies/objects/muskox_mod_time.rds")
mod_time = readRDS("./case_studies/objects/muskox_mod_time.rds")

IC_table[3, 2:4] = c(mod_time$llk, AIC(mod_time), BIC(mod_time))



diff_IC_table <- IC_table
colmin = apply(diff_IC_table[,3:4], 2, min)
diff_IC_table[,3:4] <- t(t(diff_IC_table[,3:4]) - colmin)

round(IC_table, 2)
round(diff_IC_table, 2)
