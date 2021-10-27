# 2. Manually-coded Cobb-Douglas production function estimation procedure using the Levinsohn & Petrin (2003) approach 

# Conduct the first-stage regression
cf_u_modul <- cf_u_modul %>% filter(vEnergy > 0)
cf_u_modul <- cf_u_modul[order(cf_u_modul$unr, cf_u_modul$year), ]

prod.LP <- lm(log(vOut) ~ log(Lab) + poly(cbind(log(vCap), log(vEnergy)), 
                                          degree = 4) + year + region + sector, data = cf_u_modul)

# Conduct the second-stage regression
betas <- prod.LP$coefficients[c("log(Lab)")]
xbetas <- log(as.matrix(cf_u_modul[, c("Lab")])) %*% betas
phi_hat <- predict(prod.LP, cf_u_modul) - xbetas

# Define lag function
lag <- function(x, i = cf_u_modul$unr, t = cf_u_modul$year) {
  if (length(i) != length(x) | length(i) != length(t)) {
    stop("Inputs not same length")
  }
  x.lag <- x[1:(length(x) - 1)]
  x.lag[i[1:(length(i) - 1)] != i[2:length(i)]] <- NA
  x.lag[t[1:(length(i) - 1)] + 1 != t[2:length(i)]] <- NA
  return(c(NA, x.lag))
}

# Create dataframe for second stage regression
cf_u_modul.step2 <- data.frame( 
                                lhs = log(cf_u_modul$vOut) - xbetas, 
                                k = log(cf_u_modul$vCap), 
                                M = log(cf_u_modul$vEnergy), 
                                phi_hat = phi_hat, 
                                k.lag = log(lag(cf_u_modul$vCap)), 
                                phi.lag = lag(phi_hat), 
                                M.lag = log(lag(cf_u_modul$vEnergy)) )

# Drop missing observations
cf_u_modul.step2 <- subset(cf_u_modul.step2, !apply(cf_u_modul.step2, 1, function(x) any(is.na(x))))

# Construct objective function
objective <- function(x, degree = 4) {
  
  beta_vCap = x[1]
  beta_vEnergy = x[2]
  
  prod.LP <- lm(I(lhs - beta_vCap * k - beta_vEnergy * M) ~ poly(I(phi.lag - beta_vCap * k.lag - 
                                                                beta_vEnergy * M.lag), degree), data = cf_u_modul.step2)
  return(sum(residuals(prod.LP)^2))
}

# Solve the minimization problem
opt.out <- nlm( objective, c(0.0, 0.0) )

# Predict omega
cf_u_modul["omega_LP"] <- with(cf_u_modul,
                               log(vOut) - (log(Lab) * coef(prod.LP)[2] - log(vCap) * opt.out$estimate[1]))
compute_hist(cf_u_modul, omega_LP)
