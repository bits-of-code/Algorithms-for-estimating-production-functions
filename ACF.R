# 3. Manually-coded Cobb-Douglas production function estimation procedure using the ACF (2015) correction 

# Conduct the first-stage regression
cf_u_modul <- cf_u_modul %>% filter(vEnergy > 0)
cf_u_modul <- cf_u_modul[ order(cf_u_modul$unr, cf_u_modul$year) ,]

prod.ACF <- lm(log(va) ~ poly(cbind(log(vCap), log(Lab), log(vEnergy)), 
                                            degree = 4) + year + region + sector, data = cf_u_modul)

# Conduct the second-stage regression
phi_hat <- predict(prod.ACF, cf_u_modul)

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
                               lhs = log(cf_u_modul$va),
                               k = log(cf_u_modul$vCap), 
                               l = log(cf_u_modul$Lab), 
                               phi_hat = phi_hat, 
                               k.lag = log(lag(cf_u_modul$vCap)),
                               l.lag = log(lag(cf_u_modul$Lab)),
                               phi.lag = lag(phi_hat) )

# Drop missing observations
cf_u_modul.step2 <- subset(cf_u_modul.step2, !apply(cf_u_modul.step2, 1, function(x) any(is.na(x))))

# Construct objective function
objective <- function(x, degree = 4) {
  
  beta_vCap = x[1]
  beta_Lab = x[2]
  
  prod.ACF <- lm(I(phi_hat - beta_vCap * k - beta_Lab * l) ~ poly(I(phi.lag - beta_vCap * k.lag - 
                                                                    beta_Lab * l.lag), degree), data = cf_u_modul.step2)
  return(sum(residuals(prod.ACF)^2))
}

# Solve the minimization problem
opt.out <- optim(c(0.5, 0.5), fn = objective)

# Predict omega 
cf_u_modul["omega_ACF"] <- with(cf_u_modul,
                                log(va) - (cbind(log(vCap), log(Lab)) %*% opt.out$estimate[1:2]))
compute_hist(cf_u_modul, omega_ACF)
