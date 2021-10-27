# 4. Manually-coded Cobb-Douglas production function estimation procedure using the Collard-Wexler & De Loecker (2016) correction 

# Conduct the first-stage regression
cf_u_modul <- cf_u_modul %>% filter(vInv > 0 & vEnergy > 0)
cf_u_modul <- cf_u_modul[ order(cf_u_modul$unr, cf_u_modul$year) ,]

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

# Create dataframe for first stage regression
cf_u_modul.step1 <- data.frame(
  n = cf_u_modul$unr,
  t = cf_u_modul$year,
  y = log(cf_u_modul$va),
  i = log(cf_u_modul$vInv), 
  l = log(cf_u_modul$Lab),
  d = log(cf_u_modul$vEnergy),
  i.lag = log(lag(cf_u_modul$vInv)),
  l.lag = log(lag(cf_u_modul$Lab))
  )

# Drop missing observations
cf_u_modul.step1 <- subset(cf_u_modul.step1, !apply(cf_u_modul.step1, 1, function(x) any(is.na(x))))

prod.CWD <- lm(y ~ poly(cbind(i.lag, l, d), 
                              degree = 4) + t, data = cf_u_modul.step1)

# Conduct the second-stage regression
phi_hat <- predict(prod.CWD, cf_u_modul.step1)

# Define lag function
lag <- function(x, i = cf_u_modul.step1$n, t = cf_u_modul.step1$t) {
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
  n = cf_u_modul.step1$n,
  t = cf_u_modul.step1$t,
  y = cf_u_modul.step1$y,
  i = cf_u_modul.step1$i,
  l = cf_u_modul.step1$l, 
  phi_hat = phi_hat, 
  i.lag = cf_u_modul.step1$i.lag,
  i.lag2 = lag(cf_u_modul.step1$i.lag),
  l.lag = cf_u_modul.step1$l.lag,
  phi.lag = lag(phi_hat) )

# Drop missing observations
cf_u_modul.step2 <- subset(cf_u_modul.step2, !apply(cf_u_modul.step2, 1, function(x) any(is.na(x))))

# Construct objective function
objective <- function(x, degree = 4) {
  
  beta_vCap = x[1]
  beta_Lab = x[2]
  
  prod.CWD <- lm(I(phi_hat - beta_vCap * i.lag - beta_Lab * l) ~ poly(I(phi.lag - beta_vCap * i.lag - 
                                                                      beta_Lab * l.lag), degree), data = cf_u_modul.step2)
  return(sum(residuals(prod.CWD)^2))
}

# Solve the minimization problem
opt.out <- optim(c(0.5, 0.5), fn = objective)

# Predict omega 
cf_u_modul["omega_ACF"] <- with(cf_u_modul,
                                log(va) - (cbind(log(vCap), log(Lab)) %*% opt.out$estimate[1:2]))
compute_hist(cf_u_modul, omega_ACF)
