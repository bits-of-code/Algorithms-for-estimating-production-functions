# 1. Manually-coded Cobb-Douglas production function estimation procedure using the Olley & Pakes (1996) approach 

# Drop zero values
cf_u_modul <- cf_u_modul %>% filter( vInv > 0 )

# Conduct the first-stage regression
cf_u_modul <- cf_u_modul[ order(cf_u_modul$unr, cf_u_modul$year), ]

prod.op <- lm( log( vOut ) ~ log( Lab ) + log( vMat ) +
           poly( cbind( log( vCap ), log( vInv ) ), degree = 4 ) + year + region + sector,
           data = cf_u_modul )

# Calculate fhat
betas <- prod.op$coefficients[c("log(Lab)")]
xbetas <- log(as.matrix(cf_u_modul[,c("Lab")])) %*% betas
phi_hat <- predict(prod.op, cf_u_modul) - xbetas

# Create lagged values of fhat and vCap
lag <- function(x, i = cf_u_modul$unr, t = cf_u_modul$year ) {
  if ( length(i) != length(x) | length (i) != length(t) ) {
    stop ("Inputs not same length")
  }
  x.lag <- x[1:(length(x) - 1)]
  x.lag[i [1:(length(i) - 1)] != i[2:length(i)] ] <- NA
  x.lag[t [1:(length(i) - 1) ] + 1 != t [2:length(i)] ] <- NA
  return( c(NA, x.lag) )
}

# Create dataframe for second-stage regression
cf_u_modul.step2 <- data.frame(
                           lhs = log(cf_u_modul$vOut) - xbetas, 
                           k = log(cf_u_modul$vCap), 
                           phi_hat = phi_hat, 
                           k.lag = log(lag(cf_u_modul$vCap)), 
                           phi.lag = lag(phi_hat) )

# Drop missing observations
cf_u_modul.step2 <- subset(cf_u_modul.step2, !apply(cf_u_modul.step2, 1, function(x) any(is.na(x))))

# Construct objective function
objective <- function(beta_vCap, degree = 4) {
  
  prod.op <- lm(I(lhs - beta_vCap * k) ~ poly(I(phi.lag - beta_vCap * k.lag), degree),
             data = cf_u_modul.step2)
  
  return( sum(residuals(prod.op)^2) )
}

# Solve the minimization problem
opt.out <- optim(0.0,
                 fn = objective, method = "Brent",lower = -1, upper = 1)

# Predict omega
cf_u_modul["omega_OP"] <- with(cf_u_modul,
                               log(vOut) - (log(vCap) * opt.out$par[1] + log(Lab) * coef(prod.op)[2] 
                                            + coef(prod.op)['year'] + coef(prod.op)['region']))
compute_hist(cf_u_modul, omega_OP)
