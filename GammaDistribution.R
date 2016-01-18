library(MCMCpack)


calcWeight <- function(u, alpha, beta) {
  
  
}

invGammaCumulative <- function(u, alpha) {
  x = 0
  stepSize = 0.1
  tolerance = 0.0001
  direction = 1
  integralValue = 0
  
  while(abs(integralValue - u) > tolerance){
    x = x + direction*stepSize
    if(x<0) {
      x = 0
      direction = 1
    }
    
    integral = integrate(gammaDensity, 0, x)
    integralValue = integral$value
    
    # Going left and pass the point
    if ((u > integralValue) && (direction == -1)) {
      stepSize = stepSize/2
      direction = 1
    }
    
    # Going right and pass the point
    if ((u < integralValue) && (direction == 1)) {
      stepSize = stepSize/2
      direction = -1
    }
  }
  return(x)
}

diffInvGammaCumulative <- function() {
  
}

gammaDensity <- function(x) {
  return(dgamma(x, alpha,1))  
}

alpha = 2

NUM_SAMPLES = 1000
xsamp = rep(0, NUM_SAMPLES)
for(i in 1:NUM_SAMPLES) {
  u = runif(1)
  xsamp[i] = invGammaCumulative(u, alpha)
}
hist(xsamp)
xdata = dgamma(seq(0,8, by=0.1),2)
plot(xdata)

u = runif()
