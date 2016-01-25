library(MCMCpack)


calcWeight <- function(gammaInv, diffGammaInv) {
  weight = (1/length(gammaInv))*(sum(diffGammaInv/gammaInv)) - sum(diffGammaInv)/sum(gammaInv)
  return(weight)
}

invGammaCumulative <- function(u, alpha) {
  x = 0
  stepSize = 0.1
  tolerance = 0.00001
  direction = 1
  integralValue = 0
  
  while(abs(integralValue - u) > tolerance){
    x = x + direction*stepSize
    if(x<0) {
      x = 0
      direction = 1
    }
    
    if(method == "integrate") {
      integral = integrate(gammaDensity, 0, x)
      integralValue = integral$value  
    } else if(method == "pgamma") {
      integral = pgamma(x, shape=alpha, scale=1)
      integralValue = integral
    }
    
    
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

diffInvGammaCumulative <- function(u) {
  if(u+hStep < 1) {
    firstPoint = invGammaCumulative(u, alpha)
    secondPoint = invGammaCumulative(u+hStep, alpha)
    return((secondPoint - firstPoint)/hStep)  
  } else {
    firstPoint = invGammaCumulative(u-hStep, alpha)
    secondPoint = invGammaCumulative(u, alpha)
    return((secondPoint - firstPoint)/hStep)  
  }
}

calcDerivateFunction <- function(u, alphaValue) {
 largeFInv = invGammaCumulative(u, alphaValue)
 integralPart = integrate(integralFunction, 0, largeFInv)
 return((digamma(alphaValue)*u - integralPart$value)*gamma(alphaValue)/((largeFInv^(alphaValue - 1))*exp(-largeFInv)))
}

integralFunction <- function(y) {
  return(log(y)*(y^(alphaValue - 1))*exp(-y)/gamma(alphaValue))
}

gammaDensity <- function(x) {
  return(dgamma(x, alpha,1))  
}

findAlpha <- function(s2, u) {
  alphaValue = 0.1
  stepSize = 10
  direction = 1
  tolerance = 0.0001
  tau2Value = 0
  
  while(abs(s2 - tau2Value) > tolerance) {
    alphaValue = alphaValue + direction*stepSize
    if(alphaValue<=0) {
      alphaValue = 0.1
      direction = 1
    }
    
    tau2Value = calcValueTau2(u, alphaValue)
    
    if ((s2 > tau2Value) && (direction == -1)) {
      stepSize = stepSize/2
      direction = 1
    }
    
    if ((s2 < tau2Value) && (direction == 1)) {
      stepSize = stepSize/2
      direction = -1
    }
    
  }
  return(alphaValue)
}

findBeta <- function(s1, u, alphaValue) {
  largeFInv = rep(0, length(u))
  for(i in 1:length(u)) {
    largeFInv[i] = invGammaCumulative(u[i], alphaValue)
  }
  return(s1*length(u)/(sum(largeFInv)))
}

calcValueTau2 <- function(u, alphaValue) {
  largeFInv = rep(0, length(u))
  for(i in 1:length(u)) {
    largeFInv[i] = invGammaCumulative(u[i], alphaValue)
  }
  return(length(u)*((prod(largeFInv))^(1/length(u)))/sum(largeFInv))
}

calcValueTau2Method2 <- function(x) {
  return(length(x)*((prod(x))^(1/length(x)))/sum(x))
}

alpha = 2
beta = 1
hStep = 0.01
method = "pgamma"
NUM_SAMPLES = 1000
NUM_POINTS = 100

# Generate data
gammaData = rgamma(NUM_POINTS, shape=alpha, scale = beta)
hist(gammaData)
# Calculation of statistics
s1 = sum(gammaData)/NUM_POINTS
s2 = NUM_POINTS*((prod(gammaData))^(1/NUM_POINTS))/sum(gammaData)

# Plot of tau 2 with respect to alpha
alpharange = seq(0.1 , 10, by = 0.1)
tau2 = rep(0, length(alpharange))
u = runif(NUM_POINTS)
for(i in 1:length(alpharange)) {
  tau2[i]= calcValueTau2(u, alpharange[i])
}
plot(tau2, type="l")

# Plot of derivatve of inverse cumulative distribution and the cumulative function
alphaValue = alpha
urange = seq(0.01, 0.99, by=0.01)
diff1 = rep(0, length(urange))
diff2 = rep(0, length(urange))
invF = rep(0, length(urange))
for(i in 1:length(urange)) {
  invF[i] = invGammaCumulative(urange[i], alpha)
  diff1[i] = calcDerivateFunction(urange[i], alphaValue)
  diff2[i] = diffInvGammaCumulative(urange[i])
}
plot(urange, invF, type="l", ylim=c(0, 15))
lines(urange, diff1, col="red")
lines(urange, diff2, col="green")


# Generation of new sample
u = runif(NUM_POINTS)
estAlpha = findAlpha(s2, u)
estBeta = findBeta(s1, u, estAlpha)


