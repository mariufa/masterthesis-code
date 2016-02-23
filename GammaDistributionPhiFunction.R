library(MCMCpack)


calcWeight <- function(u, alpha) {
  # Calculates weight for given u and alpha.
  # 
  # Args:
  #   u: A vector.
  #   alpha: A scalar.
  #   
  # Returns:
  #   The weight value. A scalar.
  gammaInv = rep(0, length(u))
  diffGammaInv = rep(0, length(u))
  for(i in 1:length(u)) {
    gammaInv[i] = invGammaCumulative(u[i], alpha)
    diffGammaInv[i] = diffAlphaInvGammaCumulative(u[i], alpha)
  }
  pi = getPiValue()
  weight = pi/((1/length(gammaInv))*(sum(diffGammaInv/gammaInv)) - sum(diffGammaInv)/sum(gammaInv))
  return(weight)
}

getPiValue <- function() {
  # get value of pi function to be used in calculation of weights.
  # 
  # Returns:
  #   A scalar value.
  if (piValue == "constant") {
    return(1)
  } else if (piValue == "jeffrey") {
    # Return jeffrey prior
    return(sqrt((1/(estAlpha[sampleIndex]^2)) + (1/(2-(exp(estAlpha[sampleIndex]) - exp(-estAlpha[sampleIndex]))))))
  } else if (piValue == "betaOption") {
    return(estBeta)
  } else if (piValue == "alphaOption") {
    return(estAlpha[sampleIndex])
  }
}

calcPhi <- function(u, alpha) {
  phiPoint = rep(0, length(u))
  for(i in 1:length(u)) {
    xValue = invGammaCumulative(u[i], alpha)
    phiPoint[i] = getPhiValue(xValue)
  }
  return(sum(phiPoint)/length(u))
}

calcPhiGivenX <- function(x) {
  # Calc phi value for a vector x.
  # 
  # Args:
  #   x: A vector of data.
  #   
  # Returns:
  #  A scalar.
  phiPoint = rep(0, length(x))
  for(i in 1:length(x)) {
    phiPoint[i] = getPhiValue(x[i])
  }
  return(sum(phiPoint)/length(phiPoint))
}

getPhiValue <- function(xValue) {
  return(xValue > probValue)
}

invGammaCumulative <- function(u, alpha) {
  # Finds the inverse of the cumulative gamma.
  # 
  # Args:
  #   u: Scalar value.
  #   alpha: Scalar value.
  #   
  # Returns:
  #   A scalar. The inverse value given u and alpha.
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
      integralv = integrate(gammaDensity, 0, x)
      integralValue = integral$valuev
    } else if(method == "pgamma") {
      integralv = pgamma(x, shape=alpha, scale=1)
      integralValue = integralv
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

diffAlphaInvGammaCumulative <- function(u, alpha) {
  # Derivative of gamma distribution with respect to alpha.
  # 
  # Args:
  #   u: Scalar between 0 and 1.
  #   alpha: Scalar larger than 0.
  #   
  # Returns:
  #   Scalar value. Derivative at point alpha.
  firstPoint = invGammaCumulative(u,alpha)
  secondPoint = invGammaCumulative(u, alpha + alphaHStep)
  return((secondPoint - firstPoint)/alphaHStep)
}

calcDerivateFunction <- function(u, alphaValue) {
  # Analytically find the derivative of cumulative inverse.
  # 
  # Args:
  #   u: Scalar value between 0 and 1.
  #   alphaValue: Scalar value.
  #   
  # Returns:
  #   A scalar value.
  largeFInv = invGammaCumulative(u, alphaValue)
  integralPart = integrate(integralFunction, 0, largeFInv)
  return((digamma(alphaValue)*u - integralPart$value)*gamma(alphaValue)/((largeFInv^(alphaValue - 1))*exp(-largeFInv)))
}

integralFunction <- function(y) {
  # Function to be integrated.
  return(log(y)*(y^(alphaValue - 1))*exp(-y)/gamma(alphaValue))
}

optimfindAlpha <- function() {
  if((calcValueTau2(u, alphaUpperBound) < s2) || (calcValueTau2(u, alphaLowerBound) > s2)) {
    return(-1)
  }
  return(optim(c(0.1), optimFunction, lower=alphaLowerBound, upper=alphaUpperBound, method="Brent")$par)
}

optimFunction <- function(alpha) {
  return(abs(s2-calcValueTau2(u, alpha)))
}

findAlpha <- function(s2, u) {
  # Function to calculate alpha value
  # 
  # Args:
  #   s2: Scalar value of sufficient statistic
  #   u: A vector
  #   
  # Returns:
  #   A scalar value of alpha.
  
  alphaValue = 0.1
  stepSize = 1
  direction = 1
  tolerance = 0.0001
  tau2Value = 0
  prevTau2Value = 0
  it = 0
  while(abs(s2 - tau2Value) > tolerance) {
    it = it + 1
    alphaValue = alphaValue + direction*stepSize
    if(alphaValue<=0) {
      alphaValue = 0.1
      direction = 1
    }
    prevTau2Value = tau2Value
    tau2Value = calcValueTau2(u, alphaValue)
    if ((s2 > tau2Value) && (direction == -1)) {
      stepSize = stepSize/2
      direction = 1
    }
    
    if ((s2 < tau2Value) && (direction == 1)) {
      stepSize = stepSize/2
      direction = -1
    }
    
    if(isAlphaOutsideValidInterval(alphaValue, direction)) {
      return(-1)
    }  
    if(it==100) {
      return(-1)
    }
    
  }
  
  return(alphaValue)
}

isAlphaOutsideValidInterval <- function(alphaValue, direction) {
  # Function to check if the alphaValue is outside allowed interval.
  # 
  # Args:
  #   alphaValue: A scalar value
  #   direction: A scalar value. Either 1 or -1.
  #   
  # Returns:
  #   Boolean value.
  return((alphaValue < alphaLowerBound && direction == -1) || (alphaValue > alphaUpperBound && direction == 1))
}

findBeta <- function(s1, u, alphaValue) {
  # Calculates beta
  # 
  # Args:
  #   s1: Scalar value.
  #   u: Vector of values between 0 and 1.
  #   alphaValue: Scalar value.
  #   
  # Returns:
  #   A scalar value.
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
beta = 2
hStep = 0.01
alphaHStep = 0.01
method = "pgamma"
NUM_SAMPLES = 10000
NUM_POINTS = 3
alphaUpperBound = 20
alphaLowerBound = 0.05
# Pi is used in calculation of weights
# Options are:
#   "constant"
#   "betaOption"
#   "jeffrey"
#   "alphaOption"
piValue = "jeffrey"

# Generate data
gammaData = rgamma(NUM_POINTS, shape=alpha, scale = beta)
hist(gammaData)
# Calculation of statistics
s1 = sum(gammaData)/NUM_POINTS
s2 = NUM_POINTS*((prod(gammaData))^(1/NUM_POINTS))/sum(gammaData)
# w statistic obs. Not to be used yet.
wObs = calcPhiGivenX(gammaData)

# Generation of samples
phi = rep(0, NUM_SAMPLES)
# Phi is the prob that X>probValue
probValue = 1
weightsW = rep(0, NUM_SAMPLES)
sampleIndex = 1
iterationNumber = 0
estAlpha = rep(0, NUM_SAMPLES)
estBeta = 0
while(sampleIndex <= NUM_SAMPLES) {
  u = runif(NUM_POINTS)
  estAlpha[sampleIndex] = optimfindAlpha()
  if(estAlpha[sampleIndex] != -1) {
    estBeta = findBeta(s1, u, estAlpha[sampleIndex])
    weightsW[sampleIndex] = abs(calcWeight(u, estAlpha[sampleIndex]))
    phi[sampleIndex] = calcPhi(u, estAlpha[sampleIndex])
    
    sampleIndex = sampleIndex + 1
    print(sampleIndex)
    
  }
  #print(iterationNumber)
  iterationNumber = iterationNumber + 1
}
alphaAcceptance = (sampleIndex-1)/iterationNumber
hist(weightsW, breaks = 300)
expectedPhi = sum(phi*weightsW)/sum(weightsW)
plot(estAlpha, weightsW)
unweightedExpectedPhi = sum(phi)/NUM_SAMPLES


