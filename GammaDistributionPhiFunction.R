library(MCMCpack)

weightedUniform <- function(n, numIterations) {
  # Generate sample of from a weighted uniform distribution
  # 
  # Args:
  #   n: number of points in sample
  #   
  # Returns:
  #   A vector sample
  uCurr = runif(n)
  alphaCurr = findAlpha(s2, uCurr)
  uConv = rep(0, numIterations)
  uTrace = rep(0, numIterations)
  for(i in 1:numIterations) {
    uProp = runif(n)
    alphaProp = findAlpha(s2, uProp)
    if(alphaProp != -1) {
      proportionalProp = proportionalDensity(uProp, alphaProp)
      proportionalCurr = proportionalDensity(uCurr, alphaCurr)
      alphaVal = min(1, proportionalProp/proportionalCurr)
      u = runif(1)
      print(i)
      if (u <= alphaVal) {
        uCurr = uProp
        alphaCurr = alphaProp
      }  
    }
    uConv[i] = var(uCurr)
    uTrace[i] = uCurr[1]
  }
  plot(uConv, type="l")
  plot(uTrace, type="l")
  return(uCurr)
}

proportionalDensity <- function(u, alpha) {
  # Calculate proportional density for indepence sampling
  # 
  # Args:
  #   u: A vector from an uniform distribution.
  #   alpha: A scalar value
  #   
  # Returns:
  #   A scalar value
  
  gammaInv = rep(0, length(u))
  diffGammaInv = rep(0, length(u))
  for(i in 1:length(u)) {
    gammaInv[i] = invGammaCumulative(u[i], alpha)
    diffGammaInv[i] = diffAlphaInvGammaCumulative(u[i], alpha)
  }
  return(calcWeight(gammaInv, diffGammaInv))
}

calcWeight <- function(gammaInv, diffGammaInv) {
  # Calculate weigth for a gamma distribution
  #
  # Args:
  #   gammaInv: A list of values from the inverse cumulative gamma.
  #   diffGammaInv: A list of values from the derivative of the inverse cumulative gamma.
  #   
  # Returns:
  #   A scalar value
  weight = (1/length(gammaInv))*(sum(diffGammaInv/gammaInv)) - sum(diffGammaInv)/sum(gammaInv)
  return(weight)
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
  # Derivative of gamma distribution with respect to u
  # 
  # Args:
  #   u: Scalar between 0 and 1
  #   
  # Returns:
  #   Scalar. Derivative at point u.
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
    #print(abs(s2 - tau2Value))
    #print(direction)
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
alphaHStep = 0.01
method = "pgamma"
NUM_SAMPLES = 1000
NUM_POINTS = 3
alphaUpperBound = 10
alphaLowerBound = 0.5

# Generate data
gammaData = rgamma(NUM_POINTS, shape=alpha, scale = beta)
hist(gammaData)
# Calculation of statistics
s1 = sum(gammaData)/NUM_POINTS
s2 = NUM_POINTS*((prod(gammaData))^(1/NUM_POINTS))/sum(gammaData)