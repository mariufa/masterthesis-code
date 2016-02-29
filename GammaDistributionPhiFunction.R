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
    return(sqrt((1/(estAlpha[sampleIndex]^2)) + (1/(2-(exp(estAlpha[sampleIndex]) + exp(-estAlpha[sampleIndex]))))))
  } else if (piValue == "betaOption") {
    return(estBeta)
  } else if (piValue == "alphaOption") {
    return(estAlpha[sampleIndex])
  }
}

calcPhi <- function(u, alpha) {
  xValue = rep(0, length(u))
  for(i in 1:length(u)) {
    xValue[i] = invGammaCumulative(u[i], alpha)
  }
  return(calcPhiGivenX(xValue))
}

calcPhiGivenX <- function(x) {
  # Calc phi value for a vector x.
  # 
  # Args:
  #   x: A vector of data.
  #   
  # Returns:
  #  A scalar.
  
  if(phiOption == "probValueOption") {
    phiPoint = rep(0, length(x))
    for(i in 1:length(x)) {
      phiPoint[i] = getPhiValue(x[i])
    }
    return(sum(phiPoint)/length(phiPoint))  
  } 
  return(-1)
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

gibbsSampling <- function(xInit) {
  sumX = sum(xInit)
  prodX = prod(xInit)
  NUM_ITERATIONS = 5000
  xCurrent = xInit
  for(i in 1:NUM_ITERATIONS) {
    x1 = runif(1)*sumX
    if(isValidX1Proposal(x1, sumX, prodX)) {
      roots = findRoots(x1, sumX, prodX)
      x2 = roots[1]
      x3 = roots[2]
      xProposal = c(x1, x2, x3)
      alphaMetHastings = findAlphaMetHastings(xCurrent, xProposal)
      acceptProb = runif(1)
      if(acceptProb <= alphaMetHastings) {
        xCurrent = xProposal
      }
    }
  }
  return(xCurrent)
}

isValidX1Proposal <- function(x1, sumX, prodX) {
  return((x1^3 - 2*sumX*x1^2 + (sumX^2)*x1 - 4*prodX) > 0)
}

findRoots <- function (x1, sumX, prodX) {
  root1 = ((sumX - x1) + sqrt( (sumX - x1)^2 - 4*prodX/x1 ))/2
  root2 = ((sumX - x1) - sqrt( (sumX - x1)^2 - 4*prodX/x1 ))/2
  return(c(root1, root2))
}

findAlphaMetHastings <- function(xCurrent, xProposal) {
  piProp = 1/(xProposal[1]*sqrt((sum(xProposal)  - xProposal[1])^2 - 4*prod(xProposal)/xProposal[1] ))
  piCurrent = 1/(xCurrent[1]*sqrt((sum(xCurrent)  - xCurrent[1])^2 - 4*prod(xCurrent)/xCurrent[1] ))
  return(min(1, piProp/piCurrent))
}

alpha = 2
beta = 2
hStep = 0.01
alphaHStep = 0.01
method = "pgamma"
NUM_SAMPLES = 1000
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
# Phi options:
# x larger than a: "probValueOption2
# x1 times x2 div x3: "x1x2divX4Option"
phiOption = "probValueOption"
# Phi is the prob that X>probValue
probValue = 2


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

# Gibbs sampling
NUM_GIBBS_SAMPLES = 1000
xSample = gammaData
phiGibbs = rep(0, NUM_GIBBS_SAMPLES)
gibbsObslargerWObs = 0
for(i in 1:NUM_GIBBS_SAMPLES) {
  xSample = gibbsSampling(xSample)
  phiGibbs[i] = calcPhiGivenX(xSample)
  if(phiGibbs[i] >= wObs) {
    gibbsObslargerWObs = gibbsObslargerWObs + 1
  }
  print(i)
}
gibbsPvalue = gibbsObslargerWObs/NUM_GIBBS_SAMPLES
averagePhiGibbs = sum(phiGibbs)/NUM_GIBBS_SAMPLES


