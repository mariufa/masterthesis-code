library(MCMCpack)
library(MASS)

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
  } else if(phiOption == "x1x2divX3Option") {
    return(getPhiValue(x[1]*x[2]/x[3]))
  } else if(phiOption == "x1divx2powx3Option") {
    return(getPhiValue((x[1]/x[2])^x[3]))
  } else if(phiOption == "sinusfunction") {
    return(sin(x[1]) + sin(x[2]) + sin(x[3]))
  }
  return(-1)
}


getPhiValue <- function(xValue) {
  # Calculates phi for an element of an x vector.
  # 
  # Args:
  #   xValue: A scalar value.
  #   
  # Returns:
  #   
  return(xValue > probValue)
}

optimInvGammaCumulative <- function(u, alpha) {
  result = optim(c(10.1), invGammaAbsFunction, u=u, alpha=alpha, lower=0, upper=100,  method="Brent")
  return(result)
}

invGammaAbsFunction <- function(x, u, alpha) {
  return(abs(pgamma(x, shape=alpha, scale=1) - u))
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

optimfindAlpha <- function(u, s2) {
  if((calcValueTau2(u, alphaUpperBound) < s2) || (calcValueTau2(u, alphaLowerBound) > s2)) {
    return(-1)
  }
  solution = optim(c(0.1), optimFunction, u=u, lower=alphaLowerBound, upper=alphaUpperBound, method="Brent")
  return(solution$par)
}

optimFunction <- function(alpha, u) {
  return(abs(s2-calcValueTau2(u, alpha)))
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
  NUM_ITERATIONS = 5000
  xCurrent = xInit
  for(i in 1:NUM_ITERATIONS) {
    randomXpos = sample(length(xCurrent), size=3)
    sumX = xCurrent[randomXpos[1]] + xCurrent[randomXpos[2]] + xCurrent[randomXpos[3]]
    prodX = xCurrent[randomXpos[1]]*xCurrent[randomXpos[2]]*xCurrent[randomXpos[3]]
    x1 = runif(1)*sumX
    if(isValidX1Proposal(x1, sumX, prodX)) {
      roots = findRoots(x1, sumX, prodX)
      x2 = roots[1]
      x3 = roots[2]
      if(is.nan(x2) || is.nan(x3)){
        print((x1^3 - 2*sumX*x1^2 + (sumX^2)*x1 - 4*prodX))
      }
      xProposal = xCurrent
      xProposal[randomXpos[1]] = x1
      xProposal[randomXpos[2]] = x2
      xProposal[randomXpos[3]] = x3
      alphaMetHastings = findAlphaMetHastings(c(xCurrent[randomXpos[1]], xCurrent[randomXpos[2]], xCurrent[randomXpos[3]]), c(xProposal[randomXpos[1]], xProposal[randomXpos[2]], xProposal[randomXpos[3]]))
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


GammaGibbsSampling <- function(xInit) {
  # A Gibbs sampling for a Gamma distribution
  # 
  # Args:
  #   xInit: Initial sample. A vector.
  #   
  # Return:
  #   A sample from a Gamma distribution. A vector.
  NUM_ITERATIONS = 5000
  xCurrent = xInit
  for(i in 1:NUM_ITERATIONS) {
    randomXpos = sample(length(xCurrent), size=3)
    sumX = xCurrent[randomXpos[1]] + xCurrent[randomXpos[2]] + xCurrent[randomXpos[3]]
    prodX = xCurrent[randomXpos[1]]*xCurrent[randomXpos[2]]*xCurrent[randomXpos[3]]
    x1 = runif(1)*sumX
    if(isValidX1Proposal(x1, sumX, prodX)) {
      roots = findRoots(x1, sumX, prodX)
      x2 = roots[1]
      x3 = roots[2]
      if(is.nan(x2) || is.nan(x3)){
        print((x1^3 - 2*sumX*x1^2 + (sumX^2)*x1 - 4*prodX))
      }
      xProposal = xCurrent
      xProposal[randomXpos[1]] = x1
      xProposal[randomXpos[2]] = x2
      xProposal[randomXpos[3]] = x3
      alphaMetHastings = 0
      alphaMetHastings = findGammaAlphaMetHastings(c(xCurrent[randomXpos[1]], xCurrent[randomXpos[2]], xCurrent[randomXpos[3]]), c(xProposal[randomXpos[1]], xProposal[randomXpos[2]], xProposal[randomXpos[3]]))
      acceptProb = runif(1)
      if(!isValidX1Proposal(xCurrent[randomXpos[1]], sumX, prodX)) {
        alphaMetHastings = 0
      }
      if(acceptProb <= alphaMetHastings) {
        xCurrent = xProposal
      }
    }
  }
  return(xCurrent)
  
}

findGammaAlphaMetHastings <- function(xCurrent, xProposal) {
  piProp = 1/(xProposal[1]*sqrt((sum(xProposal)  - xProposal[1])^2 - 4*prod(xProposal)/xProposal[1] ))
  piCurrent = 1/(xCurrent[1]*sqrt((sum(xCurrent)  - xCurrent[1])^2 - 4*prod(xCurrent)/xCurrent[1] ))
  if(is.nan(piCurrent)) {
    print(xCurrent)
  }
  return(min(1, piProp/piCurrent))
}


cramerVonMisesValueTest <- function(x, alpha, beta) {
  # Calculates the value for a Cramer-von Mises test.
  # 
  # Args:
  #   x: A vector sample.
  #   
  # Returns:
  #   A scalar value.
  cramerSum = 0
  for(i in 1:length(x)) {
    cramerSum = cramerSum + ((2*i - 1)/(2*length(x)) - pgamma(x[i], shape=alpha, rate=beta))^2
  }
  cramer = 12/(length(x)) + cramerSum
  return(cramer)
}

findGammaMLE <- function(x){
  solution = optim(c(1,1), negativeLogLikelihoodGamma, x=x)
  return(solution$par)
}

negativeLogLikelihoodGamma <- function(par, x) {
  # Calculates the negative log-likelihood for a gamma distribution.
  # 
  # Args:
  #   par: A vector of size 2. First element is alpha and second element is beta.
  #   x: Data to calculate log-likelihood from. A vector.
  #   
  # Returns:
  #  The log-likelihood value. A scalar
  alpha = par[1]
  beta = par[2]
  logLikelihood = -((alpha - 1)*sum(log(x)) - (1/beta)*sum(x) - length(x)*log(gamma(alpha)) - alpha*length(x)*log(beta))
  return(logLikelihood)
}

calcAveragPhiValueForData <- function(mydata) {
  sumData = sum(mydata)
  prodData = prod(mydata)
  tolerance = 0.03
  minValue = min(mydata) - 2*tolerance
  maxValue = max(mydata) + 2*tolerance
  sampleNumber = 1
  NUM_ITERATIONS = 10000
  sumPhi = 0
  while(sampleNumber <= NUM_ITERATIONS) {
    x = runif(3, max = sumData)
    if((abs(sum(x) - sumData) < tolerance) && (abs(prod(x) - prodData) < tolerance)) {
      sumPhi = sumPhi + calcPhiGivenX(x)
      sampleNumber = sampleNumber + 1  
      print(sampleNumber)
    }
  }
  return(sumPhi/(sampleNumber-1))
}

algorithm2Sampling <- function(NUM_ALG2_SAMPLES) {
  cramerNum = 0
  cramerStat = rep(0, NUM_ALG2_SAMPLES)
  vCurr = runif(NUM_POINTS)
  alphaCurr = optimfindAlpha(vCurr, s2)
  while(alphaCurr==-1) {
    vCurr = runif(NUM_POINTS)
    alphaCurr = optimfindAlpha(vCurr, s2)
  }
  piCurr = calcWeight(vCurr, alphaCurr)
  phiSum = 0
  
  for(i in 1:NUM_ALG2_SAMPLES) {
    print(i)
    vProp = runif(NUM_POINTS)
    alphaProp = optimfindAlpha(vProp, s2)
    piProp = 0
    if(alphaProp != -1) {
      piProp = calcWeight(vProp, alphaProp)
    }
    alphaMetHastings = min(1, piProp/piCurr)
    uProb = runif(1)
    if(uProb <= alphaMetHastings) {
      vCurr = vProp
      alphaCurr = alphaProp
      piCurr = piProp
    }
    betaCurr = findBeta(s1, vCurr, alphaCurr)
    xSample = rep(0, length(vCurr))
    for(j in 1:length(vCurr)) {
      xSample[j] = betaCurr*invGammaCumulative(vCurr[j], alphaCurr)
    }
    phiSum = phiSum + calcPhiGivenX(xSample)
    cramerStat[i] = cramerVonMisesValueTest(xSample, mleAlpha, mleBeta)
    if(cramerStat[i] >= cramerObs) {
      cramerNum = cramerNum + 1
    }
  }
  
  hist(cramerStat,breaks=200, main="", xlab="Cramer values", cex.lab=1.5)
  abline(v = cramerObs, col="red")
  alg2sampCramer <<- cramerStat
  
  print((cramerNum/NUM_ALG2_SAMPLES))
  alg2pvalue <<- (cramerNum/NUM_ALG2_SAMPLES)
  return(phiSum/NUM_ALG2_SAMPLES)
}

algorithm1Sampling <- function(NUM_ALG1_SAMPLES) {
  phiSum = 0
  cramerNum = 0
  cramerStat = rep(0, NUM_ALG1_SAMPLES)
  for(i in 1:NUM_ALG1_SAMPLES) {
    print(i)
    u = runif(NUM_POINTS)
    alphavalue = optimfindAlpha(u, s2)
    while(alphavalue == -1) {
      u = runif(NUM_POINTS)
      alphavalue = optimfindAlpha(u, s2)
    }
    betavalue = findBeta(s1, u, alphavalue)
    xSample = rep(0, length(u))
    for(j in 1:length(u)) {
      xSample[j] = betavalue*invGammaCumulative(u[j], alphavalue)
    }
    phiSum = phiSum + calcPhiGivenX(xSample)
    cramerStat[i] = cramerVonMisesValueTest(xSample, mleAlpha, mleBeta)
    #print(cramerStat[i])
    if(cramerStat[i] >= cramerObs) {
      cramerNum = cramerNum + 1
    }
  }
  print("here")
  print((cramerNum/NUM_ALG1_SAMPLES))
  hist(cramerStat,breaks=200, main="", xlab="Cramer values", cex.lab=1.5)
  abline(v = cramerObs, col="red")
  alg1sampCramer <<- cramerStat
  alg1pvalue <<- (cramerNum/NUM_ALG1_SAMPLES)
  return(phiSum/NUM_ALG1_SAMPLES)
}

naiveSampling2 <- function(myData, tolerance, mleAlpha, mleBeta) {
  NUM_NAIVE_SAMPLES = NUM_SAMPLES
  sumData = sum(myData)
  prodData = prod(myData)
  sampleNumber = 0
  sumPhi = 0
  iterations = 0
  cramerNum = 0
  cramerStat = rep(0, NUM_NAIVE_SAMPLES)
  while(sampleNumber<NUM_NAIVE_SAMPLES) {
    x = rgamma(3, shape = mleAlpha, rate = mleBeta)
    iterations = iterations + 1
    if((abs(sum(x) - sumData) < tolerance) && (abs(prod(x) - prodData) < tolerance)) {
      sumPhi = sumPhi + calcPhiGivenX(x)
      sampleNumber = sampleNumber + 1  
      cramerStat[sampleNumber] = cramerVonMisesValueTest(x, mleAlpha, mleBeta)
      if(cramerStat[sampleNumber] >= cramerObs) {
        cramerNum = cramerNum + 1
      }
      print(sampleNumber)
    }
  }
  
  hist(cramerStat,breaks=200, main="", xlab="Cramer values", cex.lab=1.5)
  abline(v = cramerObs, col="red")
  
  acceptRate = sampleNumber/iterations
  averagePhi = sumPhi/sampleNumber
  print((cramerNum/NUM_NAIVE_SAMPLES))
  naivePvalue <<- (cramerNum/NUM_NAIVE_SAMPLES)
  naiveCramers <<- cramerStat
  return(c(acceptRate, averagePhi))
}

alpha = 1
beta = 1
hStep = 0.01
alphaHStep = 0.01
method = "pgamma"
NUM_SAMPLES = 100000
NUM_POINTS = 3
alphaUpperBound = 200
alphaLowerBound = 0.05
# Pi is used in calculation of weights
# Options are:
#   "constant"
#   "betaOption"
#   "jeffrey"
#   "alphaOption"
piValue = "constant"

phi = rep(0, NUM_SAMPLES)
# Phi options:
# x larger than a: "probValueOption"
# x1 times x2 div x3: "x1x2divX3Option"
# x1 div x2 pow x3: "x1divx2powx3Option"
phiOption = "x1x2divX3Option"
# Phi is the prob that X>probValue
probValue = 0.5

# Data generation options:
# pgamma generated: "pgamma"
# Bo data: "bo"
# Custom data: "custom"
# Custom data2: "custom2"
dataGenOption = "custom"



# Generate data
gammaData = 0
if(dataGenOption == "pgamma") {
  gammaData = rgamma(NUM_POINTS, shape=alpha, scale = beta)  
} else if(dataGenOption == "bo") {
  NUM_POINTS = 6
  alphaUpperBound = 1.2
  alphaLowerBound = 0.8
  gammaData = c(4.399, 1.307, 0.085, 0.7910, 0.2345, 0.1915)
} else if(dataGenOption == "custom") {
  gammaData = c(0.5772030, 0.4340237, 0.4212959)
} else if(dataGenOption == "custom2") {
  gammaData = c(1.621813, 1.059797, 1.554334)
} else if(dataGenOption == "custom3") {
  gammaData = c(5, 15, 13)
} else if(dataGenOption == "ball1") {
  gammaData = c(17.88, 28.92, 33.00)
} else if(dataGenOption == "ball2") {
  gammaData = c(41.52, 42.12, 45.60)
} else if(dataGenOption == "ball3") {
  gammaData = c(48.40, 51.84, 51.96)
} else if(dataGenOption == "ball4") {
  gammaData = c(54.12, 55.56, 67.80)
} else if(dataGenOption == "ball5") {
  gammaData = c(68.64, 68.64, 68.88)
} else if(dataGenOption == "ball6") {
  gammaData = c(84.12, 93.12, 98.64)
} else if(dataGenOption == "ball7") {
  gammaData = c(105.12, 105.84, 127.92)
} else if(dataGenOption == "custom4") {
  gammaData = c(0.40, 0.42, 0.43)
} else if(dataGenOption == "custom5") {
  gammaData = c(0.72, 0.72, 0.85)
}

#17.88 28.92 33.00 41.52 42.12 45.60 48.40 51.84
#51.96 54.12 55.56 67.80 68.64 68.64 68.88 84.12
#93.12 98.64 105.12 105.84 127.92 128.04 173.40

hist(gammaData)
# Calculation of statistics
s1 = sum(gammaData)/NUM_POINTS
s2 = NUM_POINTS*((prod(gammaData))^(1/NUM_POINTS))/sum(gammaData)
# Log-likelihood
mleEstimators = fitdistr(gammaData, "gamma")
mleAlpha = mleEstimators$estimate[1]
mleBeta = mleEstimators$estimate[2]
# w statistic obs. Not to be used yet.
wObs = calcPhiGivenX(gammaData)
cramerObs = cramerVonMisesValueTest(gammaData, mleAlpha, mleBeta)


# Calc Phi value for data
#phiValue = calcAveragPhiValueForData(gammaData)
tolerance = 0.005
naiveCramers = rep(0,100000)
naivePvalue = 0
naiveSampler = naiveSampling2(gammaData, tolerance, mleAlpha, mleBeta)

## Tolerance accuracy plot
#toleranceRange = seq(0.02, 2, 0.01)
#resultAcceptance = rep(0, length(toleranceRange))
#resultValue = rep(0, length(toleranceRange))
#i = 0
#for(tol in toleranceRange) {
#  i = i+1
#  result = naiveSampling2(gammaData, tol)
#  resultAcceptance[i] = result[1]
#  resultValue[i] = result[2]
#}
#plot(toleranceRange, resultValue)
#plot(toleranceRange, resultAcceptance)

# Generation of samples. Not to be used

# weightsW = rep(0, NUM_SAMPLES)
# sampleIndex = 1
# iterationNumber = 0
# estAlpha = rep(0, NUM_SAMPLES)
# estBeta = 0
# while(sampleIndex <= NUM_SAMPLES) {
#   u = runif(NUM_POINTS)
#   estAlpha[sampleIndex] = optimfindAlpha(u, s2)
#   if(estAlpha[sampleIndex] != -1) {
#     estBeta = findBeta(s1, u, estAlpha[sampleIndex])
#     if(estBeta > alphaLowerBound && estBeta < alphaUpperBound) {
#         
#       weightsW[sampleIndex] = abs(calcWeight(u, estAlpha[sampleIndex]))
#       phi[sampleIndex] = calcPhi(u, estAlpha[sampleIndex])
#       
#       sampleIndex = sampleIndex + 1
#       print(sampleIndex)
#     }
#     
#   }
#   #print(iterationNumber)
#   iterationNumber = iterationNumber + 1
# }
# alphaAcceptance = (sampleIndex-1)/iterationNumber
# hist(weightsW, breaks = 400)
# expectedPhi = sum(phi*weightsW)/sum(weightsW)
# plot(estAlpha, weightsW)
# unweightedExpectedPhi = sum(phi)/NUM_SAMPLES

# Gibbs sampling
NUM_GIBBS_SAMPLES = NUM_SAMPLES
xSample = gammaData
phiGibbs = rep(0, NUM_GIBBS_SAMPLES)
gibbsObslargerWObs = 0
cramerNum = 0
cramerStatGibbs = rep(0, NUM_GIBBS_SAMPLES)
for(i in 1:NUM_GIBBS_SAMPLES) {
  xSample = gibbsSampling(xSample)
  phiGibbs[i] = calcPhiGivenX(xSample)
  if(phiGibbs[i] >= wObs) {
    gibbsObslargerWObs = gibbsObslargerWObs + 1
  }
  cramerStatGibbs[i] = cramerVonMisesValueTest(xSample, mleAlpha, mleBeta)
  #print(cramerStat)
  if(cramerStatGibbs[i] >= cramerObs) {
    cramerNum = cramerNum + 1
  }
  print(i)
}
gibbsS1 = sum(xSample)/NUM_POINTS
gibbsS2 = NUM_POINTS*((prod(xSample))^(1/NUM_POINTS))/sum(xSample)
gibbsPvalue = gibbsObslargerWObs/NUM_GIBBS_SAMPLES
averagePhiGibbs = sum(phiGibbs)/NUM_GIBBS_SAMPLES
## P-values
cramerPValue = cramerNum/NUM_GIBBS_SAMPLES

hist(cramerStatGibbs,breaks=200, main="", xlab="Cramer values", cex.lab=1.5)
abline(v = cramerObs, col="red")

alg2sampCramer = rep(1, 100000)
alg1sampCramer = rep(1, 100000)
alg2pvalue = 0
alg1pvalue = 0

system.time({alg1Results2 = algorithm1Sampling(NUM_SAMPLES)})
system.time({alg2Results2 = algorithm2Sampling(NUM_SAMPLES)})

#library("parallel")
#library("foreach")
#library("doParallel")
#
#cl = makeCluster(detectCores() - 1)
#registerDoParallel(cl, cores = detectCores() - 1)
#workers = 10
#stime = system.time({
#res = foreach(i=1:workers, 
#        .combine = rbind) %dopar% {
#          try({
#            result1 = algorithm1Sampling(100000/workers)
#          })
#        }
#})
#stopCluster(cl)
#
#alg1Results3 = (sum(res[,1]))/workers
#
#cl = makeCluster(detectCores() - 1)
#registerDoParallel(cl, cores = detectCores() - 1)
#workers = 10
#stime = system.time({
#  res = foreach(i=1:workers, 
#                .combine = rbind) %dopar% {
#                  try({
#                    result1 = algorithm2Sampling(100/workers)
#                  })
#                }
#})
#stopCluster(cl)
#
#alg2Results3 = (sum(res[,1]))/workers
#
#
#
## Save image
save.image(file="scen4.RData")

#print("Done")
#
## Plot of phi functions
#sumData = 4.23
#prodData = 2.67
#x1Seq = seq(1.01, 2, 0.001)
#phi2 = rep(0,length(x1Seq))
#phi3 = rep(0,length(x1Seq))
#for(i in 1:length(x1Seq)) {
#  if(isValidX1Proposal(x1Seq[i], sumData, prodData)) {
#    roots = findRoots(x1Seq[i], sumData, prodData)
#    phi2[i] = phi2Function(x1Seq[i], roots[1], roots[2])
#    phi3[i] = phi3Function(x1Seq[i], roots[1], roots[2])
#  }
#  
#}
#
#phi2Function <- function(x1, x2, x3) {
#  return((x1*x2/x3))
#}
#
#phi3Function <- function(x1, x2, x3) {
#  return((x1/x2)^x3)
#}
#
#plot(x1Seq, phi2)
#plot(x1Seq, phi3)
#
## Plot of F inverse
#u = 0.5
#alphaRange = seq(0.01, 4, 0.001) 
#fInv = rep(0, length(alphaRange))
#for(i in 1:length(fInv)) {
#  fInv[i] = invGammaCumulative(u, alphaRange[i])
#}
#plot(alphaRange, fInv, type="l")
#
#diffF = rep(0, length(fInv))
#for(i in 1:length(diffF)) {
#  diffF[i] = diffAlphaInvGammaCumulative(u, alphaRange[i])
#}
#
#plot(alphaRange, diffF, type="l")
#
#
#globalvarTest = matrix(0, 10, 10)
#
#testFunction <- function(i) {
#  globalvarTest[,i] <<- rep(i,10)
#  return(0)
#}
#
#cl = makeCluster(detectCores() - 1)
#registerDoParallel(cl, cores = detectCores() - 1)
#workers = 10
#stime = system.time({
#  res = foreach(i=1:workers, 
#                .combine = rbind) %dopar% {
#                  try({
#                    testFunction(i)
#                  })
#                }
#})
#stopCluster(cl)
#
