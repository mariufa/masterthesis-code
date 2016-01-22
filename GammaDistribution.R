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
      direction = -1
    }
    
    print(s2 - tau2Value)
    tau2Value = calcValueTau2(u, alphaValue)
    
    if ((s2 > tau2Value) && (direction == 1)) {
      stepSize = stepSize/2
      direction = -1
    }
    
    if ((s2 < tau2Value) && (direction == -1)) {
      stepSize = stepSize/2
      direction = 1
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

calcValueTau2 <- function(u, alpha) {
  largeFInv = rep(0, length(u))
  for(i in 1:length(u)) {
    largeFInv = invGammaCumulative(u[i], alpha)
  }
  return(length(u)*((prod(largeFInv))^(1/length(u)))/sum(largeFInv))
}

alpha = 2
beta = 2
hStep = 0.004
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
  print(i)
}
plot(tau2[50:length(tau2)], type="l")


# Generation of new sample
u = runif(NUM_POINTS)
estAlpha = findAlpha(s2, u)
estBeta = findBeta(s1, u, estAlpha)


# Testing
xsamp = rep(0, NUM_SAMPLES)
for(i in 1:NUM_SAMPLES) {
  u = runif(1)
  xsamp[i] = invGammaCumulative(u, alpha)
}
hist(xsamp)
xdata = dgamma(seq(0,8, by=0.1),2)
plot(xdata)

u = seq(0,1,by=0.01)
largeF = rep(0,length(u))
diffF = rep(0,length(u)-1)
for(i in 1:length(u)) {
  largeF[i] = invGammaCumulative(u[i],alpha)
  if(i<length(u)) {
    diffF[i] = diffInvGammaCumulative(u[i])  
  }
  print(i)
}
plot(u,largeF)
lines(u[1:100],diffF)

NUM_POINTS = 20
weights = rep(0, NUM_SAMPLES)
for(i in 1:NUM_SAMPLES) {
  u = runif(NUM_POINTS)
  diffF = rep(0, NUM_POINTS)
  invF = rep(0, NUM_POINTS)
  for(j in 1:NUM_POINTS) {
    invF[j] = invGammaCumulative(u[j], alpha)
    diffF[j] = diffInvGammaCumulative(u[j])
  }
  weights[i] = calcWeight(invF, diffF)
  print(i)
}

hist(weights, breaks = 200)
