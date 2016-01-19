library(MCMCpack)


calcWeight <- function(gammaInv, diffGammaInv) {
  weight = (1/length(gammaInv))*(sum(diffGammaInv/gammaInv)) - sum(diffGammaInv)/sum(gammaInv)
  return(weight)
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
    
    if(method == "integrate") {
      integral = integrate(gammaDensity, 0, x)
      integralValue = integral$value  
    } else if(method == "pgamma") {
      integral = pgamma(x, alpha)
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

alpha = 2
hStep = 0.01
method = "integrate"

NUM_SAMPLES = 1000
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
