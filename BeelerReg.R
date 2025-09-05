# BeelerReg.R
# by Russell R. Barton, September 2025
#
# Extra source code needed by testRGP to plot regression approximation to Beeler
# Apply DOE to Beeler.R f function then construct regression predictor function BeelerReg
# Uses first part of InterpModelResults.R and testRGP variables:
# fname,nlev,ylab,xlab,xlo,xhi
fr <- function(x){
  source("Beeler.R")
  xlo = c(-1,-1,0,0)
  xhi = c(1,1,2,2)
  ylab = "Cases"
  xlab = NULL
  xlab[1] = "Days"
  xlab[2] = "Infec"
  xlab[3] = "Vac"
  xlab[4] = "Quar"
  nlev = c(3,3,3,3)
  
  #first check that nlev,xlo and xhi all have length d <= 6
  d = length(nlev)
  
  #construct factorial levels, one row for each x 
  xlev = matrix(nrow=d,ncol=max(nlev))
  for (i in 1:d){
    for (j in 1:nlev[i]){
      xlev[i,j] = xlo[i] + (j-1)*(xhi[i]-xlo[i])/(nlev[i]-1)
    }
  }
 
  #construct factorial design with appropriate levels and labels x1, x2, etc.
  library(DoE.base)
  # factorial design
  DOE = data.matrix(fac.design(factor.names=xlab[1:d],nlevels=nlev))
  
  #rescale DOE values to x and compute corresponding responses, factorial
  xDOE = matrix(0,nrow=nrow(DOE),ncol=ncol(DOE))
  for (i in 1:nrow(DOE)){
    for (j in 1:d){
      xDOE[i,j] = xlo[j] + (DOE[i,j]-1)*(xhi[j]-xlo[j])/(nlev[j]-1)
    }
  }
  
 colnames(xDOE) = xlab
 

  #compute response function for each row of the xDOE
  resp = apply(xDOE,1,f)
  # order variables based on number of levels
  #attach response as last column
  DOE = as.data.frame(cbind(xDOE,resp))
  
# add model fits and characterizations here
  
# linear regression
  linRegModel = lm(resp ~ (Days + Infec + Quar + Vac)^3 +
                     I(Days^2) + I(Infec^2) + I(Quar^2) + I(Vac^2), data=DOE)
  
# compute regression predictor at x
  xPred = data.frame(Days = x[1],Infec = x[2], Quar = x[3], Vac = x[4])
  return = predict(linRegModel, newdata = xPred)
    
  }
