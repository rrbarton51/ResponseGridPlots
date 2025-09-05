# interpFit.R
# by Russell R. Barton, August 2025
#
# Given d and a function of x in R^d, construct a factorial grid
#  with number of levels nlev[1] ... nlev[d] and 
#  x limits xlo[1] ... xlo[d] and xhi[1] ... xhi[d]
# Variable names must be provided in a label vector xlab

interpFit <- function(f,fname,nlev,ylab,xlab,xlo,xhi){
  #first check that nlev,xlo and xhi all have length d <= 6
  d = length(nlev)
  if(d > 6) stop("dimension of x must be 6 or less")
  if(length(xlo) != d) stop(" dimensions for levels and limits don't match)")
  if(length(xhi) != d) stop(" dimensions for levels and limits don't match)")
  if(max(nlev) > 10) stop("maximum levels is 10")
  if(sum(nlev>5)>3) stop(" at most 3 x can have more than 5 levels")
  if(sum(nlev>3)>4) stop(" at most 4 x can have more than 3 levels")
  
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
  
  # LH design
  #set.seed(12357)
  #nruns = prod(nlev[1:d])
  #library(lhs)
  #DOE = data.matrix(create_oalhs(n = nruns, k = d, 
  #                               bChooseLargerDesign = FALSE, bverbose = FALSE))
  #colnames(DOE) = xlab[1:d]
  
  #rescale DOE values to x and compute corresponding responses, factorial
  xDOE = matrix(0,nrow=nrow(DOE),ncol=ncol(DOE))
  for (i in 1:nrow(DOE)){
    for (j in 1:d){
      xDOE[i,j] = xlo[j] + (DOE[i,j]-1)*(xhi[j]-xlo[j])/(nlev[j]-1)
    }
  }
  
  #rescale to x values oalhs
  #rescale DOE values to x and compute corresponding responses, oalhs
  #DOE = matrix(0,nrow=nrow(DOE),ncol=ncol(DOE))
  #for (i in 1:nrow(DOE)){
  #  for (j in 1:d){
  #    xDOE[i,j] = xlo[j] + (DOE[i,j])*(xhi[j]-xlo[j])
  #  }
  #}
  
 colnames(xDOE) = xlab
 
 

  #compute response function for each row of the xDOE
  resp = apply(xDOE,1,f)
  # order variables based on number of levels
  xDOE = xDOE[,order(nlev,decreasing=TRUE)]
  # alternative for oalhs design fitted to linRegModel
  #resp = predict(linRegModel,newdata=as.data.frame(xDOE))
  
  #attach response as last column
  DOE = as.data.frame(cbind(xDOE,resp))
  
  print(xDOE)
  print(DOE)
  
# add model fits and characterizations here
  
# linear regression
  linRegModel = lm(resp ~ (Days + Infec + Quar + Vac)^3 +
                     I(Days^2) + I(Infec^2) + I(Quar^2) + I(Vac^2), data=DOE)
  print(linRegModel)
  summary(linRegModel)
  # now model with just .05 significant terms
  linRegModel2 = lm(resp ~ Days + Infec + Quar + Vac +
                      + I(Infec^2) + I(Quar^2) + I(Vac^2) +
                      Days:Infec + Days:Vac + Infec:Quar + Infec:Vac +
                      Quar:Vac + Days:Infec:Quar + Days:Infec:Vac +
                      Infec:Quar:Vac, data=DOE)
  print(linRegModel2)
  summary(linRegModel2)
  
# apply rtree to the virus data
  library(caret)
  library(rpart)
  library(rpart.plot) # for plotting the tree

  #trainMeth = trainControl(method = "cv", number = 10)
  #regTreeModel <- train(resp ~ Days + Infec + Quar + Vac, data = DOE, 
  #                   method = "rpart",
  #                   trControl = trainMeth)
  #NOte: above gave best cp = .032 and R-squared of .55
  # could not plot from caret - use rpart directly
  set.seed(12357)
  regTreeModel = rpart(resp ~ Days + Infec + Quar + Vac, data = DOE,
                       method = "anova", 
                       control = rpart.control(minsplit = 3, maxdepth = 4, cp = 0.005))
  print(regTreeModel)
  summary(regTreeModel)
  rpart.plot(regTreeModel)
  #get SSTOT and SSE for R-squared
  #get predictions
  predictions <- predict(regTreeModel, newdata = DOE)
  SSEtree = sum((predictions-resp)^2)
  SStot = sum((resp - mean(resp))^2)
  treeRsquared = 1-(SSEtree/SStot)
  print(c("tree R-squared",treeRsquared))
  
  
  library(RWeka) # for M5 model
  #library(partykit) # for plot
  M5model <- M5P(resp ~ ., data = DOE, control = Weka_control(N = FALSE, M = 3))
  summary(M5model)
  print(M5model)
  #partyM5model = as.party(M5model)
  #plot(partyM5model)
  #get SSTOT and SSE for R-squared
  #get predictions
  predictions <- predict(regTreeModel, newdata = DOE)
  SSEM5 = sum((predictions-resp)^2)
  SStot = sum((resp - mean(resp))^2)
  M5Rsquared = 1-(SSEM5/SStot)
  print(c("M5 R-squared",M5Rsquared))

# Construct ICE plots for x1 - x4 about scaled (0,0,0,0), need iml
library(iml)
library(ggplot2) # included with iml but theme_minimal feature used below

  predFunc = Predictor$new(linRegModel,data = DOE)

# Compute the accumulated local effects (ALE) for the Days
aleDays = FeatureEffect$new(predFunc, feature = "Days", grid.size = 30)
p1 = aleDays$plot() + theme_minimal()
print(p1)

# Compute ICE and PDP for the Days
iceDays = FeatureEffect$new(predFunc, feature = "Days", method = "pdp+ice", grid.size = 30)
p2 = iceDays$plot()+theme_minimal()
print(p2)

# Compute the accumulated local effects (ALE) for the Infec
aleInfec = FeatureEffect$new(predFunc, feature = "Infec", grid.size = 30)
p3 = aleInfec$plot() + theme_minimal()
print(p3)

# Compute ICE and PDP for the Infec
iceInfec = FeatureEffect$new(predFunc, feature = "Infec", method = "pdp+ice", grid.size = 30)
p4 = iceInfec$plot()+theme_minimal()
print(p4)


# Compute the accumulated local effects (ALE) for the Quar
aleQuar = FeatureEffect$new(predFunc, feature = "Quar", grid.size = 30)
p5 = aleQuar$plot() + theme_minimal()
print(p5)

# Compute ICE and PDP for the Quar
iceQuar = FeatureEffect$new(predFunc, feature = "Quar", method = "pdp+ice", grid.size = 30)
p6 = iceQuar$plot()+theme_minimal()
print(p6)

# Compute the accumulated local effects (ALE) for the Vac
aleVac = FeatureEffect$new(predFunc, feature = "Vac", grid.size = 30)
p7 = aleVac$plot() + theme_minimal()
print(p7)

# Compute ICE and PDP for the Vac
iceVac = FeatureEffect$new(predFunc, feature = "Vac", method = "pdp+ice", grid.size = 30)
p8 = iceVac$plot()+theme_minimal()
print(p8)

#Compute permutation feature importance
# First, create a large training data set from the model
#predFunc2 = Predictor$new(linRegModel,data = largeDOE)

impRank = FeatureImp$new(predFunc, loss = "mse")
print(impRank)
p9 = impRank$plot()+theme_minimal()
plot(p9)

# Friedman interaction plots
# overall
interactOverall = Interaction$new(predFunc)
p10 = interactOverall$plot() + theme_minimal()
plot(p10)
# for variable with largest overall interaction
interactInfec = Interaction$new(predFunc, feature = "Infec")
p11 = interactInfec$plot() + theme_minimal()
plot(p11)
# Tornado plot
library(tornado)
tornPlot <- tornado::tornado(linRegModel, type = "ranges")
plot(tornPlot, xlabel = "Total Infections", geom_bar_control = list(width = 0.4))

}
