# ResponseGridPlots
This repository contains R files for codes used in "Response Grid Plots for Model-Agnostic Machine Learning Insight" by Russell R. Barton.

There are several code files, each are described below.

rgp.R - a code file for the function rgp(f,fname,nlev,ylab,xlab,xlo,xhi) to generate an RGP.
  f - function called by rgp, argument x a vector of dimension <= 6, result is f(x)=y,
  fname - text variable for name of the function for the plot title,
  nlev - vector of number of grid levels for each component of x,
  ylab - text variable with short name for the response,
  xlab - list of short text names for x[1], x[2] etc.,
  xlo - vector of lower limit for each component of x for the RGP,
  xhi - vector of upper limit for each component of x for the RGP.

testrgp.R - source code to exercise rgp on several functions at several grid settings.

Beeler.R - Beeler(x) codes the response for the original virus model (an ANOVA model) used in "Response Grid Plots for Model-Agnostic Machine Learning Insight".

BeelerReg.R codes the response for the polynomial regression virus model that appears in "Response Grid Plots for Model-Agnostic Machine Learning Insight".

rgp2.R - a code file for function rgp2(f1,fname1,f2,fname2,nlev,ylab,xlab,xlo,xhi) to generate a two-response RGP. 
  f1 - first response function in the two-response RGP,
  fname1 - text variable for name of the first function for the plot title,
  f2 - second response function in the two-response RGP,
  fname2 - text variable for name of the second function for the plot title,
  (other arguments the same as for rgp).

testrgp2.R - source code to exercise rgp2. For the repository, the code only contains the Khamkure arsenic desorption model, plotted for 3x3x3 grid and 5x5x5 grid.

interpFit.R - function that takes the same arguments as rgp, and creates a regression ML model and other interpretable ML models. The latter part of the code applies other explainability methods to the regression ML model.
NOTE: interpFit is not flexible code: it is only intended for the virus (Beeler) model.

testInterp.R - source code to exercise interpFit.R. It is only designed to work with the virus (Beeler) function.

