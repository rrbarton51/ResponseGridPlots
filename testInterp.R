# testInterp.R
# Russell R. Barton September 2025
# example program to define a function and levels for plotting
#  then call interpFit - for interpretative model fits and output 

# series of test functions can be place here
# Only Beeler implemented for the GitHub archive


# Beeler et al. (2012) disease infection response
source("interpFit.R")
source("Beeler.R")
fname = "Beeler"

xlo = c(-1,-1,0,0)
xhi = c(1,1,2,2)
ylab = "Cases"
xlab = NULL
xlab[3] = "Vac"
xlab[4] = "Quar"
xlab[1] = "Days"
xlab[2] = "Infec"
nlev = c(3,3,3,3)

interpFit(f,fname,nlev,ylab,xlab,xlo,xhi)

