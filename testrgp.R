# testrgp.R
# Russell R. Barton September 2025
# example program to define a function and levels for plotting
#  then call rgp - Response Grid Plot. See
#  Barton, R. R. (1998). Design-plots for factorial and fractional-factorial
#    designs. Journal of Quality Technology 30:1, 40-54. And 
#  Barton, R. R. (2025) Response Grid Plots for Model-Agnostic Machine Learning Insight

# Make sure plot window is large enough to accommodate plot legend

source("rgp.R")

# series of test functions follows

# d = 6
f <- function(x) 2 - 2*x[1]^2 - 15*x[2] - 30*x[3]*x[4]*x[5] - 
  20*x[1]*x[4] - 12*x[2]*x[4] - 25*x[1]*x[4]*x[2] + x[6]*x[5]
fname = "3^2 4^3 5^1 Quad Interaction"

xlo = -10:-5
xhi = rep(1,6)
ylab = "f"
xlab = NULL
xlab[1]="Varx1";xlab[2]="Varx2";xlab[3]="Varx3";xlab[4]="Varx4"
xlab[5]="Varx5";xlab[6]="Varx6"
nlev = c(3,3,4,5,4,4)

rgp(f,fname,nlev,ylab,xlab,xlo,xhi)


# d = 5
f <- function(x) 2 - 2*x[1]^2 - 15*x[2] - 30*x[3]*x[4]*x[5] -
  20*x[1]*x[4] - 12*x[2]*x[4] - 5*x[1]*x[4]*x[2] + x[5]^2
fname = "3^1 4^2 5^2 Quad Interaction"
xlo = rep(-5,5)
xhi = rep(1,5)
ylab = "g"
xlab = NULL
xlab[1]="Varx1";xlab[2]="Varx2";xlab[3]="Varx3";xlab[4]="Varx4"
xlab[5]="Varx5"
nlev = c(3,4,4,5,5)

rgp(f,fname,nlev,ylab,xlab,xlo,xhi)


# EOQ model 5 levels

f <- function(x){
  D = x[1]
  S = x[2]
  H = x[3]
  EOQ = sqrt(2*S*D/H)
}
fname = "EOQ 5lev"

xlo = c(35,8,4)
xhi = c(65,17,8.5)
ylab = "EOQ"
xlab = NULL
xlab[1]="D";xlab[2]="S";xlab[3]="H"
nlev = c(5,5,5)

rgp(f,fname,nlev,ylab,xlab,xlo,xhi)


# EOQ model 2 levels

f <- function(x){
  D = x[1]
  S = x[2]
  H = x[3]
  EOQ = sqrt(2*S*D/H)
}
fname = "EOQ 2lev"

xlo = c(40,10,5)
xhi = c(60,15,7.5)
ylab = "EOQ"
xlab = NULL
xlab[1]="D";xlab[2]="S";xlab[3]="H"
nlev = c(2,2,2)

rgp(f,fname,nlev,ylab,xlab,xlo,xhi)


# Beeler et al. Regression form rather than ANOVA
# (works with any nlev, here 3)

fname = "BeelerReg"

xlo = c(-1,-1,0,0)
xhi = c(1,1,2,2)
ylab = "Cases"
xlab = NULL
xlab[1] = "Days"
xlab[2] = "Infec"
xlab[3] = "Vac"
xlab[4] = "Quar"
nlev = c(3,3,3,3)
source("BeelerReg.R")

rgp(fr,fname,nlev,ylab,xlab,xlo,xhi)

# Beeler et al. Regression form rather than ANOVA
# (works with any nlev, here 4)

fname = "BeelerReg"

xlo = c(-1,-1,0,0)
xhi = c(1,1,2,2)
ylab = "Cases"
xlab = NULL
xlab[1] = "Days"
xlab[2] = "Infec"
xlab[3] = "Vac"
xlab[4] = "Quar"
nlev = c(4,4,4,4)
source("BeelerReg.R")

rgp(fr,fname,nlev,ylab,xlab,xlo,xhi)

