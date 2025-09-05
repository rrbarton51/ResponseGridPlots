# testrgp2.R
# Russell R. Barton September 2025
# example program to define a function and levels for plotting
#  then call rgp2 - Response-Scaled Design-Plots for two functions. See
#  Barton, R. R. (1998). Design-plots for factorial and fractional-factorial
#    designs. Journal of Quality Technology 30:1, 40-54. 

#  Make sure plot window is large enough to accommodate legend

source("rgp2.R")

# series of test functions follows

# Khamkure original range 4 levels
f1 <- function(x) {
  xlo = c(3,.5,.05)
  xhi = c(7,4.5,.25)
  sx = -1+2*(x-xlo)/(xhi-xlo)
f1 = 51.97-7.42*sx[1]+13.37*sx[2]*sx[3]-18.23*sx[1]^2+15.54*sx[3]^2
  }

fname1 = "AS(III)Removal"


f2 <- function(x) {
  xlo = c(3,.5,.05)
  xhi = c(7,4.5,.25)
  sx = -1+2*(x-xlo)/(xhi-xlo)
  f2 = 28.62-23.43*sx[1]+10.81*sx[2]-11.99*sx[1]*sx[2]+26.08*sx[1]^2-9.52*sx[3]^2
}
fname2 = "AS(V)Removal"

xlo = c(3,.5,.05)
xhi = c(7,4.5,.25)
ylab = "%"
xlab = NULL
xlab[1]="pH";xlab[2]="Dose";xlab[3]="AsConc"
nlev = c(4,4,4)

rgp2(f1,fname1,f2,fname2,nlev,ylab,xlab,xlo,xhi)


#Khamkure 5 levels

# Khamkure original range
f1 <- function(x) {
  xlo = c(3,.5,.05)
  xhi = c(7,4.5,.25)
  sx = -1+2*(x-xlo)/(xhi-xlo)
  f1 = 51.97-7.42*sx[1]+13.37*sx[2]*sx[3]-18.23*sx[1]^2+15.54*sx[3]^2
}

fname1 = "AS(III)Removal"


f2 <- function(x) {
  xlo = c(3,.5,.05)
  xhi = c(7,4.5,.25)
  sx = -1+2*(x-xlo)/(xhi-xlo)
  f2 = 28.62-23.43*sx[1]+10.81*sx[2]-11.99*sx[1]*sx[2]+26.08*sx[1]^2-9.52*sx[3]^2
}
fname2 = "AS(V)Removal"

xlo = c(3,.5,.05)
xhi = c(7,4.5,.25)
ylab = "%"
xlab = NULL
xlab[1]="pH";xlab[2]="Dose";xlab[3]="AsConc"
nlev = c(5,5,5)

rgp2(f1,fname1,f2,fname2,nlev,ylab,xlab,xlo,xhi)


#Khamkure 6 levels

# Khamkure original range
f1 <- function(x) {
  xlo = c(3,.5,.05)
  xhi = c(7,4.5,.25)
  sx = -1+2*(x-xlo)/(xhi-xlo)
  f1 = 51.97-7.42*sx[1]+13.37*sx[2]*sx[3]-18.23*sx[1]^2+15.54*sx[3]^2
}

fname1 = "AS(III)Removal"


f2 <- function(x) {
  xlo = c(3,.5,.05)
  xhi = c(7,4.5,.25)
  sx = -1+2*(x-xlo)/(xhi-xlo)
  f2 = 28.62-23.43*sx[1]+10.81*sx[2]-11.99*sx[1]*sx[2]+26.08*sx[1]^2-9.52*sx[3]^2
}
fname2 = "AS(V)Removal"

xlo = c(3,.5,.05)
xhi = c(7,4.5,.25)
ylab = "%"
xlab = NULL
xlab[1]="pH";xlab[2]="Dose";xlab[3]="AsConc"
nlev = c(6,6,6)

rgp2(f1,fname1,f2,fname2,nlev,ylab,xlab,xlo,xhi)


#Khamkure 3 levels

# Khamkure original range
f1 <- function(x) {
  xlo = c(3,.5,.05)
  xhi = c(7,4.5,.25)
  sx = -1+2*(x-xlo)/(xhi-xlo)
  f1 = 51.97-7.42*sx[1]+13.37*sx[2]*sx[3]-18.23*sx[1]^2+15.54*sx[3]^2
}

fname1 = "AS(III)Removal"


f2 <- function(x) {
  xlo = c(3,.5,.05)
  xhi = c(7,4.5,.25)
  sx = -1+2*(x-xlo)/(xhi-xlo)
  f2 = 28.62-23.43*sx[1]+10.81*sx[2]-11.99*sx[1]*sx[2]+26.08*sx[1]^2-9.52*sx[3]^2
}
fname2 = "AS(V)Removal"

xlo = c(3,.5,.05)
xhi = c(7,4.5,.25)
ylab = "%"
xlab = NULL
xlab[1]="pH";xlab[2]="Dose";xlab[3]="AsConc"
nlev = c(3,3,3)

rgp2(f1,fname1,f2,fname2,nlev,ylab,xlab,xlo,xhi)

