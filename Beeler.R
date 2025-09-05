# Beeler.R
# This code implements the Beeler, Aleman and Carter (WSC 2012)
#  function
f <- function(x){
  # first, map qualitative x elements into indicators
  # vaccination
  if(x[3]==0) {C2=D2=0}
  if(x[3]==1) {C2=1;D2=0}
  if(x[3]==2) {C2=0;D2=1}
  # quarantine
  if(x[4]==0) {F2=G2=0}
  if(x[4]==1) {F2=1;G2=0}
  if(x[4]==2) {F2=0;G2=1}
  # days
  if(x[1]==0) {I2=J2=0}
  if(x[1]==-1) {I2=1;J2=0}
  if(x[1]==1) {I2=0;J2=1}
  # infectiousness
  if(x[2]==0) {L2=M2=0}
  if(x[2]==-1) {L2=1;M2=0}
  if(x[2]==1) {L2=0;M2=1}
  main = 40 + (-16*C2)+(-24*D2)+(-7*F2)+(-12*G2)+
    (-13*I2)+(14*J2)+(-26*L2)+(70*M2)
  twoWay1 = (6*D2*F2)+(7*C2*G2)+(9*D2*G2)+(8*D2*I2)+
    (-6*C2*J2)+(-9*D2*J2)+(12*C2*L2)+(17*D2*L2)+
    (-50*C2*M2)+(-65*D2*M2)
  twoWay2 = (10*G2*L2)+(-42*F2*M2)+(-56*G2*M2)+
    (10*I2*L2)+(-12*J2*L2)+(-26*I2*M2)+(101*J2*M2)
  threeWay = (12*G2*I2*M2)+(-47*F2*J2*M2)+(-57*G2*J2*M2)+
    (11*C2*I2*M2)+(16*D2*I2*M2)+(-53*C2*J2*M2)+
    (-62*D2*J2*M2)+(-8*D2*G2*L2)+(43*C2*F2*M2)+
    (51*D2*F2*M2)+(54*C2*G2*M2)+(65*D2*G2*M2)
  return(main+twoWay1+twoWay2+threeWay)
}