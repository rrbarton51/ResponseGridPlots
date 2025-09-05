# rgp2.R

# by Russell R. Barton, September 2025
#
# Given d and a function of x in R^d, construct a response-scaled design
#  plot on a factorial grid with number of levels nlev[1] ... nlev[d] and 
#  x limits xlo[1] ... xlo[d] and xhi[1] ... xhi[d]
# Variable names must be provided in a label vector xlab
# Function names in fname1, fname2
# Package needed: DOE.base

# Make sure plot window is large enough to accommodate plot legend 

rgp2 <- function(f1,fname1,f2,fname2,nlev,ylab,xlab,xlo,xhi){
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
  
  DOE = data.matrix(fac.design(factor.names=xlab[1:d],nlevels=nlev))
  
  #rescale DOE values to x and compute corresponding responses
  xDOE = matrix(0,nrow=nrow(DOE),ncol=ncol(DOE))
  for (i in 1:nrow(DOE)){
    for (j in 1:d){
      xDOE[i,j] = xlo[j] + (DOE[i,j]-1)*(xhi[j]-xlo[j])/(nlev[j]-1)
    }
  }
 colnames(xDOE) = xlab
 

  #compute response functions f1 and f2 for each row of the xDOE
  resp1 = apply(xDOE,1,f1)
  resp2 = apply(xDOE,1,f2)
  # order variables based on number of levels
  xDOE = xDOE[,order(nlev,decreasing=TRUE)]
  #attach responses as last columns
  DOE = cbind(xDOE,resp1,resp2)
  
  # determine number outer, number inner variables: choose design combinations
  #   from (1,3), (2,2), (2,3), (3,3)
  # 
  # rule follows here:
   if (d == 2){outerd = 2; innerd = 0}
   if (d == 3){outerd = 3; innerd = 0}
   if (d == 3 & min(nlev)>3){outerd = 1; innerd = 2}
   if (d == 4 & min(nlev)==2){outerd = 1; innerd = 3}
   if (d == 4 & min(nlev)>2){outerd = 2; innerd = 2}
   if (d == 5){outerd = 2; innerd = 3}
   if (d == 6){outerd = 3; innerd = 3}
  
  # construct response scaled design plot for this DOE
  return = outer(DOE,outerd,innerd)
}


#==============================================================================
#source("outer.R")  # creates outer framework plot for up to 3 factors
                    #  and repeatedly 
                    #  calls inner plot routine for different outer
                    #  plot locations
# outer.R
# Russell Barton January 2025
# given a file of +/-1 scaled predictor variables in first outerd + 
#  innerd columns plus an additional response column,
#  use the first outerd columns (a,b,c) to draw outer factorial framework 
#  this will be projected in 2-D if outerd = 3 but not in a geometric
#  projection. Instead, higher values of b lead to shifts in the a and
#  c locations (this requires shifta, shiftc) 
#  a is left-right, want to shift rear points ( when b> -1) right,
#  b, front-back, used for a and b shifts 
#  c is up-down axis, want to shift rear points (b > -1) up.
# then plot a (and shifted a) and b (and shifted b)
# key notation:
#  outer levels are called a, b, c internally but passed variable names used in plot
#  number levels each are called numalev, numblev, numclev

outer <- function(unscaledDOE,outerd,innerd){
  # make sure DOE provided is numeric - do conversion
  unscaledDOE = data.matrix(unscaledDOE)
  # make sure DOE is on [-1,1] scale - but not response 
  scalpm1 = function(x){(x-(min(x)+max(x))/2)/(.5*(max(x)-min(x)))}
  DOE = apply(unscaledDOE[,1:(outerd+innerd)],2,scalpm1)
  
  # scale response to [delta,1] this will be used for plotting circles
  #   so smallest points have nonzero diameter delta, max diameter 1 
  resp1 = unscaledDOE[,(outerd+innerd+1)]
  resp2 = unscaledDOE[,(outerd+innerd+2)]
  minresp1 = min(resp1)
  maxresp1 = max(resp1)
  minresp2 = min(resp2)
  maxresp2 = max(resp2)
  delta = 0.1
  sResp1 = delta + (resp1-minresp1)/(maxresp1-minresp1)
  sResp2 = delta + (resp2-minresp2)/(maxresp2-minresp2)
  #  attach scaled response to +/-1 scaled design matrix
  DOE = cbind(DOE,sResp1,sResp2)
  
  #  need local variable for number of levels of each of a, b, c
  numalev = length(unique(DOE[,1]))
  if(outerd > 1){numblev = length(unique(DOE[,2]))} else{
    numblev = 0
  }
  if(outerd == 3){numclev = length(unique(DOE[,3]))} else{
    numclev = 0
  }
  
  #  need total dimension of outer and inner variables
  d = outerd + innerd
  
  # outer DOE points repeated multiple times for different
  #   settings of inner DOE points - extract unique outer values
  outerDOE = as.matrix(unique(DOE[,1:outerd]))
  # convert from list to matrix
  outerDOE = apply(outerDOE,2,as.numeric)
  
  # generic outer design factor names are 'a', 'b', 'c'
  DOEabc = t(outerDOE) # each column is a design point
  DOEabcS = DOEabc # Shift place holder for outerd = 3; when b > -1, 
  # a and c values shifted for perspective plot. 
  # When outerd = 1 or 2, no shift so DOEabcS = DOEabc
  
  # outer design plotted with center at zero and scale 1
  pscale = 1 # fraction of total plot range for a subplot
  # reset scale for inner subplots, ploc will differ for each inner plot
  # spscale sets size of subplot as fraction of overall plot
  if(outerd == 1){
    spscale = .5/numalev
  } else if(outerd == 2){
    spscale = 1/(numalev+numblev+2)
  } else { # outerd == 3
    spscale = 1/(numalev+numblev+numclev)
  }
  rscale = .4 # sets size of maximum circle as fraction of subplot units
  ploc = as.vector(c(0,0))
  
  # create lines - loop on columns i and j > i 
  # take only upper part so line segments not repeated
  #  set line between pair if one coordinate differ by 2 using and rest by 0
  # set placeholder matrices for start and end points of lines
  vec0mat = NULL #matrix(0,nrow = ncol(outerDOE),ncol=1)
  vec1mat = NULL #matrix(0,nrow = ncol(outerDOE),ncol=1)
  # take only upper part so line segments not repeated
  # frame line segments have change in only one coordinate and (abs) change = 1-(-1)=2
  # allow roundoff through eps
  eps = .005
  for(i in 1:(nrow(outerDOE)-1)){
    for(j in (i+1):nrow(outerDOE)){
      if(length(which(abs(outerDOE[i,]-outerDOE[j,])>2-eps))==1 &&
         length(which(abs(outerDOE[i,]-outerDOE[j,])>-eps &
                      abs(outerDOE[i,]-outerDOE[j,])<eps))==(outerd-1))
      {
        vec0mat = cbind(vec0mat,as.matrix(outerDOE[i,]))
        vec1mat = cbind(vec1mat,as.matrix(outerDOE[j,]))
      }
    }
  }
  
  # special shift of plotted lines use xshift, yshift
  # shift values for better line locations, less if more levels
  # adjustments needed both for outer levels and inner levels
  if(innerd == 3){
    sxmult = .35
  }else { # innerd ==2
    sxmult = .3
  }
  if(innerd == 3){
    symult = .35
  }else { # innerd ==2
    symult = .3
  }
  symult = ifelse(innerd ==1,.2,.3)
  xshift = sxmult*pscale/(numalev+3)
  yshift = symult*pscale/(numclev+3)
  
  # compute shifted a and c values based on b if outerd = 3
  if(outerd == 3){
    # perform shift a and c when b > -1 for DOE points
    # best values appear to be 120% (from -1 to +1 for b) for a and 70% for c
    
    shifta = 120 # out of 200 per unit increase in b
    shiftc = 70  # out of 200 per unit increase in b
    
    vec0matS = vec0mat # place holder for shifted start points for 
    #  perspective line segments
    vec1matS = vec1mat # place holder for (possibly) shifted
    #  end points of line segments
    # shift DOE points when b > -1
    for(i in 1:ncol(DOEabc)){
      DOEabcS[2,i] = DOEabc[2,i]
      lval = DOEabc[2,i]
      DOEabcS[1,i] = DOEabc[1,i] + (2-(1-lval))*shifta/200
      DOEabcS[3,i] = DOEabc[3,i] + (2-(1-lval))*shiftc/200
    }
    
    # shift line endpoints when b>-1 by shifta, shiftc
    #  a, b coordinates for c == 1 DOE points - now points will be columns not rows
    for(i in 1:ncol(vec0mat)){
      vec0matS[2,i] = vec0mat[2,i] # not used in 2D plot
      lval = vec0mat[2,i] # used to determine shift in a and c values based on b
      vec0matS[1,i] = vec0mat[1,i] + (2-(1-lval))*shifta/200
      vec0matS[3,i] = vec0mat[3,i] + (2-(1-lval))*shiftc/200
      
      vec1matS[2,i] = vec1mat[2,i]
      lval = vec1mat[2,i] # used to determine shift in a and c values based on b
      vec1matS[1,i] = vec1mat[1,i] + (2-(1-lval))*shifta/200
      vec1matS[3,i] = vec1mat[3,i] + (2-(1-lval))*shiftc/200
    }
    # locate endpoints of design line segments using scale and location
    x0 = pscale*vec0matS[1,]+ploc[1]
    y0 = pscale*vec0matS[3,]+ploc[2]
    x1 = pscale*vec1matS[1,]+ploc[1]
    y1 = pscale*vec1matS[3,]+ploc[2]
    
    x0s = pscale*vec0matS[1,]+ploc[1]+xshift
    y0s= pscale*vec0matS[3,]+ploc[2]+yshift
    x1s = pscale*vec1matS[1,]+ploc[1]+xshift
    y1s = pscale*vec1matS[3,]+ploc[2]+yshift
    # end of code when outerd = 3 ===============
    
  } else if(outerd ==2){
    x0 = pscale*vec0mat[1,]+ploc[1]
    y0 = pscale*vec0mat[2,]+ploc[2]
    x1 = pscale*vec1mat[1,]+ploc[1]
    y1 = pscale*vec1mat[2,]+ploc[2]
    # shift values for better line locations
    x0s = x0 + xshift
    y0s = y0
    x1s = x1 + xshift
    y1s = y1
    # end of code when outerd = 2
  }  else { # if outerd == 1
    x0 = pscale*vec0mat[1,]+ploc[1]
    y0 = 0+ploc[2]
    x1 = pscale*vec1mat[1,]+ploc[1]
    y1 = 0+ploc[2]  
    # shift values for better line locations
    x0s = x0 + xshift
    y0s = y0
    x1s = x1 + xshift
    y1s = y1
  }
    # end of code when outerd = 1
  
  
  
  # plot large circles to initialize figure and permit labels to fit (these are invisible)
  # plot depends on outerd
  
  if(outerd ==1){
    symbols((pscale*DOEabcS[1,1:ncol(DOEabc)])+ploc[1],
            rep(yshift,ncol(DOEabc))+ploc[2],
            circles=rep(2.5*spscale,ncol(DOEabc)),bty="n",xaxt="n",yaxt="n",
            xlab="",ylab="", bg="white",fg="white",inches = FALSE,
            xlim=c(-1.35*pscale,1.3*pscale),ylim=c(-.1*pscale,.1*pscale),
            main = paste("Functions ",fname1,fname2),asp=1)
  } else if(outerd ==2){
    symbols(1.3*pscale*DOEabcS[1,1:ncol(DOEabc)]+ploc[1]-.05,
            1.2*pscale*DOEabcS[2,1:ncol(DOEabc)]+ploc[2]-.05,
            circles=rep(.9*pscale,ncol(DOEabc)),bty="n",xaxt="n",yaxt="n",
            xlab="",ylab="", bg="white",fg="white",,inches = FALSE,
            xlim=c(-1.7*pscale,1.4*pscale),ylim=c(-1.4*pscale,1.4*pscale),
            main = paste("Functions ",fname1,fname2),asp=1)
  } else{ # outerd == 3
    symbols(1.1*pscale*DOEabcS[1,1:ncol(DOEabc)]+ploc[1]-.05,
            1.1*pscale*DOEabcS[3,1:ncol(DOEabc)]+ploc[2]-.05,
            circles=rep(.9*pscale,ncol(DOEabc)),bty="n",xaxt="n",yaxt="n",
            xlab="",ylab="", bg="white",fg="white",,inches = FALSE,
            xlim=c(-2.5*pscale,3.5*pscale),ylim=c(-1.5*pscale,2.5*pscale),
            main = paste("Functions ",fname1,fname2),asp=1)
  }
  
  # add lines of cube - background thin and black if outerd = 3
  if(outerd == 3){
    segments(x0s,y0s,x1s,y1s,col="black")
    }else{segments(x0s,y0s,x1s,y1s,col="red",lwd=2)
  }
  
  # if innerd > 0, plot large circles again to clear DOE line segments 
  #   near vertices for sub-figures -- plot depends on outerd
  if(innerd>0){
    if(outerd ==1){
      symbols(pscale*DOEabcS[1,1:ncol(DOEabc)]+ploc[1]+xshift,
              rep(yshift,ncol(DOEabc))+ploc[2],
              circles=rep(.2*pscale,ncol(DOEabc)),bty="n",xaxt="n",yaxt="n",
              xlab="",ylab="", bg="white",fg="white",inches = FALSE,
              add=TRUE, main = paste("Function ",fname))
    } else if(outerd ==2){
      symbols(pscale*DOEabcS[1,1:ncol(DOEabc)]+ploc[1]+.03,
              pscale*DOEabcS[2,1:ncol(DOEabc)]+ploc[2],
              circles=rep(.22*pscale,ncol(DOEabc)),bty="n",xaxt="n",yaxt="n",
              xlab="",ylab="", bg="white",fg="white",inches = FALSE,
              add=TRUE, main = paste("Function ",fname))
    } else{ # outerd==3
      symbols(pscale*DOEabcS[1,1:ncol(DOEabc)]+ploc[1]+xshift,
              pscale*DOEabcS[3,1:ncol(DOEabc)]+ploc[2]+yshift,
              circles=rep(.15*pscale,ncol(DOEabc)),bty="n",xaxt="n",yaxt="n",
              xlab="",ylab="", bg="white",fg="white",inches = FALSE,
              add=TRUE, main = paste("Function ",fname))
    }
  }
  
  # if outerd = 3, plot foreground lines in red and slightly thicker than
  #  background lines which are black (no background lines for outerd = 1, 2)
  # needs numalev, numblev, numclev - number of levels for each outer factor - set above
  # NOTE: since the number of levels may result in ldel being a 
  #  nonterminating decimal.strict equality cannot be checked. 
  #  Instead, determine if the absolute difference is less than eps.
  if(outerd ==3){

    #  front lines (horizontal)
    ldel = 2/(numclev-1)
    #for(i in 1:numclev){
      segments(x0s[(abs(1+x0)<eps & abs(1-x1)<eps) |
                     (abs(1-x0)<eps & abs(1+x1)<eps)],
               y0s[(abs(1+x0)<eps & abs(1-x1)<eps) |
                     (abs(1-x0)<eps & abs(1+x1)<eps)],
               x1s[(abs(1+x0)<eps & abs(1-x1)<eps) |
                     (abs(1-x0)<eps & abs(1+x1)<eps)],
               y1s[(abs(1+x0)<eps & abs(1-x1)<eps) |
                     (abs(1-x0)<eps & abs(1+x1)<eps)],col="red",lwd=2)
    #}
    
    #  front lines (vertical)
    ldel = 2/(numalev-1)
    for(i in 1:numalev){
      segments(x0s[abs(x0+1-ldel*(i-1))<eps & abs(x1+1-ldel*(i-1))<eps],
               y0s[abs(x0+1-ldel*(i-1))<eps & abs(x1+1-ldel*(i-1))<eps],
               x1s[abs(x0+1-ldel*(i-1))<eps & abs(x1+1-ldel*(i-1))<eps],
               y1s[abs(x0+1-ldel*(i-1))<eps & abs(x1+1-ldel*(i-1))<eps],col="red",lwd=2)
    }
    
    #  top lines (horizontal) - don't repeat front top horizontal
    ldel = 2/(numblev-1)
    for(i in 2:numblev){
      segments(x0s[abs(y0-(1+ldel*(i-1)*shiftc/200))<eps & abs(y1-(1+ldel*(i-1)*shiftc/200))<eps],
               y0s[abs(y0-(1+ldel*(i-1)*shiftc/200))<eps & abs(y1-(1+ldel*(i-1)*shiftc/200))<eps],
               x1s[abs(y0-(1+ldel*(i-1)*shiftc/200))<eps & abs(y1-(1+ldel*(i-1)*shiftc/200))<eps],
               y1s[abs(y0-(1+ldel*(i-1)*shiftc/200))<eps & abs(y1-(1+ldel*(i-1)*shiftc/200))<eps],col="red",lwd=2)
    }
    
    #  top lines (vertical) 
    ldel = 2/(numalev-1)
    for(i in 1:numalev){
      segments(x0s[(abs(x0+1-ldel*(i-1))<eps | abs(x1+1-ldel*(i-1))<eps)
                   & (abs(y0-(1+shiftc/100))<eps | abs(y1-(1+shiftc/100))<eps)],
               y0s[(abs(x0+1-ldel*(i-1))<eps | abs(x1+1-ldel*(i-1))<eps)
                   & (abs(y0-(1+shiftc/100))<eps | abs(y1-(1+shiftc/100))<eps)],
               x1s[(abs(x0+1-ldel*(i-1))<eps | abs(x1+1-ldel*(i-1))<eps)
                   & (abs(y0-(1+shiftc/100))<eps | abs(y1-(1+shiftc/100))<eps)],
               y1s[(abs(x0+1-ldel*(i-1))<eps | abs(x1+1-ldel*(i-1))<eps)
                   & (abs(y0-(1+shiftc/100))<eps | abs(y1-(1+shiftc/100))<eps)],
               col="red",lwd=2)
    }
    
    #  right side lines (vertical) - don't repeat front right vertical
    ldel = 2/(numblev-1)
    for(i in 2:numblev){
      segments(x0s[abs(x0-(1+ldel*(i-1)*shifta/200))<eps & abs(x1-(1+ldel*(i-1)*shifta/200))<eps],
               y0s[abs(x0-(1+ldel*(i-1)*shifta/200))<eps & abs(x1-(1+ldel*(i-1)*shifta/200))<eps],
               x1s[abs(x0-(1+ldel*(i-1)*shifta/200))<eps & abs(x1-(1+ldel*(i-1)*shifta/200))<eps],
               y1s[abs(x0-(1+ldel*(i-1)*shifta/200))<eps & abs(x1-(1+ldel*(i-1)*shifta/200))<eps],col="red",lwd=2)
    }
    
    # right side lines (horizontal) - don't repeat top line
    ldel = 2/(numclev-1)
    for(i in 1:(numclev-1)){
      segments(x0s[(abs(1-x0)<eps | abs(1-x1)<eps) & 
                     (abs(y0-(-1+ldel*(i-1)+shiftc/100))<eps | 
                                         abs(y1-(-1+ldel*(i-1)+shiftc/100))<eps)],
               y0s[(abs(1-x0)<eps | abs(1-x1)<eps) & 
                     (abs(y0-(-1+ldel*(i-1)+shiftc/100))<eps | 
                        abs(y1-(-1+ldel*(i-1)+shiftc/100))<eps)],
               x1s[(abs(1-x0)<eps | abs(1-x1)<eps) & 
                     (abs(y0-(-1+ldel*(i-1)+shiftc/100))<eps | 
                        abs(y1-(-1+ldel*(i-1)+shiftc/100))<eps)],
               y1s[(abs(1-x0)<eps | abs(1-x1)<eps) & 
                     (abs(y0-(-1+ldel*(i-1)+shiftc/100))<eps | 
                        abs(y1-(-1+ldel*(i-1)+shiftc/100))<eps)],col="red",lwd=2)
    }
    # end of outerd = 3 red highlighting
  }
  
  # label cube axes with generic names x1, x2, etc.
  plotnames = c("x1","x2","x3","x4","x5","x6")
  epsx = max(c(x0,x1))-min(c(x0,x1))
  epsy = max(c(y0,y1))-min(c(y0,y1))
  tscale = pscale*1.0
  
  if(outerd==1){ # shift x axis label below spscaled plots
    text((min(c(x0,x1))+max(c(x0,x1)))/2,-2.5*spscale,plotnames[1],cex=tscale)
  } else if(outerd==2){ #  for outerd =2 or 3 epsy not zero and can be used for shift of x axis label
    text((min(c(x0,x1))+max(c(x0,x1)))/2,min(c(y0,y1))-.2*epsy,plotnames[1],cex=tscale)
  } else { # outerd==3
    text(((min(c(x0,x1))+max(c(x0,x1)))/2)-.15*epsx,min(c(y0,y1))-.18*epsy,plotnames[1],cex=tscale)
  }
  if(outerd == 2){text(max(c(x0,x1))+.2*epsx,(min(c(y0,y1))+max(c(y0,y1)))/2,plotnames[2],srt=90,cex=tscale)}
  if(outerd == 3){text(max(c(x0,x1))-.3,min(c(y0,y1))+.05,plotnames[2],srt=40,cex=tscale)}
  if(outerd == 3){text(max(c(x0,x1))+.15*epsx,+min(c(y0,y1))+.65*epsy,plotnames[3],srt=90,cex=tscale)}
  # set legend text using colnames
  legtext = NULL
  for (ivar in 1:d){
    legtext[ivar] = paste(plotnames[ivar]," = ",colnames(DOE)[ivar])
  }
    legtext[d+1] = paste("min ",ylab," values =",signif(minresp1,digits=3),signif(minresp2,digits=3))
    legtext[d+2] = paste("max ",ylab," values =",signif(maxresp1,digits=3),signif(maxresp2,digits=3))
  if(outerd==1){
    legend(x=-1.2*pscale,y=.4*pscale,title = "Key", legend=legtext, 
           col=c(rep("white",d),rep("grey30",2)), pch=c(rep(19,d),20,19), 
           box.col="white",cex=.8)
  }else if(outerd==2){
    legend(x=-2*pscale,y=.5*pscale,title = "Key", legend=legtext, 
           col=c(rep("white",d),rep("grey30",2)), pch=c(rep(19,d),20,19), 
           box.col="white",cex=.8)
    }else{ # outerd==3
    legend(x=-1.8*pscale,y=2.6*pscale, title = "Key", inset=c(0,0), legend=legtext, 
           col=c(rep("white",d),rep("grey30",2)), pch=c(rep(19,d),20,19),
           box.col="white",cex=.8)
  }
    # add inner subplots
    
    # ploc is the 2-D position of the subplot, shifted if b > -1 
    #for(irow in 1:1){
    for(irow in 1:nrow(outerDOE)){
      outerRow = outerDOE[irow,]
      if(outerd ==3){
        ploc = twoD(outerRow, shifta, shiftc)
      } else if(outerd == 2) {
        ploc = outerRow
      } else if(outerd ==1){
        ploc = c(outerRow,0)
      }
      # if innerd = 0 just plot the circles at ploc
      if(innerd==0 & outerd<3){
        symbols(ploc[1]+xshift,ploc[2],
                circles=.5*min((2/numalev),(2/numblev))*rscale*sqrt(sResp1[irow]),bty="n",xaxt="n",yaxt="n",
                xlab="",ylab="", bg="grey30",fg="grey30", inches=FALSE, add=TRUE)
        symbols(ploc[1]+xshift+
                  .5*min((2/numalev),(2/numblev))*rscale*sqrt(sResp1[irow])+.5*min((2/numalev),(2/numblev))*rscale*sqrt(sResp2[irow]),
                  ploc[2],
                  circles=.5*min((2/numalev),(2/numblev))*rscale*sqrt(sResp2[irow]),bty="n",xaxt="n",yaxt="n",
                  xlab="",ylab="", bg="grey30",fg="grey30", inches=FALSE, add=TRUE)
      } else if(innerd==0 & outerd==3){
        symbols(ploc[1]+xshift,ploc[2]+yshift,
                  circles=.5*min((2/numalev),(2/numblev),(2/numclev))*rscale*sqrt(sResp1[irow]),bty="n",xaxt="n",yaxt="n",
                  xlab="",ylab="", bg="grey30",fg="grey30", inches=FALSE, add=TRUE)
        symbols(ploc[1]+xshift+
                  .5*min((2/numalev),(2/numblev))*rscale*sqrt(sResp1[irow])+.5*min((2/numalev),(2/numblev))*rscale*sqrt(sResp2[irow]),
                  ploc[2]+yshift,
                  circles=.5*min((2/numalev),(2/numblev),(2/numclev))*rscale*sqrt(sResp2[irow]),bty="n",xaxt="n",yaxt="n",
                  xlab="",ylab="", bg="grey30",fg="grey30", inches=FALSE, add=TRUE)
      }
      # if innerd > 0
      if(innerd>0){
        # call inner with ploc and proper subset of DOE[,outerd+1:d]
        if(outerd==1){
          innerDOE = DOE[DOE[,1:outerd] == outerRow,((outerd+1):(d+2))]
        } else {
          innerDOE = DOE[apply(DOE[,1:outerd], 1, 
                               function(y) return(all(y == outerRow))),((outerd+1):(d+2))]
        }
        cgray=ifelse(outerd==3,gray((1+outerRow[2])/3.8),"grey30") # set background to grey
        return2 = inner(innerDOE,outerd,innerd,ploc,spscale,rscale,cgray)
      } # end of if innerd > 0
      
    } # end of loop for each outerd point
}
  
# end of function outer
#==============================================================================

# twoD - function to place 3D points in two-D representation using shifta and shiftc
twoD <- function(vec3D, shifta, shiftc){
  twoDvec = rep(0,2)
  twoDvec[1] = vec3D[1] + (vec3D[2]+1)*shifta/200
  twoDvec[2] = vec3D[3] + (vec3D[2]+1)*shiftc/200
  return(twoDvec)
}
# end of function twoD
#==============================================================================
# inner.R
# Russell Barton January 2025
# given a design in two or three columns and response in last, 
#  draw points projected in 2-D (if 3 variables, i.e., innerd = 3)
#  we name the inner variables x, y, and possibly z if 3 variables
# if innerd = 3 we need projection, this requires shiftx, shiftz 
#  x is left-right, want to shift rearward points (y>-1) right,
#  y, front-back, used for x and z shifts 
#  z is up-down axis, want to shift rearward points (y>-1) up.
# if y == 1 sx = x + 2*shiftx/200 (default shiftx = 120)
# if y == 1 sz = z + 2*shiftz/200 (default shiftz = 90)
# then plot sx and sz

inner <- function(innerDOE,outerd,innerd,ploc,spscale,rscale,cgray){
  # spscale is the size of the subplot, rscale reduces max circle size to fraction
  # Assume DOE on [-1,1] scale and sResp scaled to [small nonzero, 1]
  DOE = data.matrix(innerDOE)
  sResp1 = DOE[,(innerd+1)]
  sResp2 = DOE[,(innerd+2)]
  DOExyz = t(DOE) # each column is a design point
  DOExyzS = DOExyz # Shift place holder for innerd = 3; when y > -1, 
  # x and z values shifted for perspective plot. 
  # When innerd = 1 or 2, no shift so DOExyzS = DOExyz
  
  # create lines - loop on columns i and j > i 
  # take only upper part so line segments not repeated
  #  set line between pair if one coordinate differ by 2 using and rest by 0
  # set placeholder matrices for start and end points of lines
  vec0mat = NULL #matrix(0,nrow = ncol(DOE),ncol=1)
  vec1mat = NULL #matrix(0,nrow = ncol(DOE),ncol=1)
  # extract just the DOE points (DOEpts) without response for distance
  innerDOE = DOE[,1:innerd]
  # frame line segments have change in only one coordinate and (abs) change = 1-(-1)=2
  for(i in 1:(nrow(innerDOE)-1)){
    for(j in (i+1):nrow(innerDOE)){
      if(length(which(abs(innerDOE[i,]-innerDOE[j,])==2))==1 &
         length(which(abs(innerDOE[i,]-innerDOE[j,])==0))==(innerd-1))
      {
        vec0mat = cbind(vec0mat,as.matrix(innerDOE[i,]))
        vec1mat = cbind(vec1mat,as.matrix(innerDOE[j,]))
      }
    }
  }
  
  # if innerd = 3, perform shift x and z when y > -1 for DOE points
  # best values appear to be 120% (from -1 to +1 for y) for x and 70% for z
  if(innerd == 3){
    shiftx = 120 # out of 200 per unit y
    shiftz = 90  # out of 200 per unit y
    vec0matS = vec0mat # place holder for shifted points for perspective line segments
    vec1matS = vec1mat
    
    # shift DOE points when y > -1
    for(i in 1:ncol(DOExyz)){
      DOExyzS[2,i] = DOExyz[2,i]
      lval = DOExyz[2,i]
      DOExyzS[1,i] = DOExyz[1,i] + (2-(1-lval))*shiftx/200
      DOExyzS[3,i] = DOExyz[3,i] + (2-(1-lval))*shiftz/200
    }
    
    
    # shift line endpoints when y > -1 by shiftx, shiftz factors
    #  points will be columns not rows
    for(i in 1:ncol(vec0mat)){
      vec0matS[2,i] = vec0mat[2,i] # not used in 2D plot
      lval = vec0mat[2,i] # used to determine shift in x and z values based on y
      vec0matS[1,i] = vec0mat[1,i] + (2-(1-lval))*shiftx/200
      vec0matS[3,i] = vec0mat[3,i] + (2-(1-lval))*shiftz/200
      
      vec1matS[2,i] = vec1mat[2,i]
      lval = vec1mat[2,i] # used to determine shift in x and z values based on y
      vec1matS[1,i] = vec1mat[1,i] + (2-(1-lval))*shiftx/200
      vec1matS[3,i] = vec1mat[3,i] + (2-(1-lval))*shiftz/200
    }
    
    # locate endpoints of design line segments using scale and location
    x0 = spscale*vec0matS[1,]+ploc[1]
    y0 = spscale*vec0matS[3,]+ploc[2]
    x1 = spscale*vec1matS[1,]+ploc[1]
    y1 = spscale*vec1matS[3,]+ploc[2]
    
    # end of shifting code for innerd = 3   ==============
    
  } else if(innerd ==2){
    x0 = spscale*vec0mat[1,]+ploc[1]
    y0 = spscale*vec0mat[2,]+ploc[2]
    x1 = spscale*vec1mat[1,]+ploc[1]
    y1 = spscale*vec1mat[2,]+ploc[2]
    # end of code when innerd = 2
  }  else {
    x0 = spscale*vec0mat[1,]+ploc[1]
    y0 = 0+ploc[2]
    x1 = spscale*vec1mat[1,]+ploc[1]
    y1 = 0+ploc[2]  
  }
  # end of code when innerd = 1
  
  
  # plot scaled response data (column 4) projected to x and z if innerd = 3
  #  and scaled and located - scale sResp to psResp to size circles for subplot
  if(innerd==1) {
    psResp1=.8*spscale*rscale*sqrt(sResp1);psResp2=.8*spscale*rscale*sqrt(sResp2)
  }
  if(innerd==2) {
    psResp1=.8*spscale*rscale*sqrt(sResp1);psResp2=.8*spscale*rscale*sqrt(sResp2)
  }
  if(innerd==3) {
    psResp1=spscale*rscale*sqrt(sResp1);psResp2=spscale*rscale*sqrt(sResp2)
  }
 
  if(innerd ==1){
    symbols(spscale*DOExyzS[1,1:ncol(DOExyz)]+ploc[1],rep(0,ncol(DOExyz))+ploc[2],
            circles=.5*psResp1,bty="n",xaxt="n",yaxt="n",
            xlab="",ylab="", bg=cgray,fg=cgray, inches=FALSE, add=TRUE)
    symbols(spscale*DOExyzS[1,1:ncol(DOExyz)]+ploc[1]+.5*(psResp1+psResp2),rep(0,ncol(DOExyz))+ploc[2],
            circles=.5*psResp2,bty="n",xaxt="n",yaxt="n",
            xlab="",ylab="", bg=cgray,fg=cgray, inches=FALSE, add=TRUE)
  } else if(innerd ==2){
    symbols(spscale*DOExyzS[1,1:ncol(DOExyz)]+ploc[1],spscale*DOExyzS[2,1:ncol(DOExyz)]+ploc[2],
            circles=.5*psResp1,bty="n",xaxt="n",yaxt="n",
            xlab="",ylab="", bg=cgray,fg=cgray, inches=FALSE, add=TRUE)
    symbols(spscale*DOExyzS[1,1:ncol(DOExyz)]+ploc[1]+.5*(psResp1+psResp2),spscale*DOExyzS[2,1:ncol(DOExyz)]+ploc[2],
            circles=.5*psResp2,bty="n",xaxt="n",yaxt="n",
            xlab="",ylab="", bg=cgray,fg=cgray, inches=FALSE, add=TRUE)
  } else{
    symbols(spscale*DOExyzS[1,1:ncol(DOExyz)]+ploc[1],spscale*DOExyzS[3,1:ncol(DOExyz)]+ploc[2],
            circles=.5*psResp1,bty="n",xaxt="n",yaxt="n",
            xlab="",ylab="", bg=cgray,fg=cgray, inches=FALSE, add=TRUE)
    symbols(spscale*DOExyzS[1,1:ncol(DOExyz)]+ploc[1]+.5*(psResp1+psResp2),spscale*DOExyzS[3,1:ncol(DOExyz)]+ploc[2],
            circles=.5*psResp2,bty="n",xaxt="n",yaxt="n",
            xlab="",ylab="", bg=cgray,fg=cgray, inches=FALSE, add=TRUE)
  }
  
  
  # add lines of design
  segments(x0,y0,x1,y1,col=cgray)
  
  
  # label cube axes
  epsx = max(x1)-min(x0)
  epsy = max(y1)-min(y0)
  tscale = spscale*1.6+.4
  plotnames = c("x1","x2","x3","x4","x5","x6")
  if(innerd == 1 | innerd == 2){
    text(((min(x0)+max(x0))/2)-.01*epsx,min(y0)-.3*epsy,plotnames[outerd+1],cex=1.5*tscale)
  }
  if(innerd == 2){
    text(max(x1)+.3*epsx,+min(y0)+.5*epsy,plotnames[outerd+2],srt=90,cex=1.5*tscale)
  }
  if(innerd == 3){
    text(min(x0)+.28*epsx,min(y0)-.2*epsy,plotnames[outerd+1],cex=1.2*tscale)
    text(max(x1)-.05*epsx,min(y0)-.035*epsy,plotnames[outerd+2],srt=40,cex=1.2*tscale)
    text(max(x1)+.3*epsx,+min(y0)+.6*epsy,plotnames[outerd+3],srt=90,cex=1.2*tscale)
  }
  
}
# end of function inner

