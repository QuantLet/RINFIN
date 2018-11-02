[<img src="https://github.com/QuantLet/Styleguide-and-FAQ/blob/master/pictures/banner.png" width="888" alt="Visit QuantNet">](http://quantlet.de/)

## [<img src="https://github.com/QuantLet/Styleguide-and-FAQ/blob/master/pictures/qloqo.png" alt="Visit QuantNet">](http://quantlet.de/) **RINFIN** [<img src="https://github.com/QuantLet/Styleguide-and-FAQ/blob/master/pictures/QN2.png" width="60" alt="Visit QuantNet 2.0">](http://quantlet.de/)

```yaml


Name of Quantlet:  RINFIN

 

Published in:      "RESIDUAL'S INFLUENCE INDEX (RINFIN), BAD LEVERAGE AND UNMASKING IN HIGH DIMENSIONAL L2-REGRESSION"

  

Description:       "RESIDUAL'S INFLUENCE INDEX (RINFIN), BAD LEVERAGE AND UNMASKING IN HIGH DIMENSIONAL L2-REGRESSION"

 

Keywords:          "Big Data, Data Science, Influence Function, Leverage, Masking, Residual's Influence Index (RINFIN)"



See also :         RINFIN



Author:            Yannis G. Yatracos

  

Submitted:         Tue, Oct 30 2018 by Raphael Reule

  

Datafile:          R Scripts

```

### R Code
```r

# FUNCTIONS FOR SIMULATIONS-TABLES 1 AND  2
# Note to the reader:   For Tables 1, 2,  the function # 1 creates each time a column for the data from F (Normal), with ???mesos??? (i.e. mean)=0,   ???s??? the parameter giving variance-covariance matrix,  dimension d (via J in the function as explained below), d=10,30,50,70, 90 when the sample size n=100, the number of simulations N=100.  The contaminated components follow also a normal distribution with contaminated mean \mu denoted in the function ???mesoscont??? and Covariance matrix the identity. The parameter p is the proportion of contaminated first  d*p (integer) coordinates.  When p=1, there is full contamination.
# DURATION TO OBTAIN RESULTS FOR EACH COLUMN OF THE TABLE WITH d=10,30,50,70, 90 in R-STUDIO: 10 hours; for d>90 it takes much longer.
# The function # 1) is used; the other functions are in the background, i.e. inside function 1, then inside 2 etc. 
# Parameters: n=sample size, N=number of simulations, d=data???s dimension, s=parameter giving variance-covariance matrix, mesos=mean of data,  mesoscont=contaminated mean, parameter J helps creating the dimension of the data noted d (unlike p in thee paper), =5 gives dimension d=10*(2*i-1), i=1-J,to obtain results for d=10,30,50,70,90; p is the proportion of contaminated first  d*p (integer) coordinates.   
# 
# THE FUNCTIONS
# In function # 1 you can see that the functions in it are identified  ??? # 2 FUNCTION, # 3 FUNCTION etc???.
# 
# 1)
rinsq11EX2.1ALLd (J=5,n=100,s=.5,  mesos=0,mesoscont=.5,N=100, p=1) 
{
  #p=proportion of the same contaminated coordinates
  results<-matrix(rep(0,9*J), ncol=9)
  print(results)
  for (j in 1:J)
  {
    #2 FUNCTION
    results[j,]<-        rinsq11EX2.1(n,d=10*(2*j-1),s, mesos, mesoscont, N,p) 
  } 
  cat("RESULTS","\n")
  cat("PERCENTAGE OF CONTAMINATED COORDINATES Gamma=",p, "\n")
  cat("avmisclasrinf,Contaminated Coord, d, s, mesos, mesoscont, N, n, .1*n", "\n")
  print(results)
}

#2)
rinsq11EX2.1 (n=100,d=20,s=.5,  mesos=0,mesoscont=1,N=100, p=.2) 
{
  
  { library("MASS")   
    # misproprinfstar<-rep(0,N)
    misproprinf<-rep(0,N)
    for (j in 1:N)
    {
      # 3 FUNCTION 
      mis<-rinsq11EX2.2gendatapartcont(n,d,s,mesos,mesoscont,p)
      cat("J=",j,"OUT OF",N,"ITERATIONS", "MISSCLASS=", mis,"\n")
      # NOT NEEDED IN VIEW OF AVMISCLASRINF misproprinfstar[j]<-mis[1]
      misproprinf[j]<-mis
    }
    avmisclasrinf<-mean(misproprinf)  
    
    cat("AVERAGE MISCLAS PROP RINF=",avmisclasrinf,"\n")
    cat("PARAMETERS: n=",n, "d=",d,"MEAN F=", mesos,"MEAN G=", 
        mesoscont,"SD=",s, "NUMBER OF SIMULATIONS=",N, "CONTAMINATED
        OBSERVATIONS IN THE BEGINNING OF LIST 
        ARE TEN PERCENT=", .1*n,"\n")
  }
  cat("avmisclasrinf, contaminated coordinates,d, s, mesos, mesoscont, N, n, .1*n", "\n")
  c(avmisclasrinf,d*p,d,s, mesos, mesoscont,N,n, .1*n)
  }

#3)
rinsq11EX2.2gendatapartcont (n=50,d=10,s=.5,mesos=0,mesoscont=50, p=.1) 
{
  
  {
    #################################################ADDITIONS
    bsmall<-c(1.5,.5,0,1,0,0,1.5,0,0,0,1) 
    if  (d>=11) 
    {b<-c(1.5,.5,0,1,0,0,1.5,0,0,0,1,rep(0,d-11))}
    else {b<-bsmall[1:d]}
    ##################################################
    # 4 FUNCTION 
    Sigma<-gensigma(d,s)
    x<-mvrnorm(.9*n, Sigma, mu=rep(mesos,d))
    e1<-rnorm(.9*n,mean=0, sd=1)
    cat("OBSERVATIONS OBTAINED","\n")
    #print(x)
    cat("ERRORS OBTAINED","\n")
    #print(e1)
    y<-x%*%b+e1
    #xcont=xcontaminated
    xcontvoitheia<-mvrnorm(.1*n,Sigma=diag(1,nrow=d),mu=rep(mesoscont,d))
    cat("CONTAMINATED OBSERVATIONS OBTAINED 10 PERCENT","\n")
    #y<-x%*%b+e1
    xvoitheia<-mvrnorm(.1*n, Sigma, mu=rep(mesos,d))
    e2<-rnorm(.1*n,mean=0, sd=1)
    D<-p*d
    if (D<d)
    {xcont<-cbind(xcontvoitheia[,1:D], xvoitheia[,(D+1):d])}
    else {xcont<-xcontvoitheia[,1:D]}
    cat("# HERE xcont is indeed partially contaminated only")
    ycont<-xcont%*%b+e2
    #cat("X TIMES b", "\n")
    yall<-c(ycont,y)
    xall<-rbind(xcont,x)
    cat("ALL CASES PUT TOGETHER FIRST CONTAMINATED")
    #print(yall)  #print(xall)
    #mis<-simrinfgen(yall,xall)
    # FUNCTION 5
    mis<-rinsq11EX1.3(yall,xall)
    cat("CONTAMINATED ARE THE FIRST OBSERVATIONS","\n")
    cat("d=",d, "n=", n, "(mcont-m)/s=","\n",
        (mesoscont-mesos)/s,"10 PERCENT OF DATA IS=",.1*n,"\n")
    
  }
  mis
}

#4)
gensigma (d,s) 
{
  {SIGMA<-matrix(rep(0,d^2),nrow=d)
  #cat("SIGMA", "\n")
  #print(SIGMA)
  for (i in 1:d)
  {
    for (j in i:d)
    {
      SIGMA[i,j]<-s^{j-i}
      SIGMA[j,i]<-s^{j-i}
    }        
  }
  #cat("SIGMA WITH VALUES", "\n")
  #print(SIGMA)
  }
  SIGMA
  #CHECKED VALUES TO BE ON THE SAFE SIDE
  #cat("VAR-COV MATRIX=", "\n")
  #print(SIGMA)
} 

#5) 
rinsq11EX1.3 (y,x) 
{
  
  {  
    #STARTING { NOT TO BE CONFUSED BY IT
    n<-nrow(x)
    d<-ncol(x)
    cat("NUMBER OF CASES n=", n, "\n")
    RINFIN<-rep(0,n)
    for (i in 1:n)
    { #STARTS ONE LOOP FOR CALCULATING ALL THE 
      #TAKING THE i-th CASE OUT, 
      #DOING LEAST SQUARES, FIND RESIDUAL
      # OF THE i-th CASE WITH RESPECT TO 
      #THE REGRESSION COEFFICIENT b OBTAINED
      #BELOW
      cat("TAKE OUT CASE I=", i, "\n")
      AVoitheiax<-x[-i,]
      AVoitheiay<-y[-i]
      results<-lm(AVoitheiay ~AVoitheiax) 
      b<- results$coefficients
      cat("b=", b, "\n")
      residual<-y[i]-sum(b*c(1,x[i,]))
      cat("Residual=", residual, "\n")
      
      # CALCULATING THE ESTIMATE OF E MATRIX
      #AND ITS INVERSE TO FIND INFLUENCE AND RINFIN
      # FUNCTION 6    
      INVERSE<-rin10prodmatrix(AVoitheiax)  
      INFLUENCEFUNCTIONS<-rep(0,d+1)
      cat("INFLUENCE FUNCTIONS=", INFLUENCEFUNCTIONS, "\n")
      residualtimesvector<-residual*c(1,x[i,])
      INFLUENCEFUNCTIONS<-INVERSE%*%residualtimesvector
      cat("INFLUENCE FUNCTIONS FOR CASE i=",i, "\n" )
      print(INFLUENCEFUNCTIONS)
      
      #CALCULATING INFLUENCE IN A WAY THAT HELPS CALCULATING
      #DERIVATIVES AND RINFIN
      
      xvoith<-c(1,x[i,])
      INFLUENCEFUNCTIONSMANUAL<-rep(0,d+1)
      COMPONENT<-rep(0,d+1)
      for (j in 1:(d+1))
      {COMPONENT[j]<-sum (INVERSE[j,]*xvoith)
      INFLUENCEFUNCTIONSMANUAL[j]<- residual*COMPONENT[j]
      #residual*sum (INVERSE[j,]*xvoith)
      }
      #cat("INFLUENCE FUNCTIONS MANUAL FOR CASE i=",i, "\n" )
      #print(INFLUENCEFUNCTIONSMANUAL)
      #MATRIX OF DERIVATIVES HAS (d+1) COLUMNS AS THE NUMBER OF 
      #INFLUENCE FUNCTIONS AND THERE ARE d PARTIAL DERIVATIVES
      DERIVATIVESOFINFLUENCEFUNCTIONS<-matrix(rep(0,d*(d+1)),ncol=d+1) 
      for (m in 1:d)
      {
        for (j in 1:(d+1))
        {DERIVATIVESOFINFLUENCEFUNCTIONS[m,j]<--b[m+1]*COMPONENT[j]
        +residual*INVERSE[j,m+1]
        
        }
        
        #cat("DERIVATIVES OF INFLUENCE FUNCTIONS")
        #print(DERIVATIVESOFINFLUENCEFUNCTIONS)
        
        
      }
      cat("DERIVATIVES OF INFLUENCE FUNCTIONS", "\n")
      print(DERIVATIVESOFINFLUENCEFUNCTIONS)
      INFLUENCEITHCOORDINATEOFX<-rep(0,d)
      for(k in 1:d)
      {
        INFLUENCEITHCOORDINATEOFX[k]<-INFLUENCEFUNCTIONS[k+1]+
          sum(DERIVATIVESOFINFLUENCEFUNCTIONS[k,]*xvoith)
        # CLOSES THE "FOR" 
        cat("INFLUENCEITHCOORDINATEOFX[", k,"]", "\n") 
        print(INFLUENCEITHCOORDINATEOFX[k])
      }    
      
      RINFIN[i]    <-sum((INFLUENCEITHCOORDINATEOFX)^2/n)
      cat("RINFINI")
      print(RINFIN[i])
      #cat("INFLUENCEITHCOORDINATEOFX[k]")
      #print(RINFIN[i])
      
    }
    
    #cat("CASE")
    #print(cbind(1:n, RINFIN))
    colorder<-c(1:n)
    orinf<-order(RINFIN)
    #cat("THESE ARE THE ORDERED RINFIN SCORES", "\n")
    #print(cbind(colorder[orinf], RINFIN[orinf]))
    
    #####MAYBE MEAT STARTS HERE 
    
    #o2mean<-order(RINFIN)
    
    #orinf<-order(ABSINDEXRINF)
    #colorder<-c(1:nrows)
    
    
    n<-length(y)
    #cat("RINFIN* 10% OF n", "\n")
    print(colorder[orinf][-c(1:(.9*n))])
    highestrinfstar<-colorder[orinf][-c(1:(.9*n))]
    #cat("HIGHEST RINFIN*","\n")
    #print(highestrinfstar)
    #cat("WHICH ARE SMALLER THAN .1*n","\n")
    #print(which(highestrinfstar <=.1*n))
    lengthrinfstar<-length(which(highestrinfstar <=.1*n))
    proprinfstar<-lengthrinfstar/(.1*n)
    #cat("SUCCESS PROPORTION RINF*=",proprinfstar,"\n")
    
    nrows<-n
    cat("\n", "RINFIN RESULTS: CASE  AND TIMES BY 1/n, n=", nrows,"\n")
    print(cbind(colorder[orinf],
                RINFIN[orinf]))
    
    ###THE MEAT BELOW#######################################
    ###########MEAT##############################################
    cat("RINFIN 10% OF n", "\n")
    print(colorder[orinf][-c(1:(.9*n))])
    highestrinf<-colorder[orinf][-c(1:(.9*n))]
    cat("HIGHEST RINFIN","\n")
    print(highestrinf)
    cat("WHICH ARE SMALLER THAN .1*n","\n")
    print(which(highestrinf <=.1*n))
    lengthrinf<-length(which(highestrinf <=.1*n))
    proprinf<-lengthrinf/(.1*n)
    cat("SUCCESS PROPORTION RINF=",proprinf,"\n")
    
    #misclassproprinfstar<-1-proprinfstar
    misclassproprinf<-1-proprinf
    
    cat("MISCLASS PROPORTION RINF=",misclassproprinf,"\n")
    
  }
  mis<-misclassproprinf
  
  
}

#6) 
rin10prodmatrix  (A) 
{
  {   
    print(A)
    n<-nrow(A)
    d<-ncol(A)
    #A<-matrix(A,ncol=d)
    #print(A)
    ESTIMATESOFE<-matrix(rep(0,d^2),
                         nrow=d)
    #cat("ESTIMATESOFE", "\n")
    #print(ESTIMATESOFE)
    for (i in 1:d)
    {ESTIMATESOFE[i,i]<-rin10prodequalvector(A[,i])
    }
    
    for(i in 1:(d-1))
    {u<-A[,i]
    
    for(j in (i+1):d)
    {
      v<-A[,j]
      ESTIMATESOFE[i,j]<-rin10prod(u,v)
      ESTIMATESOFE[j,i]<-ESTIMATESOFE[i,j]
      #    cat("ESTIMATESOFE[i,j]=",
      #        ESTIMATESOFE[i,j], "i=", i, "j=", j, "\n")   
    }
    }
    
    
    cat("A", "\n")
    print(A)
    cat("ESTIMATESOFE", "\n")
    print(ESTIMATESOFE)
    means<-rep(0,d)
    for (i in 1:d)
    {means[i]<-mean(A[,i])
    }
    cat("Means 1st line=", means)
    ESTIMATESOFE1<-rbind(means,ESTIMATESOFE)
    cat("ESTIMATESOFE1", "\n")
    print(ESTIMATESOFE1)
    ESTIMATESOFE2<-cbind(c(1,means),ESTIMATESOFE1)
    cat("ESTIMATESOFE2", "\n")
    print(ESTIMATESOFE2)
    
    
    cat("INVERSE MATRIX THAT GIVES INFLUENCE", "\n")
    INVERSE<-solve(ESTIMATESOFE2)
    cat("INVERSE MATRIX", "\n")
    print(INVERSE)
    #INF<-rep(0,d+1) 
    #INF<-INVERSE%*%r*c(1,x[])
  }
  INVERSE
}


# FUNCTIONS FOR TABLES 3 AND 4 - MICRO ARRAY DATA
# Note to the reader:   For Tables 3 (IN A) ), 4 (in B)),  the function # 1  is used and creates the results. The other functions are in the background, i.e. inside function 1, then inside 2 etc. 
# The results are the TOTAL RINFIN for each case and then the ordered values from the smallest to the largest. It can be properly adjusted for other data sets.  Note that the last column in the micro-array data set used is the y-column.
# 
# You may use the code with data sets of size (1+multiples of 100). D=regression data with last column the y-vector of responses. Blocks of 100 covariates are used because n>100. The block size ???100??? can be changed.
# 
# DURATION TO OBTAIN THE RESULTS FOR THIS DATA SET IN R-STUDIO: 3 days, if I recall correctly for the results of Table 3. Similar duration for Table 4.
# 
# A) FUNCTIONS FOR TABLE 3-THE RINFIN RESULTS

# 1)
rinsq11EX3.1  (D) 
{    
  
  # This gives 
  lastcol<-ncol(D)
  cat("LAST COL =", lastcol, "\n")
  print(lastcol)
  n<-nrow(D)
  print(n)
  cat("NUMBER OF ROWS =", n, "\n")
  rinftotscore<-rep(0,n)
  y<-D[,lastcol]
  cat("LAST COLUMN=",lastcol,"\n")
  cat("NUMBER OF ROWS=",n, "\n")
  ###THE LINES BEFORE COULD BE ADJUSTED FOR BLOCK OF COLUMNS
  ####DIFFERENT THAN 100 TO CALCULATE SEVERAL RINFIN VALUES
  REPET<-(lastcol-1)/100
  SCORES<-rep(0,n)
  for (j in 1:REPET)
  {
    L<-1+100*(j-1)
    U<-100*j
    # FUNCTION 2
    x<-rinsq11EX3.2(y, D[,L:U])
    SCORES<-SCORES+x
    
  }
  cat("WHERE ARE THE SCORES?", "\n")
  print(cbind(1:n,SCORES))
  
  colorder<-c(1:n)
  orinf<-order(SCORES)
  cat("THESE ARE THE ORDERED SCORES FOR MICROARRAY-DATA", "\n")
  print(cbind(colorder[orinf], SCORES[orinf]))
  
}

#2) 
rinsq11EX3.2 (y,x) 
{
  
  
  {
    
    
    #STARTING { NOT TO BE CONFUSED BY IT
    n<-nrow(x)
    d<-ncol(x)
    cat("n=", n, "\n")
    RINFIN<-rep(0,n)
    for (i in 1:n)
    { #STARTS ONE LOOP FOR CALCULATING ALL THE 
      #TAKING THE i-th CASE OUT, 
      #DOING LEAST SQUARES, FIND RESIDUAL
      # OF THE i-th CASE WITH RESPECT TO 
      #THE REGRESSION COEFFICIENT b OBTAINED
      #BELOW
      cat("ITERATION I=", i, "\n")
      AVoitheiax<-x[-i,]
      AVoitheiay<-y[-i]
      results<-lm(AVoitheiay ~AVoitheiax) 
      b<- results$coefficients
      cat("b=", b, "\n")
      residual<-y[i]-sum(b*c(1,x[i,]))
      cat("Residual=", residual, "\n")
      
      # CALCULATING THE ESTIMATE OF E MATRIX
      #AND ITS INVERSE TO FIND INFLUENCE AND RINFIN
      # FUNCTION 3  
      INVERSE<-rin10prodmatrix(AVoitheiax)  
      INFLUENCEFUNCTIONS<-rep(0,d+1)
      cat("INFLUENCE FUNCTIONS=", INFLUENCEFUNCTIONS)
      residualtimesvector<-residual*c(1,x[i,])
      INFLUENCEFUNCTIONS<-INVERSE%*%residualtimesvector
      cat("INFLUENCE FUNCTIONS FOR CASE i=",i, "\n" )
      print(INFLUENCEFUNCTIONS)
      
      #CALCULATING INFLUENCE IN A WAY THAT HELPS CALCULATING
      #DERIVATIVES AND RINFIN
      
      xvoith<-c(1,x[i,])
      INFLUENCEFUNCTIONSMANUAL<-rep(0,d+1)
      COMPONENT<-rep(0,d+1)
      for (j in 1:(d+1))
      {COMPONENT[j]<-sum (INVERSE[j,]*xvoith)
      INFLUENCEFUNCTIONSMANUAL[j]<-
        residual*COMPONENT[j]
      # residual*sum (INVERSE[j,]*xvoith)
      }
      cat("INFLUENCE FUNCTIONS MANUAL FOR CASE i=",i, "\n" )
      print(INFLUENCEFUNCTIONSMANUAL)
      #MATRIX OF DERIVATIVES HAS (d+1) COLUMNS AS THE NUMBER OF 
      #INFLUENCE FUNCTIONS AND THERE ARE d PARTIAL DERIVATIVES
      DERIVATIVESOFINFLUENCEFUNCTIONS<-matrix(rep(0,d*(d+1)),ncol=d+1) 
      for (m in 1:d)
      {
        for (j in 1:(d+1))
        {DERIVATIVESOFINFLUENCEFUNCTIONS[m,j]<--b[m+1]*COMPONENT[j]
        +residual*INVERSE[j,m+1]
        
        }
        
        #cat("DERIVATIVES OF INFLUENCE FUNCTIONS")
        #print(DERIVATIVESOFINFLUENCEFUNCTIONS)
        
        
      }
      cat("DERIVATIVES OF INFLUENCE FUNCTIONS")
      print(DERIVATIVESOFINFLUENCEFUNCTIONS)
      INFLUENCEITHCOORDINATEOFX<-rep(0,d)
      for(k in 1:d)
      {
        INFLUENCEITHCOORDINATEOFX[k]<-INFLUENCEFUNCTIONS[k+1]+
          sum(DERIVATIVESOFINFLUENCEFUNCTIONS[k,]*xvoith)
        # CLOSES THE "FOR" 
        cat("INFLUENCEITHCOORDINATEOFX[k]") 
        print(INFLUENCEITHCOORDINATEOFX[k])
      }    
      
      RINFIN[i]    <-sum((INFLUENCEITHCOORDINATEOFX)^2/n)
      cat("RINFINI")
      print(RINFIN[i])
      #cat("INFLUENCEITHCOORDINATEOFX[k]")
      #print(RINFIN[i])
      
    }
    
    EXTRACTRINFIN<-RINFIN
    
    #cat("RINFIN VALUES FOR EACH CASE")
    #colorder<-c(1:n)
    #orinf<-order(RINFIN)
    #cat("THESE ARE THE RINFIN ORDERED SCORES", "\n")
    #print(cbind(colorder[orinf], RINFIN[orinf]))
  }
  
}

#3)
rin10prodmatrix(A) 
{
  {
    
    
    print(A)
    n<-nrow(A)
    d<-ncol(A)
    #A<-matrix(A,ncol=d)
    #print(A)
    ESTIMATESOFE<-matrix(rep(0,d^2),
                         nrow=d)
    #cat("ESTIMATESOFE", "\n")
    #print(ESTIMATESOFE)
    for (i in 1:d)
      # FUNCTION 4 
    {ESTIMATESOFE[i,i]<-rin10prodequalvector(A[,i])
    }
    
    for(i in 1:(d-1))
    {u<-A[,i]
    
    for(j in (i+1):d)
    {
      v<-A[,j]
      # FUNCTION 5 
      ESTIMATESOFE[i,j]<-rin10prod(u,v)
      ESTIMATESOFE[j,i]<-ESTIMATESOFE[i,j]
      #    cat("ESTIMATESOFE[i,j]=",
      #        ESTIMATESOFE[i,j], "i=", i, "j=", j, "\n")   
    }
    }
    
    
    cat("A", "\n")
    print(A)
    cat("ESTIMATESOFE", "\n")
    print(ESTIMATESOFE)
    means<-rep(0,d)
    for (i in 1:d)
    {means[i]<-mean(A[,i])
    }
    cat("Means 1st line=", means)
    ESTIMATESOFE1<-rbind(means,ESTIMATESOFE)
    cat("ESTIMATESOFE1", "\n")
    print(ESTIMATESOFE1)
    ESTIMATESOFE2<-cbind(c(1,means),ESTIMATESOFE1)
    cat("ESTIMATESOFE2", "\n")
    print(ESTIMATESOFE2)
    
    
    cat("INVERSE MATRIX THAT GIVES INFLUENCE", "\n")
    INVERSE<-solve(ESTIMATESOFE2)
    cat("INVERSE MATRIX", "\n")
    print(INVERSE)
    #INF<-rep(0,d+1) 
    #INF<-INVERSE%*%r*c(1,x[])
  }
  INVERSE
  
  
}

#4)
rin10prodequalvector (u) 
{
  { n<-length(u)
  #pr<-matrix(rep(0, n^2), nrow=n)
  # cat("length of vector=",n,"\n")
  # cat("PRODUCT MATRIX", "\n")
  #  print(pr)
  MEQVECTOR<-mean(u^2)
  cat("ESTIMATE OF SECOND MOMENT=",MEQVECTOR, "\n" )
  # cat("uv", "\n")
  #print(cbind(u,v))
  }
  MEQVECTOR
  
}
5) 
rin10prod  (u,v) 
{
  { n<-length(u)
  #pr<-matrix(rep(0, n^2), nrow=n)
  # cat("length of vector=",n,"\n")
  # cat("PRODUCT MATRIX", "\n")
  #  print(pr)
  pr<-rep(0,n)
  for (i in 1:n)
  {pr[i]<- u[i]*sum(v)}
  # cat("pr","\n")
  # print(pr)
  S<-sum(pr)
  cat("SUM OF PRODUCTS=", S, "\n")
  M<-sum(pr)/n^2
  cat("ESTIMATE OF MEAN PRODUCT=", M,"\n")
  # cat("uv", "\n")
  #print(cbind(u,v))
  }
  M
  
}


#B) FUNCTIONS FOR TABLE 4- THE RINFINABS RESULTS

#1)
rin10EX3.1 (D) 
{    
  
  # This gives 
  lastcol<-ncol(D)
  cat("LAST COL =", lastcol, "\n")
  print(lastcol)
  n<-nrow(D)
  print(n)
  cat("NUMBER OF ROWS =", n, "\n")
  rinftotscore<-rep(0,n)
  y<-D[,lastcol]
  cat("LAST COLUMN=",lastcol,"\n")
  cat("NUMBER OF ROWS=",n, "\n")
  ###THE LINES BEFORE COULD BE ADJUSTED FOR BLOCK OF COLUMNS
  ####DIFFERENT THAN 100 TO CALCULATE SEVERAL RINFIN VALUES
  REPET<-(lastcol-1)/100
  SCORES<-rep(0,n)
  for (j in 1:REPET)
  {
    L<-1+100*(j-1)
    U<-100*j
    # FUNCTION 2 
    x<-rin10EX3.2(y, D[,L:U])
    SCORES<-SCORES+x
    
  }
  cat("WHERE ARE THE SCORES?", "\n")
  print(cbind(1:n,SCORES))
  
  colorder<-c(1:n)
  orinf<-order(SCORES)
  cat("THESE ARE THE ORDERED SCORES FOR MICROARRAY-DATA", "\n")
  print(cbind(colorder[orinf], SCORES[orinf]))
  
}



#2)
rin10EX3.2 (y,x) 
{  
  {       
    #STARTING { NOT TO BE CONFUSED BY IT
    n<-nrow(x)
    d<-ncol(x)
    cat("n=", n, "\n")
    RINFIN<-rep(0,n)
    for (i in 1:n)
    { #STARTS ONE LOOP FOR CALCULATING ALL THE 
      #TAKING THE i-th CASE OUT, 
      #DOING LEAST SQUARES, FIND RESIDUAL
      # OF THE i-th CASE WITH RESPECT TO 
      #THE REGRESSION COEFFICIENT b OBTAINED
      #BELOW
      cat("ITERATION I=", i, "\n")
      AVoitheiax<-x[-i,]
      AVoitheiay<-y[-i]
      results<-lm(AVoitheiay ~AVoitheiax) 
      b<- results$coefficients
      cat("b=", b, "\n")
      residual<-y[i]-sum(b*c(1,x[i,]))
      cat("Residual=", residual, "\n")
      
      # CALCULATING THE ESTIMATE OF E MATRIX
      #AND ITS INVERSE TO FIND INFLUENCE AND RINFIN
      # FUNCTION 3   
      INVERSE<-rin10prodmatrix(AVoitheiax)  
      INFLUENCEFUNCTIONS<-rep(0,d+1)
      cat("INFLUENCE FUNCTIONS=", INFLUENCEFUNCTIONS)
      residualtimesvector<-residual*c(1,x[i,])
      INFLUENCEFUNCTIONS<-INVERSE%*%residualtimesvector
      cat("INFLUENCE FUNCTIONS FOR CASE i=",i, "\n" )
      print(INFLUENCEFUNCTIONS)
      
      #CALCULATING INFLUENCE IN A WAY THAT HELPS CALCULATING
      #DERIVATIVES AND RINFIN
      
      xvoith<-c(1,x[i,])
      INFLUENCEFUNCTIONSMANUAL<-rep(0,d+1)
      COMPONENT<-rep(0,d+1)
      for (j in 1:(d+1))
      {COMPONENT[j]<-sum (INVERSE[j,]*xvoith)
      INFLUENCEFUNCTIONSMANUAL[j]<-
        residual*COMPONENT[j]
      # residual*sum (INVERSE[j,]*xvoith)
      }
      cat("INFLUENCE FUNCTIONS MANUAL FOR CASE i=",i, "\n" )
      print(INFLUENCEFUNCTIONSMANUAL)
      #MATRIX OF DERIVATIVES HAS (d+1) COLUMNS AS THE NUMBER OF 
      #INFLUENCE FUNCTIONS AND THERE ARE d PARTIAL DERIVATIVES
      DERIVATIVESOFINFLUENCEFUNCTIONS<-matrix(rep(0,d*(d+1)),ncol=d+1) 
      for (m in 1:d)
      {
        for (j in 1:(d+1))
        {DERIVATIVESOFINFLUENCEFUNCTIONS[m,j]<--b[m+1]*COMPONENT[j]
        +residual*INVERSE[j,m+1]
        
        }
        
        #cat("DERIVATIVES OF INFLUENCE FUNCTIONS")
        #print(DERIVATIVESOFINFLUENCEFUNCTIONS)
        
        
      }
      cat("DERIVATIVES OF INFLUENCE FUNCTIONS")
      print(DERIVATIVESOFINFLUENCEFUNCTIONS)
      INFLUENCEITHCOORDINATEOFX<-rep(0,d)
      for(k in 1:d)
      {
        INFLUENCEITHCOORDINATEOFX[k]<-INFLUENCEFUNCTIONS[k+1]+
          sum(DERIVATIVESOFINFLUENCEFUNCTIONS[k,]*xvoith)
        # CLOSES THE "FOR" 
        cat("INFLUENCEITHCOORDINATEOFX[k]") 
        print(INFLUENCEITHCOORDINATEOFX[k])
      }    
      
      RINFIN[i]    <-sum(abs(INFLUENCEITHCOORDINATEOFX)/n)
      cat("RINFINI")
      print(RINFIN[i])
      #cat("INFLUENCEITHCOORDINATEOFX[k]")
      #print(RINFIN[i])
      
    }  
    EXTRACTRINFIN<-RINFIN 
    #cat("RINFIN VALUES FOR EACH CASE")
    #colorder<-c(1:n)
    #orinf<-order(RINFIN)
    #cat("THESE ARE THE RINFIN ORDERED SCORES", "\n")
    #print(cbind(colorder[orinf], RINFIN[orinf]))
  }
}

#3) 
#rin10prodmatrix  IS ALREADY IN FOR TABLE 3 SO THE OTHER FUNCTIONS ARE IN.




```

automatically created on 2018-11-02