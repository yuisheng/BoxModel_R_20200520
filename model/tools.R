### =======================================================================
### = calC13mass
### = Zhen Zhang
### = 10/30/2017
### =----------------------------------------------------------------------
### = NOTES
### =  ( 1): functions for calculating deltaC13 mass and C12 mass 
### =----------------------------------------------------------------------
### = INPUTS
### =  ( 1): ems     -- CH4 Emission source in ppb
### =  ( 2): deltaC  -- measured 13C/12C ratio 
### =----------------------------------------------------------------------
### = OUTPUTS
### =  ( 1): deltaCmass -- 13C Mass (ppb)
### =======================================================================
#Constant Parameter
VPDB    <- 0.0112372               # Vienna-PDB from Sapart et al., (2012)

calC13mass <- function(ems,deltaC){
  U <- (deltaC/1000 + 1)*VPDB
  return(U/(1+U)*ems)
}
calC12mass <- function(ems,deltaC){
  U <- (deltaC/1000 + 1)*VPDB
  return(1/(1+U)*ems)  
}
### =======================================================================
### = calC13ratio
### = Zhen Zhang
### = 10/30/2017
### =----------------------------------------------------------------------
### = NOTES
### =  ( 1): A function for calculating 13C/12C ratio
### =----------------------------------------------------------------------
### = INPUTS
### =  ( 1): C12mass  -- 12CH4 Mass in ppb
### =  ( 2): C13mass  -- 13CH4 Mass in ppb
### =----------------------------------------------------------------------
### = OUTPUTS
### =  ( 1): R -- 13C:12C Ratio (per mille)
### =======================================================================
calC13ratio <- function(C12mass,C13mass){
  return(((C13mass/C12mass/VPDB - 1) * 1e3))
}

### =======================================================================
### = generateIC
### = Zhen Zhang
### = 10/30/2017
### =----------------------------------------------------------------------
### = NOTES
### =  ( 1): A function for generating a input vector of Initial conditions
### =----------------------------------------------------------------------
### = INPUTS
### =  ( 1): CH4_conc_NH  -- CH4 concentration in NH (ppb)
### =  ( 2): CH4_conc_SH  -- CH4 concentration in SH (ppb)
### =  ( 3): CH4_sig_NH  -- CH4 signature in NH (per mille)
### =  ( 3): CH4_sig_SH  -- CH4 signature in SH (per mille)
### =----------------------------------------------------------------------
### = OUTPUTS
### =  ( 1): R -- 13C:12C Ratio (per mille)
### =======================================================================
generateIC <- function(CH4_conc_NH,CH4_conc_SH,CH4_sig_NH,CH4_sig_SH){
  nh_12ch4 <- calC12mass(CH4_conc_NH,CH4_sig_NH)     #unit: ppb
  sh_12ch4 <- calC12mass(CH4_conc_SH,CH4_sig_SH)     #unit: ppb
  nh_13ch4 <- calC13mass(CH4_conc_NH,CH4_sig_NH)     #unit: ppb
  sh_13ch4 <- calC13mass(CH4_conc_SH,CH4_sig_SH)     #unit: ppb
  return(c(nh_12ch4,sh_12ch4,nh_13ch4,sh_13ch4))
}

KIEtoEpsilon <- function(KIE){
  return((1/KIE -1)*1000)
}

EpsilontoKIE <- function(epsilon){
  return(1/(epsilon/1000 + 1))
}

### =======================================================================
### = bisection
### = Zhen Zhang
### = 10/12/2018
### =----------------------------------------------------------------------
### = NOTES
### =  ( 1): A bisection tool for find the best solution for parameters
### =----------------------------------------------------------------------
### = INPUTS
### =  ( 1): f  -- R function
### =  ( 2): a  -- starting range
### =  ( 3): b  -- ending range
### =  ( 3): n  -- number of iteration
### =----------------------------------------------------------------------
### = OUTPUTS
### =  ( 1): y -- root
### =======================================================================
bisection <- function(f, a, b, n, tol=1e-8){
  # # If the signs of the function at the evaluated points, a and b, stop the function and return message.
  # if (!(f(a) < 0) && (f(b) > 0)) {
  #   stop('signs of f(a) and f(b) differ')
  # }else if ((f(a) > 0) && (f(b) < 0)){
  #   stop('signs of f(a) and f(b) differ')
  # }
  
  for(i in 1:n){
    c <- (a + b) / 2  #calculate midpoint
    
    # If the function equals 0 at the midpoint or the midpoint is below the desired tolerance, stop the 
    # function and return the root.
    if( f(c) == 0 || (abs(b - a)/2 < tol) ){
      return(c)
    }
    
    # If another iteration is required,
    # check the signs of the function at the points c and a and reassign
    # a or b accordingly as the midpoint to be used in the next iteration.
    ifelse(sign(f(c)) == sign(f(a)),
           a <- c,
           b <- c)
  }
  # If the max number of iterations is reached and no root has been found,
  # return message and end function.
  print('Too many iterations')
}
### =======================================================================
### = END
### =======================================================================
