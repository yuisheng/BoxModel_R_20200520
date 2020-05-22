### =======================================================================
### = getParameters
### = Zhen Zhang
### = 10/30/2017
### =----------------------------------------------------------------------
### = NOTES
### =  ( 1): Computes the parameters that are needed for the box model.
### =----------------------------------------------------------------------
### = INPUTS
### =  ( 1): St -- Our time vector.
### =----------------------------------------------------------------------
### = OUTPUTS
### =  ( 1): params -- Structure with the parameters for the box model.
### =======================================================================
source('model/tools.R')
getParameters <- function(St){
  DaysToS <- 60 * 60 * 24;         # Days to Seconds
  YrToDay <- 365.25;               # Years to Days
  tau_NS  <- 1.0;                  # Interhemispheric exchange rate (years)
  
  ### Box model parameters and unit conversions
  atm.airmass        <- 5.1352e21;                #Trenberth 2004; Rough: Dry air mass of atmosphere (g)
  #atm.airmass        <- 5.1480e21;               #Trenberth 2004; Mean air mass of atmosphere from  (equal to 5.22371*1E15*985.5 mbar)
  #atm.airmass        <- (5.22371*1e15*960)*1000; #Trenberth 2004; Mass of atmosphere from Fung 1991 using Trenberth conversion of 960 mbar from Fung
  #atm.airmass        <- 5.15e21;                 #Turner 2017; Rough: Total mass of atmosphere in g
  #m_air    <- 28.8;                              #Turner 2017; Average molar mass of the atmosphere in g/mol (mostly N2 and O2)
  m_air    <- 28.966;                             #Standard accepted
  m_ch4   <- 16.04;                               # Molar mass of CH4
  m_mcf   <- 133.4;                               # Molar mass of methyl-chloroform
  mConv   <- atm.airmass/m_air/1e12*1e-9;	        # Factor to convert mixing ratios to Tg
   
  #mm_ch4  <- 2.861111;           # Turner 2017; Convert CH4 mixing ratios to Tg
  mm_ch4  <- 2.7476               # Prather 2012; Convert CH4 mixing ratios to Tg
  #mm_ch4  <- 2.845               # Prbir's estimate that considers vertical profile
  mm_mcf  <- m_mcf * mConv;       # Convert MCF mixing ratios to Tg
  
  KIE.OH     <- 1.0054             # -5.4% from Sapart 2012; -3.9 Kinetic isotope effect for OH sink (k12/k13) from Burkholder et al. (2015) (not used)
  KIE.Soil   <- 1.020              # -20 Kinetic isotope effect for Soil removal (k12/k13) from Lassey 2007 ;-18 from Sapart et al., (2012) (not used)
  KIE.Stra   <- 1.013              # -12.8 Kinetic isotope effect for Statospheric loss (k12/k13) from Sapart et al., (2012) 
  KIE.Cl     <- 1.066              # -61 Kinetic isotope effect for Cl loss (k12/k13) from Crowley 1999; 1.064 from Lassey 2007
  iKIE.OH    <- round(1/KIE.OH,3);    # Inverse of KIE for OH (for later)
  iKIE.Soil  <- round(1/KIE.Soil,3)   # Inverse of KIE for Soil removal (for later)
  iKIE.Stra  <- round(1/KIE.Stra,3)   # Inverse of KIE for Stratospheric loss (for later)
  iKIE.Cl    <- round(1/KIE.Cl,3)     # Inverse of KIE for stratospheric loss (for later)
  iKIE.test <- 1/1.008                # 
  #Calculate  the sink-weighted fractionation factor of all sink processes
  #k_12ch4 <- 3.395e-15;            # reaction rate of OH with 12CH4 (cm3/molec/s) Turner et al., 2016
  
  k_12ch4 <- 3.105e-15;            # reaction rate of OH with 12CH4 (cm3/molec/s) from Sander SP et al., (2011)
  k_13ch4 <- k_12ch4 * iKIE.test;  # reaction rate of OH with 13CH4 (cm3/molec/s)
  k_mcf_A <- 5.66e-15;             # reaction rate with MCF (computed such that lifetime is about 5.5 years, Talukdar et al 1992, taken at 273K)
  k_mcf_B <- 6.05e-15;             # reaction rate determined by AJT
  k_ch3d  <- 5.66e-15;             #reaction rate with CH3D (need update)
  
  frac_ch4nh <- 0.67                # fraction of CH4 emissions from northern hemisphere (not used)
  frac_mcfnh <- 0.7                # fraction of MCF emissions from northern hemisphere
  total_ch4  <- 550                # assumed average total CH4 emissions. Unit: Tg CH4/yr
  
  ### Initial conditions for the box model
  nh_ch4 <- 1662                        # ppb
  sh_ch4 <- 1512                        # ppb
  nh_deltac <- -47.5  # per Mille
  sh_deltac <- -47.5 # per Mille
  nh_mcf   <- 85                          # ppt
  sh_mcf   <- 75                          # ppt
  
  # Assemble the ICs into a vector
  IC <- c(generateIC(nh_ch4,sh_ch4,nh_deltac,sh_deltac),nh_mcf,sh_mcf)
  
  ### Make a structure with the parameters
  # Unit conversions
  #params <- list(mm_ch4,mm_mcf,YrToDay,k_12ch4,k_13ch4,k_mcf,tau_NS,Tspan,IC,opts)
  params <- list(mm_ch4=mm_ch4,mm_mcf=mm_mcf,DaysToS=DaysToS,YrToDay=YrToDay,k_12ch4=k_12ch4,k_13ch4=k_13ch4,iKIE=iKIE.test,k_mcf_A=k_mcf_A,k_mcf_B=k_mcf_B,
                 tau_NS=tau_NS,IC=IC,frac_ch4nh=frac_ch4nh,frac_mcfnh=frac_mcfnh,
                 total_ch4=total_ch4)
  return(params)
}


