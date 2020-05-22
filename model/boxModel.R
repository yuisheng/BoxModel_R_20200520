### =======================================================================
### = boxModel
### = Zhen Zhang
### = 10/13/2017
### =----------------------------------------------------------------------
### = NOTES
### =  ( 1): 2-box model for methane, delta13C, and deltaD, (C2H6 for future).
### =  ( 2): Adapted from A. Turner's box model.
### =----------------------------------------------------------------------
### = INPUTS
### =  ( 1): St       -- (Vector) contains time stamp.
### =  ( 2): S        -- (Matrix) Emission sources for the box model (6 columns contains 12CH4, 13CH4, and CH3D for NH and SH).
### =  ( 3): OH       -- (Vector) time series of OH concentration (molec/cm3)
### =  ( 4): params   -- (List)   Structure with parameters for the box model.
### =  ( 5): conc     -- (Vector) time series of CH4 concentration
### =  ( 6): inits    -- (Vector) initial conditions of atmospheric concentraion, isotopic ratio (status from spinup)
### =  ( 7): runMode  -- String: inverse mode for calculating OH (E and C are known) or Emissions (OH and C are known)
### =  ( 8): states   -- (Vector) state variables from restart file
### =----------------------------------------------------------------------
### = OUTPUTS
### =  ( 1): dy -- (Matrix) Changes in the box model concentrations (forward mode).
### =  ( 2): OH -- (Vector) OH concentration (OH_inverse mode) 
### =======================================================================
source('model/tools.R')

# The model is designed to be running in forward or inverse mode.
# - For forward mode, the run starts in the specified startyear to endyear with
#   inputs of emissions (E) and time series of OH. Output is simulated atmospheric
#   CH4 concentration (C) (including isotopes)
# - For OH_inverse mode, the run starts from startyear to end with inputs of emissions (E)
#   and observed atmospheric CH4 concentration (C). Output is retrived time series of OH using a bisection method
# - For E-Inverse mode, the run starts from startyear to endyear with inputs of time series of OH
#   and observed atmospheric CH4 concentration (C). Output is reconstucted Emissions (E)

boxModel <-function(St,S,OH,params,obs,inits,iKIE,runMode,...){
  
  
  if(runMode=='forward'){
    #Set output matrix
    # - Column 1: 12CH4 in the Northern Hemisphere (out[,1])  Unit: ppb
    # - Column 2: 12CH4 in the Southern Hemisphere (out[,2])  Unit: ppb
    # - Column 3: 13CH4 in the Northern Hemisphere (out[,3])  Unit: ppb
    # - Column 4: 13CH4 in the Southern Hemisphere (out[,4])  Unit: ppb
    # - Column 5: CH3D  in the Northern Hemisphere (out[,5])  Unit: ppt
    # - Column 6: CH3D  in the Southern Hemisphere (out[,6])  Unit: ppb
    
    #forward mode
    out <- matrix(0,nrow=length(St),ncol=6)
    
    OHfac   <- OH * params$DaysToS * params$YrToDay
    k_12ch4 <- OHfac * params$k_12ch4 * 0.5
    k_13ch4 <- OHfac * params$k_12ch4 * iKIE * 0.5
    k_ch3d  <- OHfac * params$k_ch3d
    
    # Calculate C13mass from total CH4 emissions and deltaCH4 
    Sc      <- matrix(0,nrow=length(St),ncol=6) # Emissions of 12CH4, 13CH4, and CH3D in ppb
    Sc[,3] <- calC13mass(S[,1]/params$mm_ch4,S[,3])
    Sc[,4] <- calC13mass(S[,2]/params$mm_ch4,S[,4])
    Sc[,1] <- S[,1]/params$mm_ch4 - Sc[,3]
    Sc[,2] <- S[,2]/params$mm_ch4 - Sc[,4]
    
    
    #Convert emission (Unit: Tg CH4/yr) to ppb by using factor from params$mm_ch4
    for(i in 1:length(St)){
      if(i == 1){
        # CH412_NH
        out[i,1] <- inits[1] + Sc[i,1] - k_12ch4[i]*inits[1] + (inits[2] - inits[1])/2
        #CH412_SH
        out[i,2] <- inits[2] + Sc[i,2] - k_12ch4[i]*inits[2] + (inits[1] - inits[2])/2
        # CH413_NH
        out[i,3] <- inits[3] + Sc[i,3] - k_13ch4[i]*inits[3] + (inits[4] - inits[3])/2
        #CH413_SH
        out[i,4] <- inits[4] + Sc[i,4] - k_13ch4[i]*inits[4] + (inits[3] - inits[4])/2
        #CH3D_NH
        out[i,5] <- inits[5] + Sc[i,5] - k_ch3d[i]*inits[5] + (inits[6] - inits[5])/2
        #CH3D_SH
        out[i,6] <- inits[6] + Sc[i,6] - k_ch3d[i]*inits[6] + (inits[5] - inits[6])/2
      }else{
        # CH412_NH
        out[i,1] <- out[i-1,1] + Sc[i,1] - k_12ch4[i]*out[i-1,1] + (out[i-1,2] - out[i-1,1])/2
        #CH412_SH
        out[i,2] <- out[i-1,2] + Sc[i,2] - k_12ch4[i]*out[i-1,2] + (out[i-1,1] - out[i-1,2])/2
        # CH413_NH
        out[i,3] <- out[i-1,3] + Sc[i,3] - k_13ch4[i]*out[i-1,3] + (out[i-1,4] - out[i-1,3])/2
        #CH413_SH
        out[i,4] <- out[i-1,4] + Sc[i,4] - k_13ch4[i]*out[i-1,4] + (out[i-1,3] - out[i-1,4])/2
        #CH3D_NH
        out[i,5] <- out[i-1,5] + Sc[i,5] - k_ch3d[i]*out[i-1,5] + (out[i-1,6] - out[i-1,5])/2
        #CH3D_SH
        out[i,6] <- out[i-1,6] + Sc[i,6] - k_ch3d[i]*out[i-1,6] + (out[i-1,5] - out[i-1,6])/2
      }
      
    }
    
    return(out)
    
  }else if(runMode=='OH_inverse'){
    #OH inverse mode
    #global average CH4 concencetration (obs) and Emissions (S) are known 
    # we use a bisection method to get an estimate for OH at each time step
    # we don't distingush the hemispheric CH4 state but assume a one box condition
    # Out(t) = Out(t-1) + S(t) - k * Out(t-1)
    #we get k (t) = 1 + (S(t) - Out(t))/Out(t-1)

    #Set input matrix (same nrows and ncols as output mateix in forward mode)
    # - Column 5:   CH3D in the Northern Hemisphere (out[,5])  Unit: ppt
    # - Column 6:   CH3D in the Southern Hemisphere (out[,6])  Unit: ppb
    out_obs <- matrix(0,nrow=length(St),ncol=6) # observed 12CH4, 13CH4, and CH3D at each time step
    out_sim <- matrix(0,nrow=length(St),ncol=6) # simulated 12CH4, 13CH4, and CH3D at each time step
    
    #out_obs value are known from observations (obs)
    out_obs[,1] <- calC12mass(obs[,1],obs[,3])  #Column 1: 12CH4 in the Northern Hemisphere (out[,1])  Unit: ppb
    out_obs[,2] <- calC12mass(obs[,2],obs[,3])  #Column 2: 12CH4 in the Southern Hemisphere (out[,2])  Unit: ppb
    out_obs[,3] <- obs[,1] - out_obs[,1]            #Column 3: 13CH4 in the Northern Hemisphere (out[,3])  Unit: ppb
    out_obs[,4] <- obs[,2] - out_obs[,2]            #Column 4: 13CH4 in the Southern Hemisphere (out[,4])  Unit: ppb
    
    
    # Calculate C13mass from total CH4 emissions and deltaCH4 
    Sc      <- matrix(0,nrow=length(St),ncol=6) # Emissions of 12CH4, 13CH4, and CH3D in ppb
    Sc[,3] <- calC13mass(S[,1]/params$mm_ch4,S[,3])
    Sc[,4] <- calC13mass(S[,2]/params$mm_ch4,S[,4])
    Sc[,1] <- S[,1]/params$mm_ch4 - Sc[,3]
    Sc[,2] <- S[,2]/params$mm_ch4 - Sc[,4]
    
    
    #Set output vector
    k <- rep(1,length(vec.St))
    
    for(i in 1:length(vec.St)){
      if(i ==1){
        #create a function with k as the only var
        calk <- function(x){ # x is k[i]
        #   return(p12 + p13 + e - k * p12 - k * iKIE * p13 - o)
           return(  (inits[1] + Sc[i,1] - x*inits[1] + (inits[2] - inits[1])/2 - out_obs[i,1]) +
                    (inits[2] + Sc[i,2] - x*inits[2] + (inits[1] - inits[2])/2 - out_obs[i,2]) +
                    (inits[3] + Sc[i,3] - x*iKIE*inits[3] + (inits[4] - inits[3])/2 - out_obs[i,3]) +
                    (inits[4] + Sc[i,4] - x*iKIE*inits[4] + (inits[3] - inits[4])/2 - out_obs[i,4])
                  )
         }
        k[i] <- bisection(calk,0.1,0.02,n=1e5)
        
        #calculate in forward mode for nect step
        out_sim[i,1] <- inits[1] + Sc[i,1] - k[i]*inits[1] + (inits[2] - inits[1])/2
        out_sim[i,2] <- inits[2] + Sc[i,2] - k[i]*inits[2] + (inits[1] - inits[2])/2
        out_sim[i,3] <- inits[3] + Sc[i,3] - k[i]*iKIE*inits[3] + (inits[4] - inits[3])/2
        out_sim[i,4] <- inits[4] + Sc[i,4] - k[i]*iKIE*inits[4] + (inits[3] - inits[4])/2
      }else{
        #create a function with k as the only var
        calk <- function(x){ # x is k[i]
          #   return(p12 + p13 + e - k * p12 - k * iKIE * p13 - o)
          return(  (out_sim[i-1,1] + Sc[i,1] - x*out_sim[i-1,1] + (out_sim[i-1,2] - out_sim[i-1,1])/2 - out_obs[i,1]) +
                     (out_sim[i-1,2] + Sc[i,2] - x*out_sim[i-1,2] + (out_sim[i-1,1] - out_sim[i-1,2])/2 - out_obs[i,2]) +
                     (out_sim[i-1,3] + Sc[i,3] - x*iKIE*out_sim[i-1,3] + (out_sim[i-1,4] - out_sim[i-1,3])/2 - out_obs[i,3]) +
                     (out_sim[i-1,4] + Sc[i,4] - x*iKIE*out_sim[i-1,4] + (out_sim[i-1,3] - out_sim[i-1,4])/2 - out_obs[i,4])
          )
        }
        k[i] <- bisection(calk,0.1,0.02,n=1e5)
        
        #calculate in forward mode for next step
        out_sim[i,1] <- out_sim[i-1,1] + Sc[i,1] - k[i]*out_sim[i-1,1] + (out_sim[i-1,2] - out_sim[i-1,1])/2
        out_sim[i,2] <- out_sim[i-1,2] + Sc[i,2] - k[i]*out_sim[i-1,2] + (out_sim[i-1,1] - out_sim[i-1,2])/2
        out_sim[i,3] <- out_sim[i-1,3] + Sc[i,3] - k[i]*iKIE*out_sim[i-1,3] + (out_sim[i-1,4] - out_sim[i-1,3])/2
        out_sim[i,4] <- out_sim[i-1,4] + Sc[i,4] - k[i]*iKIE*out_sim[i-1,4] + (out_sim[i-1,3] - out_sim[i-1,4])/2
      }
    }
    return(k/params$k_12ch4*2/params$DaysToS/params$YrToDay)
    
  }

}
