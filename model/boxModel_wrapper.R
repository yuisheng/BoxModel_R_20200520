### =======================================================================
### = boxModel_wrapper.R
### = rewriten by Zhen Zhang
### = 10/15/2017
### =======================================================================
### =----------------------------------------------------------------------
### = NOTES
### =  ( 1): Wrapper of the box model
### =  ( 2): Export the outputs as a data.frame
### =  ( 3): run from restart only works for forward mode
### =  ( 4): Output includes: CH4 emissions, CH4 concentration, delta13CH4 in emissions, delta13CH4 in atmospheric, OH concentration
### =----------------------------------------------------------------------
### = INPUTS
### =  ( 1): vec.St       -- Vector: Time stamp.
### =  ( 2): mat.ch4      -- Matrix: Time series of CH4 Emission (E)
### =  ( 3): vec.OH       -- Vector: Time series of OH concentration (OH)
### =  ( 3): vec.OH       -- Vector: atmospheric CH4 concentration (C)
### =  ( 5): IC           -- Vector: initial conditions for mass of 12CH4_NH,12CH4_SH,13CH4_NH,13CH4_SH,CH3D_NH,CH3d_SH
### =  ( 6): runMode      -- String: inverse mode for calculating OH (E and C are known) or Emissions (OH and C are known)
### =  ( 7): inRestart    -- String: input restart file path (Optional)
### =  ( 8): outRestart   -- String: output retstart file path (Optional)

### =----------------------------------------------------------------------
### = OUTPUTS
### =  ( 1): out -- Data frame with the designated output variables
### =======================================================================

boxModel_wrapper <- function(vec.St,mat.ch4,vec.OH,mat.obs,IC,runMode,...){
  
  ################################################
  isfromrestart <- F
  writeRestart <- F
  #Initial condition: Vector
  params <- getParameters(vec.St)
  
  val_iKIE <- params$iKIE
  
  if( length(list(...)) ){
    Lst <- list(...)
    
    if (!is.null(Lst$inRestart)){
      inRestart <- Lst$inRestart
      isfromrestart <- T
    }
    if( !is.null(Lst$outRestart)){
      outRestart <- Lst$outRestart
      writeRestart <- T
    }
    if( !is.null(Lst$iKIE)){  #if iKIE is not specified, use default from params
      val_iKIE <- Lst$iKIE
    }else{
      val_iKIE <- params$iKIE
    }
  }

  if(isfromrestart){
    # set up the state variables from reading in restart file
    vec.inits <- as.numeric(read.table(inRestart))
    
  }else{

    if(missing(IC)){
      vec.inits <- params$IC  
    }else{
      vec.inits <- IC
    }
  }
  
  
  # Input Matrix
  # - Column 1: 12CH4 in the Northern Hemisphere (mat.S[,1])  Unit: Tg CH4/yr
  # - Column 2: 12CH4 in the Southern Hemisphere (mat.S[,2])  Unit: Tg CH4/yr
  # - Column 3: 13CH4 in the Northern Hemisphere (mat.S[,3])  Unit: per mille
  # - Column 4: 13CH4 in the Southern Hemisphere (mat.S[,4])  Unit: per mille
  # - Column 5: CH3D in the Southern Hemisphere (mat.S[,5])   Unit: per mille
  # - Column 6: CH3d in the Southern Hemisphere (mat.S[,6])   Unit: per mille
  mat.S <- matrix(0,nrow=length(vec.St),ncol=6)
  
  mat.S[,1:4] <-mat.ch4
  
  #CH3D input, to be implemented
  mat.S[,5] <- rep(0,length(vec.St))
  mat.S[,6] <- rep(0,length(vec.St))  
  
  if(runMode=='forward'){
    ################################################
    ### Run the box model
    output.box <- boxModel(St = vec.St,S = mat.S,OH = vec.OH,params = params,obs = mat.obs,inits = vec.inits,iKIE=val_iKIE,runMode = runMode)
    
    #export output data.frame
    df.out <- data.frame(Year=vec.St)
    #CH4 emissions
    df.out[['ECH4_NH']] <- mat.S[,1]
    df.out[['ECH4_SH']] <- mat.S[,2]
    #delta13CH4 in emissions
    df.out[['DCH4E_NH']] <- mat.S[,3]
    df.out[['DCH4E_SH']] <- mat.S[,4]
    #CH4 concentration
    df.out[['CH4_NH']]  <- output.box[,1]+output.box[,3]  #CH4 concentraion in NH
    df.out[['CH4_SH']]  <- output.box[,2]+output.box[,4]  #CH4 concentration in SH
    df.out[['CH4_GL']]  <- rowSums(output.box[,1:4])/2
    #delta13CH4 in atmospheric
    df.out[['DCH4_NH']] <- calC13ratio(output.box[,1],output.box[,3]) #13CH4 signature in NH
    df.out[['DCH4_SH']] <- calC13ratio(output.box[,2],output.box[,4]) #13CH4 signature in SH
    df.out[['DCH4_GL']] <- (calC13ratio(output.box[,1],output.box[,3])+calC13ratio(output.box[,2],output.box[,4]))/2
    #CH4 sink, and OH concentration
    df.out[['CH4_SINK']] <- (rowSums(mat.S[,1:2])/params$mm_ch4 - (df.out$CH4_GL - c(sum(params$IC[1:4])/2, df.out$CH4_GL[1:length(vec.St)-1])) )* params$mm_ch4
    df.out[['OH']] <- vec.OH
    
    if(writeRestart){
      write.table(tail(output.box,1),outRestart)
    }
    return(df.out)
  }else if(runMode=='OH_inverse'){
    out.OH <- boxModel(St = vec.St,S = mat.S,OH = vec.OH,params = params,obs = mat.obs,inits = vec.inits,iKIE=val_iKIE,runMode = runMode)
    return(out.OH) # radicals/cm3 
  }
  

}
