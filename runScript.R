### =======================================================================
### = runScript.R
### = Zhen Zhang
### = Created on 09/30/2019
### =----------------------------------------------------------------------
### = NOTES:
### = 
### = This is the  script for running the 2-box model methane in forward mode
### = and for plotting the time series of CH4 concentration and atmospheric delta13CH4
### = for the emission scenarios in Zhang et al., (2020).
### = The emission scenarios were derived from combining different bottom-
### = up CH4 estimates. For more details, check Supporting Information in 
### = Zhang et al., (2020)
### = The inputs are stored as Rdata in 'data/emission_scenarios'. The description
### = about strcture of Rdata inputs can be found in README.md
### = The plot of time series for emission scenarios are stored in 'plots/'
### =======================================================================
          
###Set up the worksapce to the source code folder
#setwd(".../BoxModel_R_20200520")
setwd("/Users/zzhang88/Research/Github/BoxModel_R_20200520")

##Load source code
library(foreach)
library(zoo)
#load  code and parameters
source('model/tools.R')
source('model/boxModel.R')
source('model/getParameter.R')
source('model/boxModel_wrapper.R')


################################################

#Set up Box model inputs
startyear <- 1980
endyear <- 2017
nyear <- endyear - startyear + 1
vec.St <- seq.Date(from=as.Date(paste(startyear,'-01-01',sep='')),to=as.Date(paste(endyear,'-01-01',sep='')),by='years')
path.in <- "data/emission_scenarios/"
if.varOH <- TRUE

#Read in Observations
load('data/CH4_obs_boxmodel_1980-2017.RData')
################################################################################
################################################################################
file.list <- dir(path.in)
for(i in 1:length(file.list)){  # Index of Emission scenarios
  #Read emission scenario in Rdata
  load(paste(path.in,file.list[i],sep=''))
  
  #matrix of input, CH4 emissions and delta13CH4
  mat.ch4 <- matrix(0,nrow=length(vec.St),ncol=4)
  mat.ch4[,1] <- list.es$EM_NH$Total - list.es$SU_NH
  mat.ch4[,2] <- list.es$EM_SH$Total - list.es$SU_SH
  mat.ch4[,3] <- list.es$EM_D13_NH$Total
  mat.ch4[,4] <- list.es$EM_D13_NH$Total
  
  #######################################################
  #Run box model in the forward mode
  #outputs is a data.frame
  df.out <- boxModel_wrapper(vec.St = vec.St,mat.ch4 = mat.ch4,vec.OH = list.es$OH,
                             mat.obs = df.obs_ch4[1:nyear,c(2,3,8)],IC = IC, runMode='forward',
                             inRestart='data/spinup_restart.txt',iKIE=1/EpsilontoKIE(list.es$EPSILON))
  
  #run box model with the 1000 sets of delta13C-CH4 inputs generated from Monte Carlo
  list.df.out <- foreach(irun=1:1000) %dopar%{
    mat.ch4[,1] <- list.es$EM_NH$Total - list.es$SU_NH
    mat.ch4[,2] <- list.es$EM_SH$Total - list.es$SU_SH
    mat.ch4[,3] <- list.es$EM_D13_MC_NH[[irun]]$Total
    mat.ch4[,4] <- list.es$EM_D13_MC_SH[[irun]]$Total
    
    df.out <- boxModel_wrapper(vec.St = vec.St,mat.ch4 = mat.ch4,vec.OH = list.es$OH,
                               mat.obs = df.obs_ch4[1:nyear,c(2,3,8)],IC = IC, runMode='forward',
                               inRestart='data/spinup_restart.txt',iKIE=1/EpsilontoKIE(list.es$EPSILON))
    return(df.out)
  }
  
  #######################################################
  #Plotting the global time seris in CH4 under varying OH sceanrio
  old.par <- par( no.readonly = TRUE )
  pdf(paste('plots/',file.list[i],'.pdf',sep=''),width=5,height=8)
  par(mfcol=c(2,1),oma=c(0,0,2,0))
  par(mar = c(2,4,2,1)+0.1)
  
  plot(df.obs_ch4$Year,df.obs_ch4$CH4_GL,
       type='p',xlab='Year',ylab='CH4 conc (ppbv)',ylim=c(1500,2000),pch=3,cex=0.5)
  title("CH4 concentration",line=-1)
  lines(startyear:endyear,df.out$CH4_GL,type='l',cex=0.5)
  legend('bottomright',col=c('black'),lty=c(NA,1),c('obs',"model"), bty='n',
         pch=c(3,NA),cex=0.7,x.intersp = 1,y.intersp=0.8,inset=c(0.04,0.01))
  
  #Sanity check of CH4 concentration for the 1000 sets of Monte Carlo runs
  for(irun in 1:1000){
    lines(startyear:endyear,list.df.out[[irun]]$CH4_GL,type='l',cex=0.5,col='black')
  }
  
  #Plotting the global time seris in CH4 and C13-CH4 under varying OH sceanrio
  plot(df.obs_ch4$Year,df.obs_ch4$DCH4_GL,
       type='p',xlab='',ylab='delta13C-CH4 (per mille)',ylim=c(-48.0,-46.0),pch=3,cex=0.5)
  title("delta13C-CH4 ",line=-1)
  lines(startyear:endyear,list.es$OUT_VAR_OH$DCH4_GL,type='l',cex=0.5)
  legend('bottomright',col=c('black'),lty=c(NA,1),c('obs',"model"), bty='n',
         pch=c(3,NA),cex=0.7,x.intersp = 1,y.intersp=0.8,inset=c(0.04,0.01))
  
  for(irun in 1:1000){
    lines(startyear:endyear,list.df.out[[irun]]$DCH4_GL,type='l',cex=0.5,col='#D9D9D996')
  }
  lines(startyear:endyear,list.es$OUT_VAR_OH$DCH4_GL,col=1,lwd=1.5)
  
  mtext(paste("ES: ",strsplit(file.list[i],split = ".Rdata"),sep=''),outer=T,cex=1,line=-1)
  dev.off()
  par(old.par)
  #######################################################
  print(i)
  rm(list.es)
}
