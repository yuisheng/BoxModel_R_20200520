# BoxModel_R_20200520
The two-box atmospheric model from Zhang et al. (2020)

Zhen Zhang
May 20, 2020


# Info
 * The corresponding paper is Zhang et al., (2020): "Domintant contribution of anthropogenic emissions to the rise of atmospheric methane".
 * This box model code is a rewrite of Matlab code from Turner et al., 2017: https://doi.org/10.1073/pnas.1616020114.
 * An example of running the box model can be found in the "runScript.R" file. The script also contains code for ploting results as shown in 'plots/' folder 
 * The bottom-up estimates use publicly available datasets (details can be found in the Data Availability of the paper). As such, we have not included those datasets in the repository.
 * Users could manually download the data by contacting the PIs for those datasets.
 * The inputs and outputs of the box model results are provided in the 'data' folder in the format of R data.
 * Individual Rdata file in the 'data' folder represents one emission scenario, which contains the inputs of
 * CH4 emissions, global weighted-average delta13-CH4 signature in source, and OH time series derived from inverse mode.
 * Emission scenario Rdata naming:
 *   - Industrial fossil Fuel: IFF1:EDGARv4.2; IFF2:EDGAR4.3.2; IFF3:GAINS; IFF4: Schewitzke et al., (2016)
 *   - Agricutural and waste: AGW1:EDGARv4.2; AGW2: EDGAR4.3.2; AGW3: Wolf et al., (2017)
 *   - Wetlands: WET1: REN; WET2: CRU
 *   - Biomass Burning: BB1: GFED; BB2: Worden et al., (2017)
 *   - Geologic: GEO1: High geologic; GEO2: Low geologic
 * Rdata structure (type:description):
 * ├── EM_D13_MC_NH (list: 1000 sets of Monte Carlo EM_D13_NH [per mille])
 * ├── EM_D13_MC_SH (list: 1000 sets of Monte Carlo EM_D13_SH [per mille])
 * ├── EM_D13_NH (data frame: global average delta13-CH4 in source from Northern Hemispheresoil [per mille])
 * ├── EM_D13_SH (data frame: global average delta13-CH4 in source from Southern Hemisphere [per mille])
 * ├── EM_NH (data frame: CH4 emissions from Northern Hemisphere [Tg CH4 yr-1])
 * ├── EM_SH (data frame: CH4 emissions from Southern Hemisphere [Tg CH4 yr-1])
 * ├── EPSILON (kinetic isotope effect value)
 * ├── OH (vector: varying OH level derived from running box model in inverse mode)
 * ├── OUT_CON_OH (data frame: outputs from box model run assuming constant OH since 1992)
 * └── OUT_VAR_OH (data frame: outputs from box model run using varying OH)
 * ├── OUT_MC_CON_OH (list:outputs from Monte Carlo run assuming constant OH since 1992)
 * ├── OUT_MC_VAR_OH (list:outputs from Monte Carlo box model run using varying OH )
 * ├── SU_NH (vector: climatology of soil uptake in Northern Hemisphere [Tg CH4 yr-1])
 * ├── SU_SH (vector: climatology of soil uptake in Southern Hemisphere [Tg CH4 yr-1])
 