##############################################################################################.
################################### HEADER ###################################################.
####       HELCOM CORE indicator 'Mean Size Total Stock' (MSTS) R-script                  ####.
####       A tool for the assessment of MSTS from ICES ZP community data extract          ####.        
####                                                  or aggregated annual means          ####.
#                   #                   #                   #                     #         ##.
# Developed within the project "Baltic Data Flows" Activity 6: Biodiversity.                ##.
# Project funded by the EU Innovation and Networks Agency (INEA) via CEF funding instrument ##.
# More information on the project: https://balticdataflows.helcom.fi/                       ##.
#                   #                   #                   #                     #         ##.
##############################################################################################.
##############################################################################################.
rm(list=ls())

#####.
# Set working directory where /input and /output folders are located!
#####.
#setwd("C:/ ... ")

#local setwd()
setwd("C:/Users/Astra/Dropbox/PROJEKTI/2020-2023_BalticDataFlows/A6-biodiversity/MSTSautocalc_v2")

## Can be ignored but then a full path has to be specified for input_msts, stations_for_msts,
## addPar_for_msts,species_list_for_msts and path_input_msts 
## (see 'SPECIFY INPUT FILES/VALUES')

###############################################################################################.
### load/install required packages ####
if (!require("data.table")) install.packages("data.table")
if (!require("tidyr")) install.packages("tidyr")
if (!require("magrittr")) install.packages("magrittr")
if (!require("dplyr")) install.packages("dplyr") 
if (!require("ggplot2")) install.packages("ggplot2") # for plots
if (!require("tidyquant")) install.packages("tidyquant") # for plots
if("ggpubr" %in% rownames(installed.packages()) == FALSE) {install.packages("ggpubr")}  ## Kendall's tau via stat_cor()
options(dplyr.summarise.inform = FALSE)
# * - see output of sessionInfo() in the last section!
cat("\nRequired packages added!")

################################################################################################.
### SPECIFY INPUT FILES/VALUES:                                                             ####
################################################################################################.

input_msts <- "./input/ICES_HELCOM_zp.txt"      # ICES zooplankton extract (data from data.ices.dk)
stations_for_msts <- "./input/MSTS_statn.txt"   # Defined stations for reference period (if relevant)
addPar_for_msts <- "./input/MSTS_addPar.txt"    # Thresholds and other infoParameters see input file
path_input_msts <- "./input/MSTS_input_data.txt"  # Input data for GES calculations 
                                                  # (in aggregated form: yearly average)

### !!! is redundant in case when data includes only MSTS-relevant taxa
species_list_for_msts <- "./input/MSTS_ZPspecies.txt" ### Species list

### the beginning of the assessment period
first_year_of_the_assessment_period <- 2016

cat("\nInput files defined!")

###############################################################################################.
#### ZOOPLANKTON HISTORICAL DATA for Ref periods/ threshold/lambda calculations           #####
###############################################################################################.
MSTS_infoPar <- fread(addPar_for_msts, encoding = "UTF-8")

## identify HELCOM_subbasins where historical data are defined in outer source
hist_data_input_Hsubbas <-MSTS_infoPar[MSTS_infoPar$Hist_data_source == "input",]$HELCOM_subbasin
cat(paste0("\nHistorical data are defined in external source ('MSTS_hist_data.txt') for: ", 
           hist_data_input_Hsubbas, "\n"))

## identify HELCOM_subbasins where historical data are from COMBINE extracts
hist_data_COMBINE_Hsubbas <-
  MSTS_infoPar[MSTS_infoPar$Hist_data_source == "COMBINE",]$HELCOM_subbasin
cat(paste0("\nHistorical data are taken from ICES data for: ", 
           hist_data_COMBINE_Hsubbas, "\n"))

## identify HELCOM_subbasins where historical data input is not defined
hist_data_NA_Hsubbas <-MSTS_infoPar[MSTS_infoPar$Hist_data_source == "",]$HELCOM_subbasin
cat(paste0("\nHistorical data input is not defined: ", hist_data_NA_Hsubbas))

## give a message where to define 'input sources'
if(length(hist_data_NA_Hsubbas) > 0)
  {message(paste0("\n --> Please see 'MSTS_addPar.txt' in input folder and fill the empty cells in column Hist_data_source by defining either 'COMBINE' (if historical data should be taken from ICES extract) or 'input' (if defined in 'MSTS_hist_data.txt') to include ", 
                  hist_data_NA_Hsubbas, " in the assessment!\n"))
  }else{ cat("Historical data input sources are defined for all HELCOM subbasins!")
}


##############################################################################################.
### DATA AGGREGATION for MSTS requirements (still in progress)                            ####
##############################################################################################.
## aggregation script is still in progress / 
## when ready should be merged with 'MSTS_GEScalc.R' outcome (see Line 260)
####
## 
source('./addSCRIPTS/ICES_ZP_data_aggregation.R') #output object: MSTS_data


## 'ICES_ZP_data_aggregation.R' aggregates data considering:

###       1. All Baltic Sea subbasins except the ones with unidentified historical data source
###          (see input file 'MSTS_addPar.txt' or object 'hist_data_NA_Hsubbas'in Global 
###          Environment)

###       2. Data obtained following HELCOM COMBINE or earlier methodological recommendations
###          identified by the use of WP-2 net (SMTYP == "WP-2") or noting the use of BMP or  
###          COMB approach in MPROG variable ( grepl("BMP", MPROG) | grepl("COMB", MPROG) )

###       3. Relevant season/period - defined in input file 'MSTS_addPar.txt' (information stored 
###          in variable MSTS_infoPar$Per_start[month] and .$Per_end[month])

###       4. Relevant species - defined in input file 'MSTS_ZPspecies.txt' (1 - relevant; 
###          0 - not used). See object 'species_out' in the Global Environment to see cases
###          filtered out during the aggregation process. THIS WILL BECOME REDUNDANT in case when 
###          ICES data extract will include only MSTS-relevant taxa

###       5. Only sampling events representing the whole defined water column are included.  
###          See defined depth limit for every subbasin in input file 'MSTS_addPar.txt'  
###          information stored in variable MSTS_infoPar$Depth_minLimit[m]. 
###          See object 'layers_out' in the Global Environment to see samples filtered out
###          during the aggregation process.

###       6. Only STNNO that have both parameters reported (ABUNDNR, BMWETWT) are included.
###          See object 'Samples_OnePar_out' in the Global Environment to see samples filtered 
###          out during the aggregation process. 

###       7. SMVOL where reported in different units m3 and L. A threshold of 100 were used to identify
###          reported units. If SMVOL >= 100 then assumed that it is reported in L; if < 100 then - in m3. 
###          SMVOL where then adjusted accordingly to represent m3.

###       8. MUNITS_Final are harmonized to include only mg/m3 and ind/m3. See 'SMVOL_rep_zero'
###          to see cases where SMVOL reported '0' yielding 'Inf' for Fin_value. 'SMVOL_rep_zero'
###          where filtered out.

###       9. Mean from parallel samples (unique layer+Date+STNNO combination) are calculated.
###          Values are calculated to represent the whole water column. Sampling from layers 
###          have been considered and harmonized with samples representing the whole column. 

###        * Final values in the outcome are total MSTS-relevant zooplankton abundance (ind/m3) 
###          and biomass (mg/m3) representing the whole water column sampled (and the depth of the 
###          sampled column is equal or higher than defined in Depth_minLimit[m] for the subbasin).



### save data after aggregation in output folder 
suppressWarnings( dir.create(file.path(getwd(), "output")) )   
suppressWarnings( dir.create(file.path(paste0(getwd(),"/output"), "Aggregated data")) )  

write.csv(MSTS_data, paste0("./output/Aggregated data/MSTS_aggr_output_",Sys.Date(),".txt"),
          row.names=FALSE)
cat(paste0("\nFile: ", paste0("./output/Aggregated data/MSTS_aggr_output_",Sys.Date(),".txt"),
           "   saved in output folder"))

##############################################################################################.
### MSTS Calculations                                                                     ####
##############################################################################################.

### A) REFERENCE PERIODS / THRESHOLD VALUES  ####

## Type I: Historic Data from ICES ZP extract ####
## from the output of 'ICES_ZP_data_aggregation.R'
## Data saved in object 'MSTS_data' in Global Environment.

## Still in progress - some cases fall through layer-maxDepth filter (!!)

#### A-1.1 STATION lIST                    #####
## relevant only if historical data from ICES ZP extract 
## and for threshold value calculations

MSTS_StatnRef <- fread(stations_for_msts, encoding = "UTF-8") %>% dplyr::filter(MSTS_calc == 1)

MSTS_refData_COMBINE_input <- MSTS_data %>%
  dplyr::filter(HELCOM_subbasin %in% hist_data_COMBINE_Hsubbas) %>%
  dplyr::filter(STATN %in% MSTS_StatnRef$STATN) %>%

#####.
  ###  A-1.2 AVERAGING ####
  ###  Monthly averages
  # calculates monthly averages per STATN
  mutate(MONTH = month(Date)) %>%
  group_by(Latitude, Longitude,MYEAR, STATN, HELCOM_subbasin,
          ICES_area, PARAM, MUNIT_msts, Layer_considered, MONTH) %>%
  summarize(month_ave = mean(col_val)) %>%
  
  # monthly average per HELCOM_subbasin
  ungroup() %>% 
  group_by(MYEAR, MONTH, HELCOM_subbasin, ICES_area, PARAM, MUNIT_msts) %>%
  summarize(month_ave_subb = mean(month_ave)) %>% 

  ### Annual averages 
  ungroup() %>% 
  group_by(MYEAR,HELCOM_subbasin, ICES_area, PARAM, MUNIT_msts) %>%
  summarize(annual_ave = mean(month_ave_subb)) %>%
  select(-MUNIT_msts) %>%
  spread(PARAM, annual_ave) %>%
  mutate(MS = BMWETWT*1000/ABUNDNR) %>%
  select(-ABUNDNR)

#####.
  ###  A-1.3 HARMONIZATION ####

## harmonize the output with 'MSTS_input_data.txt' format for merging

MSTS_refData_COMBINE_input <- MSTS_refData_COMBINE_input %>%
  left_join(., MSTS_infoPar[,1:3], by = "HELCOM_subbasin") %>% ungroup()

# identify ref periods for RefCHL and RefFISH per subbasin
s_names <- unique(MSTS_refData_COMBINE_input$HELCOM_subbasin)
refFISH_Comb <- vector("list", length(s_names))
refCHL_Comb <- vector("list", length(s_names))
                                              
for(subb in 1:length(s_names)) {
  refFISH_Comb[[subb]] <- 
    (as.numeric(strsplit(unique(MSTS_refData_COMBINE_input$RefFISH), "[:]")[[subb]][1])):
    (as.numeric(strsplit(unique(MSTS_refData_COMBINE_input$RefFISH), "[:]")[[subb]][2]))
}
names(refFISH_Comb) <- s_names
 
for(subb in 1:length(s_names)) {
  refCHL_Comb[[subb]] <- 
    (as.numeric(strsplit(unique(MSTS_refData_COMBINE_input$RefCHL), "[:]")[[subb]][1])):
    (as.numeric(strsplit(unique(MSTS_refData_COMBINE_input$RefCHL), "[:]")[[subb]][2]))
}
names(refCHL_Comb) <- s_names

# mark reference years with 1 and non-ref years with 0
for(subb in 1:length(s_names)) {
  for(n in 1:nrow(MSTS_refData_COMBINE_input)){
    if(MSTS_refData_COMBINE_input$HELCOM_subbasin[n] == s_names[subb]){
    
      if(MSTS_refData_COMBINE_input$MYEAR[n] %in% c(refCHL_Comb[[subb]]))
      {MSTS_refData_COMBINE_input$RefCHL[n] =1    
      }else{MSTS_refData_COMBINE_input$RefCHL[n]= 0}
}}}
for(subb in 1:length(s_names)) {
  for(n in 1:nrow(MSTS_refData_COMBINE_input)){
    if(MSTS_refData_COMBINE_input$HELCOM_subbasin[n] == s_names[subb]){
    
    if(MSTS_refData_COMBINE_input$MYEAR[n] %in% c(refFISH_Comb[[subb]]))
    {MSTS_refData_COMBINE_input$RefFISH[n] =1    
    }else{MSTS_refData_COMBINE_input$RefFISH[n]= 0}
}}}

# harmonize column names with MSTS_input_data.txt' format for merging

colnames(MSTS_refData_COMBINE_input)[1] <- "year"
colnames(MSTS_refData_COMBINE_input)[4] <- c("TZB")
MSTS_refData_COMBINE_input$RefFISH <- as.numeric(MSTS_refData_COMBINE_input$RefFISH)
MSTS_refData_COMBINE_input$RefCHL<- as.numeric(MSTS_refData_COMBINE_input$RefCHL)
MSTS_refData_COMBINE_input$year<- as.numeric(MSTS_refData_COMBINE_input$year)

#write.csv(MSTS_refData_COMBINE_input, paste0("./output/MSTS_aggr_yearlyRef_output_",Sys.Date(),".txt"))

  ### Type II: Historic Data predefined ####
### Historic Data predefined in: 
MSTS_data_input <- fread(path_input_msts, encoding = "UTF-8") %>%
  dplyr::filter(HELCOM_subbasin %in% hist_data_input_Hsubbas) 
#%>% full_join(., MSTS_refData_COMBINE_input) # when the data gaps will be resolved 
                                              # MSTS_refData_COMBINE_input can be included in further analysis
                                              # for now its commented out

##
suppressWarnings(
  source('./addSCRIPTS/MSTS_GEScalc.R') )

## 'MSTS_GEScalc.R' does MSTS assessment:

###       1. Evaluates Normality using fBasics::dagoTest, shapiro.test, ks.test
###       2. Does Box-Cox transformation if Normality is not met using EnvStats::boxcox
###       3. calculate GES value per subbasin

cat("\n GES values set for ", paste0(MSTS_GES$HELCOM_subbasin, sep = ", "))

suppressWarnings( dir.create(file.path(paste0(getwd(),"/output"), "Fin Asessment Data")) )
write.csv(MSTS_GES, paste0("./output/Fin Asessment Data/MSTS_GES_",Sys.Date(),".txt"), 
          row.names=FALSE)
cat(paste0("\nFile: ", paste0("./output/Fin Asessment Data/MSTS_GES_",Sys.Date(),".txt"),
           "   saved in output folder"))

#####.
### B) MSTS ASSESSMENT                                         ####
#####.

# MSTS_data - an output of 'ICES_ZP_data_aggregation.R' is used as input values
# see description in Lines 92:134 for more information on input object

### Confidence parameters                                      ####
  # Number of observations per year (temporal coverage)
MSTS_conf_ObPerYear <- MSTS_data %>%
  select(HELCOM_subbasin, MYEAR, STATN, Date) %>% unique %>%
  count(HELCOM_subbasin, MYEAR)
colnames(MSTS_conf_ObPerYear)[3] <- "n_ObPerYear"
  # Number of stations observed per year (spatial coverage)
MSTS_conf_STATNobserved <- MSTS_data %>%
  select(HELCOM_subbasin, MYEAR, STATN) %>% unique %>%
  count(HELCOM_subbasin, MYEAR)
colnames(MSTS_conf_STATNobserved)[3] <- "n_STATNobserved"

   # merge into one object:
MSTS_conf <- MSTS_conf_ObPerYear %>%
  full_join(., MSTS_conf_STATNobserved) %>%
  filter(MYEAR >= first_year_of_the_assessment_period)
rm(MSTS_conf_ObPerYear, MSTS_conf_STATNobserved)  

### B0. Averaging                                             ####
  ###  Monthly averages
  # calculates monthly averages per STATN
MSTS_data <- MSTS_data %>%
  mutate(MONTH = month(Date)) %>%
  group_by(Latitude, Longitude,MYEAR, STATN, HELCOM_subbasin,
           ICES_area, PARAM, MUNIT_msts, Layer_considered, MONTH) %>%
  summarize(month_ave = mean(col_val)) %>%

  # monthly average per HELCOM_subbasin
  ungroup() %>% 
  group_by(MYEAR, MONTH, HELCOM_subbasin, ICES_area, PARAM, MUNIT_msts) %>%
  summarize(month_ave_subb = mean(month_ave)) %>% 
  
  ### Annual averages 
  ungroup() %>% 
  group_by(MYEAR,HELCOM_subbasin, ICES_area, PARAM, MUNIT_msts) %>%
  summarize(annual_ave = mean(month_ave_subb)) %>%
  select(-MUNIT_msts) %>%
  spread(PARAM, annual_ave) %>%
  mutate(MS = BMWETWT*1000/ABUNDNR) %>%
  select(-ABUNDNR)

colnames(MSTS_data)[4] <- "TZB"
colnames(MSTS_data)[1] <- "year"

# filter HELCOM_subbasin that have defined GES values
MSTS_assess <- MSTS_data_input %>% 
  dplyr::filter(year < first_year_of_the_assessment_period) %>%
  full_join(., MSTS_data[MSTS_data$year>=first_year_of_the_assessment_period,]) %>%
  dplyr::filter(HELCOM_subbasin %in% MSTS_GES$HELCOM_subbasin) %>%
  select(-Notes, -`Data from`, -Reference) %>%
  mutate(TZB.bc = -1, MS.bc = -1, Z_tzb = 99, Z_ms = 99,  
         LC_tzb = 99, LC_ms = 99, Alarm_tzb = 99, Alarm_ms = 99) %>%
  group_by(HELCOM_subbasin) %>%
  group_split()

names(MSTS_assess) <- MSTS_GES$HELCOM_subbasin
### B1. Data transformation                                 ####
## TZB
for(subb in seq_along(MSTS_assess)){
  if(transf_tzb[subb] == "Passed"){
    MSTS_assess[[subb]]$TZB.bc <- MSTS_assess[[subb]]$TZB
  }else{
    MSTS_assess[[subb]]$TZB.bc <- BCTransform(MSTS_assess[[subb]]$TZB, lambda = l_tzb[[subb]])
  }
}
## MS
for(subb in seq_along(MSTS_assess)){
  if(transf_ms[subb] == "Passed"){
    MSTS_assess[[subb]]$MS.bc <- MSTS_assess[[subb]]$MS
  }else{
    MSTS_assess[[subb]]$MS.bc <- BCTransform(MSTS_assess[[subb]]$MS, lambda = l_ms[[subb]])
  }
}
  

### B2. Z-scores                                            ####
  
# identify which mean and sd (RefCHL or RefFISH) has to be chosen for every subbasin per parameter
MSTS_parGES <- MSTS_ref_outputs %>% dplyr::filter(parameter == "GES")
# TZB
RefP_tzb <- c(rep("", length(MSTS_GES$HELCOM_subbasin)))
names(RefP_tzb) <- MSTS_GES$HELCOM_subbasin

for(subb in 1:length(MSTS_GES$HELCOM_subbasin)) {
  for(n in 1:length(MSTS_parGES$value)){
    if(!is.na(MSTS_parGES$value[n])){
      if(round(MSTS_parGES$value[n], 0) == round(MSTS_GES$GES_tzb[[subb]],0)){
        if(MSTS_parGES$HELCOM_subbasin[n] == MSTS_GES$HELCOM_subbasin[subb]){
          RefP_tzb[[subb]] <- MSTS_parGES$RefP[n]
        }
      } 
    }
  } 
} 

# identify which mean and sd (RefCHL or RefFISH) has to be chosen for every subbasin per parameter
# MS
RefP_ms <- c(rep("", length(MSTS_GES$HELCOM_subbasin)))
names(RefP_ms) <- MSTS_GES$HELCOM_subbasin

for(subb in 1:length(MSTS_GES$HELCOM_subbasin)) {
  for(n in 1:length(MSTS_parGES$value)){
    if(!is.na(MSTS_parGES$value[n])){
      if(round(MSTS_parGES$value[n], 1) == round(MSTS_GES$GES_ms[[subb]],1)){
        if(MSTS_parGES$HELCOM_subbasin[n] == MSTS_GES$HELCOM_subbasin[subb]){
          RefP_ms[[subb]] <- MSTS_parGES$RefP[n]
        }
      } 
    }
  } 
} 

MSTS_Zsc <- MSTS_ref_outputs[MSTS_ref_outputs$RefP == "chl" | MSTS_ref_outputs$RefP == "fish",] %>%
  unite("var", c("variable", "parameter"), sep = "_") %>%
  spread(key = var, value = value) %>%
  group_by(HELCOM_subbasin) %>% group_split()
names(MSTS_Zsc) <- MSTS_GES$HELCOM_subbasin

for(subb in seq_along(MSTS_assess)){
  if(RefP_tzb[[subb]] == "chl"){
    MSTS_assess[[subb]]$Z_tzb <- 
      (MSTS_assess[[subb]]$TZB.bc - MSTS_Zsc[[subb]][MSTS_Zsc[[subb]]$RefP == "chl",]$tzb_mean)/
      MSTS_Zsc[[subb]][MSTS_Zsc[[subb]]$RefP == "chl",]$tzb_sd
  }
  
  if(RefP_tzb[[subb]] == "fish"){
    MSTS_assess[[subb]]$Z_tzb <- 
      (MSTS_assess[[subb]]$TZB.bc - MSTS_Zsc[[subb]][MSTS_Zsc[[subb]]$RefP == "fish",]$tzb_mean)/
      MSTS_Zsc[[subb]][MSTS_Zsc[[subb]]$RefP == "fish",]$tzb_sd
  }
  
  if(RefP_ms[[subb]] == "chl"){
    MSTS_assess[[subb]]$Z_ms <- 
      (MSTS_assess[[subb]]$MS.bc - MSTS_Zsc[[subb]][MSTS_Zsc[[subb]]$RefP == "chl",]$ms_mean)/
      MSTS_Zsc[[subb]][MSTS_Zsc[[subb]]$RefP == "chl",]$ms_sd
  }
  
  if(RefP_ms[[subb]] == "fish"){
    MSTS_assess[[subb]]$Z_ms <- 
      (MSTS_assess[[subb]]$MS.bc - MSTS_Zsc[[subb]][MSTS_Zsc[[subb]]$RefP == "fish",]$ms_mean)/
      MSTS_Zsc[[subb]][MSTS_Zsc[[subb]]$RefP == "fish",]$ms_sd
  }
}

### B3. CuSum values #####

for(subb in seq_along(MSTS_assess)){
  
  for (i in 1:nrow(MSTS_assess[[subb]])) {
    if(i==1){
      MSTS_assess[[subb]]$LC_tzb[i] = min(0, MSTS_assess[[subb]]$Z_tzb[i]+0.5)
      MSTS_assess[[subb]]$LC_ms[i] = min(0, MSTS_assess[[subb]]$Z_ms[i]+0.5)
    
    }else{
      MSTS_assess[[subb]]$LC_tzb[i] = min(0, MSTS_assess[[subb]]$Z_tzb[i]+0.5+MSTS_assess[[subb]]$LC_tzb[i-1])
      MSTS_assess[[subb]]$LC_ms[i] = min(0, MSTS_assess[[subb]]$Z_ms[i]+0.5+MSTS_assess[[subb]]$LC_ms[i-1])
    }
  }
}

### B4. ALARM #####

for(subb in seq_along(MSTS_assess)){
  for(i in 1:nrow(MSTS_assess[[subb]])) {
    if(MSTS_assess[[subb]]$LC_tzb[i] > -5){
      MSTS_assess[[subb]]$Alarm_tzb[i] = 0
    }else{
      MSTS_assess[[subb]]$Alarm_tzb[i] = 1
    }
    
    if(MSTS_assess[[subb]]$LC_ms[i] > -5){
      MSTS_assess[[subb]]$Alarm_ms[i] = 0
    }else{
      MSTS_assess[[subb]]$Alarm_ms[i] = 1
    }
    
  }  
}

rm(list = c("s_names", "refP_fish_GES", "refP_chl_GES", "refFISH_Comb", "refCHL_Comb", 
            "norm_ms", "norm_tzb", "MSTS_parGES"))

### OUTPUTS #####

### Final data set
# transform lists to df
MSTS_assess_df <- data.frame(Reduce(rbind, MSTS_assess))

# save df
write.csv(MSTS_assess_df, paste0("./output/Fin Asessment Data/MSTS_assessment_df_", Sys.Date(),".txt"), 
          row.names=FALSE)
cat(paste0("\nFile: ", paste0("./output/Fin Asessment Data/MSTS_assessment_df_",Sys.Date(),".txt"),
           "   saved in output folder"))
write.csv(MSTS_conf, paste0("./output/Fin Asessment Data/MSTS_confidence_", Sys.Date(),".txt"), 
          row.names=FALSE)
cat(paste0("\nFile: ", paste0("./output/Fin Asessment Data/MSTS_confidence_",Sys.Date(),".txt"),
           "   saved in output folder"))

### Plots
# create folder for Plots
suppressWarnings( dir.create(file.path(paste0(getwd(),"/output"), "Plots")) )
suppressWarnings( dir.create(file.path(paste0(getwd(),"/output/Plots"), "trend_TZB")) )
suppressWarnings( dir.create(file.path(paste0(getwd(),"/output/Plots"), "trend_MS")) )
suppressWarnings( dir.create(file.path(paste0(getwd(),"/output/Plots"), "trend_CuSum")) )
suppressWarnings( dir.create(file.path(paste0(getwd(),"/output/Plots"), "MSTS")) )


### Long-term trends of TZB ####
## scales min-max - for all subbasins equal

plot_list_tzb = list()
for (subb in seq_along(MSTS_assess)){

  p = ggplot(MSTS_assess[[subb]], aes(year, TZB))+
  geom_rect(aes(xmin=first_year_of_the_assessment_period-0.3, 
                xmax=first_year_of_the_assessment_period + 5.5, 
                ymin=-Inf, ymax=MSTS_GES[subb,]$GES_tzb), 
            fill="lightpink", alpha=0.01)+
  geom_rect(aes(xmin=first_year_of_the_assessment_period-0.3, 
                xmax=first_year_of_the_assessment_period + 5.5, 
                ymin=MSTS_GES[subb,]$GES_tzb, ymax=Inf), 
            fill="palegreen", alpha=0.01)+
  geom_hline(aes(yintercept = MSTS_GES[subb,]$GES_tzb),
             linetype="solid", color="red", 
             size=0.8, alpha=0.6)+
  geom_point(size = 2, shape = 19, alpha = 0.7, color = "grey30")+
  tidyquant::geom_ma(ma_fun = SMA, n = 2, linetype = "solid", size = 0.5, color = "grey30")+
  geom_smooth(method = "lm", se=FALSE, linetype = "dotted", colour = "black", size = 0.5)+
  ggpubr::stat_cor(cor.coef.name = "tau", p.accuracy = 0.001)+
  theme_classic(base_size = 12)+
  theme(axis.title.x = element_blank(), 
        plot.title = element_text(hjust = 0.5, face = "plain"))+
  labs(y = expression("Biomass (TZB),"~mg~m^-3), 
      title = paste0(MSTS_GES[subb,]$HELCOM_subbasin))+
  xlim(min(MSTS_assess_df$year), first_year_of_the_assessment_period + 5.5) + 
  ylim(min(MSTS_assess_df$TZB), max(MSTS_assess_df$TZB))
  
  plot_list_tzb[[subb]] <- p
}

# Save plots to tiff. Makes a separate file for each plot.
for (subb in seq_along(MSTS_assess)) {
  file_name = paste("./output/Plots/trend_TZB/trend_TZB_", MSTS_GES[subb,]$HELCOM_subbasin, ".tiff", sep="")
  tiff(file_name, width = 15, height = 10, units = "cm", res=300, compression = 'lzw')
  print(plot_list_tzb[[subb]])
  dev.off()
}

# Another option: create pdf where each page is a separate plot.
pdf("./output/Plots/trend_TZB/trend_TZB.pdf")
for (subb in seq_along(MSTS_assess)) {
    print(plot_list_tzb[[subb]])
}
dev.off()

## Y autoscale - Y scale different for subbasins

plot_list_tzb = list()
for (subb in seq_along(MSTS_assess)){
  
  p = ggplot(MSTS_assess[[subb]], aes(year, TZB))+
    geom_rect(aes(xmin=first_year_of_the_assessment_period-0.3, 
                  xmax=first_year_of_the_assessment_period + 5.5, 
                  ymin=-Inf, ymax=MSTS_GES[subb,]$GES_tzb), 
              fill="lightpink", alpha=0.01)+
    geom_rect(aes(xmin=first_year_of_the_assessment_period-0.3, 
                  xmax=first_year_of_the_assessment_period + 5.5, 
                  ymin=MSTS_GES[subb,]$GES_tzb, ymax=Inf), 
              fill="palegreen", alpha=0.01)+
    geom_hline(aes(yintercept = MSTS_GES[subb,]$GES_tzb),
               linetype="solid", color="red", 
               size=0.8, alpha=0.6)+
    geom_point(size = 2, shape = 19, alpha = 0.7, color = "grey30")+
    tidyquant::geom_ma(ma_fun = SMA, n = 2, linetype = "solid", size = 0.5, color = "grey30")+
    geom_smooth(method = "lm", se=FALSE, linetype = "dotted", colour = "black", size = 0.5)+
    ggpubr::stat_cor(cor.coef.name = "tau", p.accuracy = 0.001)+
    theme_classic(base_size = 12)+
    theme(axis.title.x = element_blank(), 
          plot.title = element_text(hjust = 0.5, face = "plain"))+
    labs(y = expression("Biomass (TZB),"~mg~m^-3), 
         title = paste0(MSTS_GES[subb,]$HELCOM_subbasin))+
    xlim(min(MSTS_assess_df$year), first_year_of_the_assessment_period + 5.5)
  
  plot_list_tzb[[subb]] <- p
}

# Save plots to tiff. Makes a separate file for each plot.
for (subb in seq_along(MSTS_assess)) {
  file_name = paste("./output/Plots/trend_TZB/trend_TZB_autoY_", MSTS_GES[subb,]$HELCOM_subbasin, ".tiff", sep="")
  tiff(file_name, width = 15, height = 10, units = "cm", res=300, compression = 'lzw')
  print(plot_list_tzb[[subb]])
  dev.off()
}

# Another option: create pdf where each page is a separate plot.
pdf("./output/Plots/trend_TZB/trend_TZB_autoY.pdf")
for (subb in seq_along(MSTS_assess)) {
  print(plot_list_tzb[[subb]])
}
dev.off()


### Long-term trends of MS ####
## scales min-max - for all subbasins equal
plot_list_ms = list()
for (subb in seq_along(MSTS_assess)){
  
  p = ggplot(MSTS_assess[[subb]], aes(year, MS))+
    geom_rect(aes(xmin=first_year_of_the_assessment_period-0.3, 
                  xmax=first_year_of_the_assessment_period + 5.5, 
                  ymin=-Inf, ymax=MSTS_GES[subb,]$GES_ms), 
              fill="lightpink", alpha=0.01)+
    geom_rect(aes(xmin=first_year_of_the_assessment_period-0.3, 
                  xmax=first_year_of_the_assessment_period + 5.5, 
                  ymin=MSTS_GES[subb,]$GES_ms, ymax=Inf), 
              fill="palegreen", alpha=0.01)+
    geom_hline(aes(yintercept = MSTS_GES[subb,]$GES_ms),
               linetype="solid", color="red", 
               size=0.8, alpha=0.6)+
    geom_point(size = 2, shape = 19, alpha = 0.7, color = "grey30")+
    tidyquant::geom_ma(ma_fun = SMA, n = 2, linetype = "solid", size = 0.5, color = "grey30")+
    geom_smooth(method = "lm", se=FALSE, linetype = "dotted", colour = "black", size = 0.5)+
    ggpubr::stat_cor(cor.coef.name = "tau", p.accuracy = 0.001)+
    theme_classic(base_size = 12)+
    theme(axis.title.x = element_blank(), 
          plot.title = element_text(hjust = 0.5, face = "plain"))+
    labs(y = expression("Mean Size (MS), \U00B5g"~ind^-1), 
          title = paste0(MSTS_GES[subb,]$HELCOM_subbasin))+
    xlim(min(MSTS_assess_df$year), first_year_of_the_assessment_period + 5.5) + 
    ylim(min(MSTS_assess_df$MS), max(MSTS_assess_df$MS))
  
  plot_list_ms[[subb]] <- p
}

# Save plots to tiff. Makes a separate file for each plot.
for (subb in seq_along(MSTS_assess)) {
  file_name = paste("./output/Plots/trend_MS/trend_MS_", MSTS_GES[subb,]$HELCOM_subbasin, ".tiff", sep="")
  tiff(file_name, width = 15, height = 10, units = "cm", res=300, compression = 'lzw')
  print(plot_list_ms[[subb]])
  dev.off()
}

# Another option: create pdf where each page is a separate plot.
pdf("./output/Plots/trend_MS/trend_MS.pdf")
for (subb in seq_along(MSTS_assess)) {
  print(plot_list_ms[[subb]])
}
dev.off()

## Y autoscale - Y scale different for subbasins
plot_list_ms = list()
for (subb in seq_along(MSTS_assess)){
  
  p = ggplot(MSTS_assess[[subb]], aes(year, MS))+
    geom_rect(aes(xmin=first_year_of_the_assessment_period-0.3, 
                  xmax=first_year_of_the_assessment_period + 5.5, 
                  ymin=-Inf, ymax=MSTS_GES[subb,]$GES_ms), 
              fill="lightpink", alpha=0.01)+
    geom_rect(aes(xmin=first_year_of_the_assessment_period-0.3, 
                  xmax=first_year_of_the_assessment_period + 5.5, 
                  ymin=MSTS_GES[subb,]$GES_ms, ymax=Inf), 
              fill="palegreen", alpha=0.01)+
    geom_hline(aes(yintercept = MSTS_GES[subb,]$GES_ms),
               linetype="solid", color="red", 
               size=0.8, alpha=0.6)+
    geom_point(size = 2, shape = 19, alpha = 0.7, color = "grey30")+
    tidyquant::geom_ma(ma_fun = SMA, n = 2, linetype = "solid", size = 0.5, color = "grey30")+
    geom_smooth(method = "lm", se=FALSE, linetype = "dotted", colour = "black", size = 0.5)+
    ggpubr::stat_cor(cor.coef.name = "tau", p.accuracy = 0.001)+
    theme_classic(base_size = 12)+
    theme(axis.title.x = element_blank(), 
          plot.title = element_text(hjust = 0.5, face = "plain"))+
    labs(y = expression("Mean Size (MS), \U00B5g"~ind^-1), 
         title = paste0(MSTS_GES[subb,]$HELCOM_subbasin))+
    xlim(min(MSTS_assess_df$year), first_year_of_the_assessment_period + 5.5)
  
  plot_list_ms[[subb]] <- p
}

# Save plots to tiff. Makes a separate file for each plot.
for (subb in seq_along(MSTS_assess)) {
  file_name = paste("./output/Plots/trend_MS/trend_MS_autoY_", MSTS_GES[subb,]$HELCOM_subbasin, ".tiff", sep="")
  tiff(file_name, width = 15, height = 10, units = "cm", res=300, compression = 'lzw')
  print(plot_list_ms[[subb]])
  dev.off()
}

# Another option: create pdf where each page is a separate plot.
pdf("./output/Plots/trend_MS/trend_MS_autoY.pdf")
for (subb in seq_along(MSTS_assess)) {
  print(plot_list_ms[[subb]])
}
dev.off()

### CuSum trend ####
MSTS_cusum_df <- MSTS_assess_df %>% gather(cusum, value, LC_tzb:LC_ms) 
MSTS_cusum <- MSTS_cusum_df %>%
  mutate(cusum = factor(cusum, 
                        levels = c("LC_ms", "LC_tzb"), 
                        labels = c("Lower CuSum for MS", 
                                   "Lower CuSum for TZB"))) %>%
  group_by(HELCOM_subbasin) %>%
  group_split()

plot_list_cusum = list()
for (subb in seq_along(MSTS_cusum)){
  
  p = ggplot(MSTS_cusum[[subb]], aes(x=year, value))+
      geom_rect(aes(xmin=first_year_of_the_assessment_period-0.3, 
                    xmax=first_year_of_the_assessment_period+ 5.5, 
                    ymin=-5, ymax=-Inf), fill="lightpink", alpha=0.01)+
      geom_rect(aes(xmin=first_year_of_the_assessment_period-0.3, 
                    xmax=first_year_of_the_assessment_period+ 5.5, 
                    ymin=Inf, ymax=-5), fill="palegreen", alpha=0.01)+
      geom_hline(yintercept = 0, colour="grey20", size=0.3, alpha = 0.2)+
      geom_hline(yintercept = -5, colour="red", size=0.8, alpha = 0.6)+
      scale_colour_manual(values=c("black", "grey50", "grey80"))+
      geom_line(aes(colour=cusum, linetype=cusum), size=0.5)+
      geom_point(aes(fill=cusum, shape=cusum), size=2)+
      scale_linetype_manual(values=c("solid", "solid", "solid"))+
      scale_fill_manual(values=c("black", "white", "black"))+
      scale_shape_manual(values=c(21, 21, 4))+
      theme_classic(base_size = 12)+
      theme(panel.grid.minor = element_blank(),
            panel.grid.major = element_blank(),
            legend.position = "top",
            legend.title = element_blank(),
            plot.margin=unit(c(0.5,0.5,0,0),"cm"),
            axis.title.x = element_blank(),
            title = element_text(face = "plain", hjust = 0.5))+
      labs(y="Lower CuSum",
           title = paste0(MSTS_GES[subb,]$HELCOM_subbasin))+
    xlim(min(MSTS_cusum_df$year), first_year_of_the_assessment_period + 5.5) + 
    ylim(min(MSTS_cusum_df$value), 0.3)
  
  plot_list_cusum[[subb]] <- p
  
}  
      
# Save plots to tiff. Makes a separate file for each plot.
for (subb in seq_along(MSTS_cusum)) {
  file_name = paste("./output/Plots/trend_CuSum/trend_cusum_", MSTS_GES[subb,]$HELCOM_subbasin, ".tiff", sep="")
  tiff(file_name, width = 15, height = 10, units = "cm", res=300, compression = 'lzw')
  print(plot_list_cusum[[subb]])
  dev.off()
}

# Another option: create pdf where each page is a separate plot.
pdf("./output/Plots/trend_CuSum/trend_cusum.pdf")
for (subb in seq_along(MSTS_cusum)) {
  print(plot_list_cusum[[subb]])
}
dev.off()

rm(list = c("MSTS_cusum", "MSTS_cusum_df"))

### MSTS plot ####
plot_list_msts = list()
for (subb in seq_along(MSTS_assess)){
  
     p =  ggplot(MSTS_assess[[subb]], aes(x=TZB, y=MS,
                                   shape = year<first_year_of_the_assessment_period,
                                   color = year<first_year_of_the_assessment_period))+
        geom_point(size=3, alpha=0.75)+
        scale_color_manual(values=c("black", "grey50"))+
        scale_shape_manual(values=c(19, 46))+
        geom_text(aes(label=year), hjust=0.5, vjust=-0.5, 
                  size=3, show.legend = FALSE)+
        geom_hline(aes(yintercept = MSTS_GES[subb,]$GES_ms),
                   linetype="dashed", color="green4", 
                   size=0.8, alpha=0.5)+
        geom_vline(aes(xintercept = MSTS_GES[subb,]$GES_tzb), 
                   linetype="dashed", color="green4",
                   size=0.8, alpha=0.5)+
        theme_light(base_size = 12)+
        theme(panel.grid.minor = element_blank(),
              panel.grid.major = element_blank(),
              legend.position = "none",
              legend.title=element_blank(),
              plot.title = element_text(face = "plain", hjust = 0.5))+
        labs(
          x=expression("TZB,"~mg~m^-3),
          y=expression("MS,"~"\U00B5g"~ind^-1),
          title = paste0(MSTS_GES[subb,]$HELCOM_subbasin))

     plot_list_msts[[subb]] <- p
}

# Save plots to tiff. Makes a separate file for each plot.
for (subb in seq_along(MSTS_assess)) {
  file_name = paste("./output/Plots/MSTS/MSTS_", MSTS_GES[subb,]$HELCOM_subbasin, ".tiff", sep="")
  tiff(file_name, width = 15, height = 10, units = "cm", res=300, compression = 'lzw')
  print(plot_list_msts[[subb]])
  dev.off()
}

# Another option: create pdf where each page is a separate plot.
pdf("./output/Plots/MSTS/MSTS.pdf")
for (subb in seq_along(MSTS_assess)) {
  print(plot_list_msts[[subb]])
}
dev.off()

#### SessionInfo()  ##############################################################################
  
# R version 4.0.3 (2020-10-10)
# Platform: x86_64-w64-mingw32/x64 (64-bit)
# Running under: Windows 10 x64 (build 19042)
# 
# Matrix products: default
# 
# locale:
# [1] LC_COLLATE=Latvian_Latvia.1257  LC_CTYPE=Latvian_Latvia.1257   
# [3] LC_MONETARY=Latvian_Latvia.1257 LC_NUMERIC=C                   
# [5] LC_TIME=Latvian_Latvia.1257    
# system code page: 1252
# 
# attached base packages:
# [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
# [1] ggrepel_0.9.1              tidyquant_1.0.3           
# [3] quantmod_0.4.18            TTR_0.24.2                
# [5] PerformanceAnalytics_2.0.4 xts_0.12.1                
# [7] zoo_1.8-8                  lubridate_1.7.10          
# [9] ggplot2_3.3.5              dplyr_1.0.7               
# [11] magrittr_2.0.1            tidyr_1.1.4               
# [13] data.table_1.14.2         
# 
#####################################################################################################.
