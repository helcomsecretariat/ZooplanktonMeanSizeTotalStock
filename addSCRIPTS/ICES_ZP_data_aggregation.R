# source file for 'MAIN_MSTS_calc.R'
### Input files have been specified in 'MAIN_MSTS_calc.R'

##############################################################################################.
###         DATA AGGREGATION                                                              ####
##############################################################################################.

### 1. DATA FILTERING ####
### pick out data relevant for MSTS calculation which Final_value is not NULL
MSTS_data <- fread(input_msts, encoding = "UTF-8") %>% 

  #recalculate Value to get Final_value, considering also SMVOL units 
  # Assumption: if SMVOL > 100 then it is reported as L, if SMVOL < 100 then it is m3
  mutate(SMVOL_adj = case_when(SMVOL <= 100 ~ SMVOL,
                               SMVOL > 100 ~ SMVOL/1000)) %>%
  mutate(Fin_value = case_when(MUNIT == "nr/m3" ~ Value
                                , MUNIT == "mg/m3" ~ Value
                                , MUNIT == "g" ~ (Value*1000)*(NPORT/(CPPORT*SMVOL_adj))
                                , MUNIT == "nr" ~ Value*(NPORT/(CPPORT*SMVOL_adj))
                                , MUNIT == "ug/m3" ~ Value/1000),
         
         MUNIT_msts = case_when(MUNIT == "nr/m3" ~ MUNIT
                                , MUNIT == "mg/m3" ~ MUNIT
                                , MUNIT == "g" ~ "mg/m3"
                                , MUNIT == "nr" ~ "nr/m3"
                                , MUNIT == "ug/m3" ~ "mg/m3"))

SMVOL_rep_zero <- MSTS_data[MSTS_data$Fin_value == Inf,] 
  
#filter out cases where SMVOL reported zero yielding Inf as Fin_value
MSTS_data <- MSTS_data %>%
  dplyr::filter(HELCOM_subbasin %in% hist_data_COMBINE_Hsubbas |
          HELCOM_subbasin %in% hist_data_input_Hsubbas ) %>%
  dplyr::filter(!Fin_value == Inf)%>%
  #####.
  ### 1.1. Sampling method ####
### Filter only data obtained using WP-2 net (SMTYP == WP-2) or MPROG reported as BMP or COMB
  mutate(COMB_method = 
         case_when(SMTYP == "WP-2" | grepl("BMP", MPROG) | grepl("COMB", MPROG) ~ 1,
                               TRUE ~ 0)) %>%
  #####.    
  ### 1.2. Relevant period/season ####
### filter only period relevant for MSTS as defined in input file 'MSTS_addPar.txt'
### information stored in MSTS_infoPar$Per_start[month] and .$Per_end[month]
  full_join(MSTS_infoPar[,c("HELCOM_subbasin","Per_start[month]","Per_end[month]", 
                            "Depth_minLimit[m]",  "Depth_maxLimit[m]")], 
          by = "HELCOM_subbasin") %>%
  mutate(MONTH = month(Date)) %>%
  mutate(period = 
           case_when(MONTH >= `Per_start[month]` & MONTH <= `Per_end[month]` ~ "include")) %>%
  dplyr::filter(period == "include") %>%

#############################################################################################.
### !!! This section is redundant in case when data includes only MSTS-relevant taxa    ##
############################################################################################.
### SPECIES LIST - filtering only relevant species defined as '1' in MSTS_SpeciesList$MSTS_use

#####.
### 1.3. Relevant species ####
### merge information about MSTS-relevant species from MSTS_SpeciesList
  left_join(., fread(species_list_for_msts, encoding = "UTF-8")
            [, c("SPECI","MSTS_use")], by="SPECI")

#save cases left out
species_out <- MSTS_data %>% dplyr::filter (MSTS_use == 0)

# filter cases with MSTS-relevant species/taxa
MSTS_data <- MSTS_data %>% dplyr::filter(MSTS_use == 1)

##########################################################################################.
##########################################################################################.

#####.
### 1.4. Layers/whole column ####
### Sampling from layers (not the whole column)
## identify information on the sampled layer and compare to depth considered for the station (see 'MSTS_StatnRef.txt')
b <- MSTS_data %>%
  select(HELCOM_subbasin, Date,STATN, Latitude, Longitude, STNNO, SMPNO, 
         `Depth_minLimit[m]`, `Depth_maxLimit[m]`,MNDEP, MXDEP) %>%
  # return all sampled layers for every unique sampling event
  # create groups for splitting dataset into every unique station visit (including all layers sampled)
  group_by(HELCOM_subbasin, Date,Latitude, Longitude, STNNO,  `Depth_minLimit[m]`, `Depth_maxLimit[m]`) %>%
  summarise(layer_from = MNDEP, layer_to = MXDEP, SMPNO = SMPNO, STATN = STATN) %>%
  # create variable that will hold the information whether the layer should be included or not
  # value is set to 99 (to be different from the resulting values [i.e. 1 or] and to set variable to numeric) 
  unique() %>% mutate(layer_use = 99) %>%
  # split df into list by groups
  group_split()

## Loop looks for samples that do not represent the considered water column (see 'MSTS_StatnRef.txt'), e.g. are deeper
### Variation of sampled depth vs defined `BotDepth_considered[m]` depth has been considered (-30% and +20% variation)
for(n in seq_along(b)){
  
  for(L in 1:nrow(b[[n]])){
    if(b[[n]]$layer_to[L] %in% b[[n]]$layer_from){
      b[[n]]$layer_use[L] = 1
    }else if(as.numeric(b[[n]]$layer_to[L]) >= b[[n]]$`Depth_minLimit[m]`[L])
    {b[[n]]$layer_use[L] = 1
    }else{b[[n]]$layer_use[L] = 0}
  }
  
  ##
  #### look for cases when the maximal depth was smaller then defined `Depth_minLimit[m]`
  #### in order to be able to remve all sampled layers for the STNNO from further calculations for these cases  
  if(max(as.numeric(b[[n]]$layer_to)) < b[[n]]$`Depth_minLimit[m]`[L]){
    b[[n]]$layer_use = 0
  }
  if(!is.na(b[[n]]$`Depth_maxLimit[m]`[L])){
    if(as.numeric(b[[n]]$layer_to[L]) > b[[n]]$`Depth_maxLimit[m]`[L]){
    b[[n]]$layer_use[L] = 0
  }}
} 
## .$layer_use == 1 if the sampling event covers all of the considered water column and 
### have information for all successive layers from top to bottom

## .$layer_use == 0 if the sample is collected deeper (relevant for Landsort Deep where top 30 meters are considered)
### samples representing 30-60 meters are marked with 0
### also samples that do not have the next successive layer sampled 
### (e.g. samples 0-10, 10-41, 10-60, 60-200: sample from layer 10-41 is excluded .$out == 0)

######################################################################################.
## make df from the list!
b <- data.frame(Reduce(rbind, b)) %>%
  select(-layer_from, -layer_to, -`Depth_minLimit.m.`)

### merge with MSTS_data
MSTS_data <- MSTS_data %>%
  left_join(., b, by = c("HELCOM_subbasin", "Date", "Latitude", "Longitude", "STNNO", "SMPNO", "STATN"))

### remove b, MSTS_data_lay, n, L from the Global Environment. 
### Information needed has been merged into MSTS_data
rm(b, n, L)

### save information about the layers/samples left out after layer-selection
layers_out <- MSTS_data %>%
  select(HELCOM_subbasin, Date, STATN, STNNO, SMPNO, Latitude, Longitude, MNDEP, MXDEP, layer_use) %>%
  dplyr::filter(layer_use == 0) %>% 
  group_by(HELCOM_subbasin, Date, STATN, STNNO, SMPNO, Latitude, Longitude, MNDEP, MXDEP) %>%
  unique()

#######################################################################.
### DATA AVERAGING ####
######################################################################.
### 2. CALCULATIONS (TZA, TZB, MS) ####
### calculate zooplankton parameters for the whole column
#####.
### 2.1. Final MUNITS ####

###check if Final_MUNIT == mg/m3 and ind/m3
if(all(unique(MSTS_data$MUNIT_msts) %in% c("mg/m3", "nr/m3")) == TRUE){
  cat("\nAll Final units are mg/m3 and nr/m3, no unit transformation is needed!")
}else{
  stop("\nMSTS-message: Final units .$MUNIT_msts need to be checked and transformed
        to include only mg/m3 (for PARAM == BMWETWT) and nr/m3 (for PARAM == ABUNDNR)")
}

### Filter only cases which correspond to considered depth and 
### (if sampled in layers) have consecutive layers representing the whole column
MSTS_data <- MSTS_data %>% dplyr::filter(layer_use == 1) %>%

#####.
### 2.2.  Calculations ####
## calculate sum for species (summing all size groups/stages together to represent one species/taxa)
group_by(HELCOM_subbasin, CNTRY,CRUIS,STNNO, Latitude, Longitude,MYEAR, Date,MNDEP,MXDEP,STATN,
         ICES_area, PARAM, SPECI_name, MUNIT_msts) %>%
  summarise(speci_sum = sum(Fin_value)) %>%
  
  ## calculate mean value if there are two parallel samples for the same layer+Date+STNNO
  group_by(HELCOM_subbasin, CNTRY,CRUIS,STNNO, Latitude, Longitude,MYEAR, Date,MNDEP,MXDEP,STATN,
           ICES_area, PARAM, SPECI_name, MUNIT_msts) %>%
  summarise(speci_sum_mean = mean(speci_sum)) %>%
  
  ## calculate total zooplankton abundance and biomass per sampled layer
  group_by(HELCOM_subbasin, CNTRY,CRUIS,STNNO, Latitude, Longitude,MYEAR, Date,MNDEP,MXDEP,STATN,
           ICES_area, PARAM, MUNIT_msts) %>%
  summarise(Value_msts_sum = sum(speci_sum_mean))

#####.
### 2.3. Presence of both parameters #####
####################################################################################.
### !!! This section  is redundant in case when:                                  ##
### !!!               ICES calculate BMWETWT centralized for all reported ABUNDNR ##
### !!!               and no cases of only ABUNDNR reported are present           ##
####################################################################################.
### identify reported cases when only one of two 
## required parameters ("BMWETWT", "ABUNDNR") is reported 
MSTS_data_onePARAM <- MSTS_data %>%
  ungroup() %>% select(HELCOM_subbasin, CRUIS, Date,STATN,STNNO, Latitude, Longitude, PARAM) %>%
  unique() %>% count(HELCOM_subbasin, CRUIS, STATN, STNNO, Latitude, Longitude, Date)

Samples_OnePar_out <- MSTS_data_onePARAM %>% filter(n == 1)

### remove cases when only one of two 
## required parameters ("BMWETWT", "ABUNDNR") is reported 
MSTS_data <- MSTS_data %>%
  left_join(., MSTS_data_onePARAM, by = c("HELCOM_subbasin", "CRUIS","STNNO", "Date", "STATN", 
                                          "Latitude", "Longitude")) %>%
  dplyr::filter(n == 2) %>%
  select(-n)

rm(MSTS_data_onePARAM)

####################################################################################.
####################################################################################.
#####.
### 2.4. Layers to whole column #####

sampled_depth <- MSTS_data %>%
  group_by(CNTRY, CRUIS,STNNO, Latitude, Longitude,Date,STATN,HELCOM_subbasin, ICES_area, PARAM, MUNIT_msts) %>%
  summarise(MXdepth_sampled = max(MXDEP))

MSTS_data <- MSTS_data %>%
  
  ## 
  left_join(., sampled_depth, 
            by = c("HELCOM_subbasin","CNTRY", "CRUIS","STNNO", "Latitude", "Longitude", "Date", "STATN", 
                   "ICES_area", "PARAM", "MUNIT_msts")) %>%
  group_by(CNTRY, CRUIS,STNNO, Latitude, Longitude,MYEAR, Date,STATN,HELCOM_subbasin,
           ICES_area, MXdepth_sampled, PARAM, MUNIT_msts) %>%
  summarise(MNDEP = MNDEP, MXDEP = MXDEP, value = Value_msts_sum) %>%
  mutate(lay_value = (MXDEP-MNDEP)*value) %>%
  ungroup() %>% 
  group_by(CNTRY, CRUIS,STNNO, Latitude, Longitude,MYEAR, Date,STATN,HELCOM_subbasin,
           ICES_area, MXdepth_sampled, PARAM, MUNIT_msts) %>%
  summarise(col_val = sum(lay_value)/MXdepth_sampled) %>% unique() %>%
  mutate(Layer_considered = paste0("0--", MXdepth_sampled)) %>% ungroup() %>%
  select(-MXdepth_sampled)

## remove variable 'sampled_depth'; information merged with MSTS_data
rm(sampled_depth)
