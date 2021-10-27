
if("EnvStats" %in% rownames(installed.packages()) == FALSE) {install.packages("EnvStats")} ## for normalization boxcox approach
if("fBasics" %in% rownames(installed.packages()) == FALSE) {install.packages("fBasics")}  ## dagoTest()

################################################################################################.
# source file for 'MAIN_MSTS_calc.R'
# performs calculation on aggregated data for MSTS indicator (data without NAs)

### Input data ####
################################################################################################.

h_df <- MSTS_data_input %>%
  select(-Notes, -`Data from`, -Reference) %>%
  group_by(HELCOM_subbasin, ICES_area) %>%        #group by subbasins
  group_split()                                   #separate by groups (subbasins) into lists

# name list object as HELCOM_subbasin
namesList <- rep("", length(h_df))
for (subb in seq_along(h_df)) {
  namesList[[subb]] <- unique(h_df[[subb]]$HELCOM_subbasin) 
}


names(h_df) <- namesList
rm(namesList)

## Normality check ####
################################################################################################.
# create empty list for storing Normality check results
norm_tzb <- vector("list", length(h_df))
names(norm_tzb) <- names(h_df)
norm_ms <- vector("list", length(h_df))
names(norm_ms) <- names(h_df)

# create empty list for storing box-cox lambda for every subbasin
l_tzb <- vector("list", length(h_df))
names(l_tzb) <- names(h_df)
l_ms <- vector("list", length(h_df))
names(l_ms) <- names(h_df)

for (subb in seq_along(h_df)){
  ## check normality for TZB
  norm_tzb[[subb]]$dagoTest <- fBasics::dagoTest(h_df[[subb]]$TZB)@test   # D'Agostino & Pearson omnibus normality test
  norm_tzb[[subb]]$shapiroTest <- shapiro.test(h_df[[subb]]$TZB)          # Shapiro-Wilk normality test
  norm_tzb[[subb]]$ksTest <- (ks.test(h_df[[subb]]$TZB, "pnorm",          # Kolmogorov-Smirnov normality test
                                      mean = mean(h_df[[subb]]$TZB), 
                                      sd = sd(h_df[[subb]]$TZB)))
  ## check normality for MS
  norm_ms[[subb]]$dagoTest <- fBasics::dagoTest(h_df[[subb]]$MS)@test     # D'Agostino & Pearson omnibus normality test
  norm_ms[[subb]]$shapiroTest <- shapiro.test(h_df[[subb]]$MS)            # Shapiro-Wilk normality test
  norm_ms[[subb]]$ksTest <- (ks.test(h_df[[subb]]$MS, "pnorm",            # Kolmogorov-Smirnov normality test
                                     mean = mean(h_df[[subb]]$MS), 
                                     sd = sd(h_df[[subb]]$MS)))
}  

# create empty list for logical variable identifying if box-cox transformation is needed
transf_tzb <- rep("", length(h_df))
names(transf_tzb) <- names(h_df)
transf_ms <- rep("", length(h_df))
names(transf_ms) <- names(h_df)

for (subb in seq_along(h_df)){
  if((norm_tzb[[subb]]$dagoTest$p.value[1] > 0.1 &
      norm_tzb[[subb]]$shapiroTest$p.value > 0.1 &
      norm_tzb[[subb]]$ksTest$p.value > 0.1)) {
    transf_tzb[subb] <- "Passed" 
  }else{transf_tzb[subb] <- "Failed"}
  
  if((norm_ms[[subb]]$dagoTest$p.value[1] > 0.1 &
      norm_ms[[subb]]$shapiroTest$p.value > 0.1 &
      norm_ms[[subb]]$ksTest$p.value > 0.1)) {
    transf_ms[subb] <- "Passed" 
  }else{transf_ms[subb] <- "Failed"}
} 

## TRUE for Box-Cox calculations
# lambda calculation NAs ignore
for (subb in seq_along(h_df)){
  if(transf_tzb[subb] == "Failed"){
    l_tzb[[subb]] = round(EnvStats::boxcox(h_df[[subb]]$TZB, optimize = TRUE, 
                                           objective.name = "Log-Likelihood")$lambda, 3)
  }
}
for (subb in seq_along(h_df)){
  if(transf_ms[subb] == "Failed"){
    l_ms[[subb]] = round(EnvStats::boxcox(h_df[[subb]]$MS, optimize = TRUE, 
                                          objective.name = "Log-Likelihood")$lambda, 3)
  }
}

### Data transformation ####
################################################################################################.
#### Define Box-Cox transformation

BCTransform <- function(y, lambda=0) {
  if (lambda == 0L) { log(y) }
  else { (y^lambda - 1) / lambda }
}

#### Define Box-Cox inverse transformation
BCTransformInverse <- function(yt, lambda=0) {
  if (lambda == 0L) { exp(yt) }
  else { exp(log(1 + lambda * yt)/lambda) }
}

## transform data accordingly to the reported lambdas: TZB
for (subb in seq_along(h_df)){
  if(transf_tzb[subb] == "Passed"){
    h_df[[subb]]$TZB <- h_df[[subb]]$TZB
  }else{
    h_df[[subb]]$TZB <- BCTransform(h_df[[subb]]$TZB, lambda = l_tzb[[subb]])
  }
}

## transform data accordingly to the reported lambdas: MS
for (subb in seq_along(h_df)){
  if(transf_ms[subb] == "Passed"){
    h_df[[subb]]$MS <- h_df[[subb]]$MS
  }else{
    h_df[[subb]]$MS <- BCTransform(h_df[[subb]]$MS, lambda = l_ms[[subb]])
  }
}

MSTS_normData <- data.frame(Reduce(rbind, h_df))
rm(h_df)

### Ref mean and SD ####
####################################################################################.

#### CHL Ref period ####
refP_chl_GES <- MSTS_normData[MSTS_normData$RefCHL == 1,] %>%
  select(-RefFISH) %>%
  group_by(HELCOM_subbasin, ICES_area) %>%
  
  #calculate mean and sd for the ref period CHL  
  summarize(chl_tzb_mean = mean(TZB), 
            chl_ms_mean = mean(MS), 
            chl_tzb_sd = sd(TZB), 
            chl_ms_sd = sd(MS),
            n_yearsCHL = sum(RefCHL)) %>%
  
  #calculate Lower 99% confidence level for the ref period CHL
  mutate(chl_tzb_L99CL = chl_tzb_mean - (qnorm(0.995)*chl_tzb_sd/sqrt(n_yearsCHL)),
         chl_ms_L99CL = chl_ms_mean - (qnorm(0.995)*chl_ms_sd/sqrt(n_yearsCHL)),
         l_tzb = 0,
         l_ms = 0,
         chl_tzb_GES = -1,
         chl_ms_GES = -1)
  
  # inverse box-cox transformation for the calculated Lower 99% confidence level for the ref period CHL 
  #  to obtain absolute GES values
for (subb in seq_along(l_tzb)){
  if (is.null(l_tzb[[subb]])){
    #RefCHL TZB / no retransformation needed
    refP_chl_GES[refP_chl_GES$HELCOM_subbasin == names(l_tzb[subb]),]$chl_tzb_GES <- 
      refP_chl_GES[refP_chl_GES$HELCOM_subbasin == names(l_tzb[subb]),]$chl_tzb_GES
    refP_chl_GES[refP_chl_GES$HELCOM_subbasin == names(l_tzb[subb]),]$l_tzb <- NA
  }else{
    #RefCHL TZB / retransformation
    refP_chl_GES[refP_chl_GES$HELCOM_subbasin == names(l_tzb[subb]),]$chl_tzb_GES <-
        BCTransformInverse(refP_chl_GES[refP_chl_GES$HELCOM_subbasin == names(l_tzb[subb]),]$chl_tzb_L99CL,
                           lambda = l_tzb[[subb]])
    refP_chl_GES[refP_chl_GES$HELCOM_subbasin == names(l_tzb[subb]),]$l_tzb = l_tzb[[subb]]
    }
}
for (subb in seq_along(l_ms)){
  if (is.null(l_ms[[subb]])){
    #RefCHL MS / no retransformation needed
    refP_chl_GES[refP_chl_GES$HELCOM_subbasin == names(l_ms[subb]),]$chl_ms_GES <- 
      refP_chl_GES[refP_chl_GES$HELCOM_subbasin == names(l_ms[subb]),]$chl_ms_L99CL
    refP_chl_GES[refP_chl_GES$HELCOM_subbasin == names(l_ms[subb]),]$l_ms <- NA
  }else{
    #RefCHL MS / retransformation
    refP_chl_GES[refP_chl_GES$HELCOM_subbasin == names(l_ms[subb]),]$chl_ms_GES <-
      BCTransformInverse(refP_chl_GES[refP_chl_GES$HELCOM_subbasin == names(l_ms[subb]),]$chl_ms_L99CL,
                         lambda = l_ms[[subb]])
    refP_chl_GES[refP_chl_GES$HELCOM_subbasin == names(l_ms[subb]),]$l_ms = l_ms[[subb]]
  }
}

refP_chl_GES <- refP_chl_GES %>%
  gather(abbrev, value,chl_tzb_mean:chl_ms_GES)

#### FISH Ref period ####

refP_fish_GES <- MSTS_normData[MSTS_normData$RefFISH == 1,] %>%
  select(-RefCHL) %>%
  group_by(HELCOM_subbasin, ICES_area) %>%
  
  #calculate mean and sd for the ref period FISH  
  summarize(fish_tzb_mean = mean(TZB), 
            fish_ms_mean = mean(MS), 
            fish_tzb_sd = sd(TZB), 
            fish_ms_sd = sd(MS),
            n_yearsFISH = sum(RefFISH)) %>%
  
  #calculate Lower 99% confidence level for the ref period CHL
  mutate(fish_tzb_L99CL = fish_tzb_mean - (qnorm(0.995)*fish_tzb_sd/sqrt(n_yearsFISH)),
         fish_ms_L99CL = fish_ms_mean - (qnorm(0.995)*fish_ms_sd/sqrt(n_yearsFISH)),
         l_tzb = 0,
         l_ms = 0,
         fish_tzb_GES = -1,
         fish_ms_GES = -1)

# inverse box-cox transformation for the calculated Lower 99% confidence level for the ref periodFISH 
#  to obtain absolute GES values
for (subb in seq_along(l_tzb)){
  if (is.null(l_tzb[[subb]])){
    #RefFISH TZB / no retransformation needed
    refP_fish_GES[refP_fish_GES$HELCOM_subbasin == names(l_tzb[subb]),]$fish_tzb_GES <- 
      refP_fish_GES[refP_fish_GES$HELCOM_subbasin == names(l_tzb[subb]),]$fish_tzb_GES
    refP_fish_GES[refP_fish_GES$HELCOM_subbasin == names(l_tzb[subb]),]$l_tzb <- NA
  }else{
    #RefFISH TZB / retransformation
    refP_fish_GES[refP_fish_GES$HELCOM_subbasin == names(l_tzb[subb]),]$fish_tzb_GES <-
      BCTransformInverse(refP_fish_GES[refP_fish_GES$HELCOM_subbasin == names(l_tzb[subb]),]$fish_tzb_L99CL,
                         lambda = l_tzb[[subb]])
    refP_fish_GES[refP_fish_GES$HELCOM_subbasin == names(l_tzb[subb]),]$l_tzb = l_tzb[[subb]]
  }
}
for (subb in seq_along(l_ms)){
  if (is.null(l_ms[[subb]])){
    #RefFISH MS / no retransformation needed
    refP_fish_GES[refP_fish_GES$HELCOM_subbasin == names(l_ms[subb]),]$fish_ms_GES <- 
      refP_fish_GES[refP_fish_GES$HELCOM_subbasin == names(l_ms[subb]),]$fish_ms_L99CL
    refP_fish_GES[refP_fish_GES$HELCOM_subbasin == names(l_ms[subb]),]$l_ms <- NA
  }else{
    #RefFISH MS / retransformation
    refP_fish_GES[refP_fish_GES$HELCOM_subbasin == names(l_ms[subb]),]$fish_ms_GES <-
      BCTransformInverse(refP_fish_GES[refP_fish_GES$HELCOM_subbasin == names(l_ms[subb]),]$fish_ms_L99CL,
                         lambda = l_ms[[subb]])
    refP_fish_GES[refP_fish_GES$HELCOM_subbasin == names(l_ms[subb]),]$l_ms = l_ms[[subb]]
  }
}

refP_fish_GES <- refP_fish_GES %>%
  gather(abbrev, value,fish_tzb_mean:fish_ms_GES)


MSTS_ref_outputs <- rbind(refP_chl_GES, refP_fish_GES) %>%
    separate("abbrev", c("RefP", "variable", "parameter"), sep = "_" ) %>%
    unique()



#### GES ######

MSTS_GES_tzb <- MSTS_ref_outputs %>%
  group_by(HELCOM_subbasin, ICES_area) %>%
  dplyr::filter(variable == "tzb", parameter == "GES") %>%
  summarize(GES_tzb = max(value))

MSTS_GES <- MSTS_ref_outputs %>%
  group_by(HELCOM_subbasin, ICES_area) %>%
  dplyr::filter(variable == "ms", parameter == "GES") %>%
  summarize(GES_ms = max(value)) %>%
  full_join(., MSTS_GES_tzb, by = c("HELCOM_subbasin", "ICES_area"))

rm(MSTS_GES_tzb)

############################################################################################.

