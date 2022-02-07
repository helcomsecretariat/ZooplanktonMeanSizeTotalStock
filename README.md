# ZooplanktonMeanSizeTotalStock
### Brief intro to the indicator

The HELCOM core indicator Zooplankton mean size and total stock (MSTS) is used to assess the status of pelagic habitat in the HELCOM Holistic Assessment of the Baltic Sea. The indicator report with information on indicator concept and results can be found on the HELCOM website, [here](https://helcom.fi/baltic-sea-trends/indicators/).

### Indicator data

The indicator use data collected as part of the [HELCOM COMBINE monitoring program](https://helcom.fi/action-areas/monitoring-and-assessment/monitoring-guidelines/combine-manual/). Data are reported to the COMBINE database maintained by [ICES](https://data.ices.dk/). The indicator calculation script here use data extracts from the COMBINE database, as well as it allows to use aggregated data (yearly means) as input. 

### Inputs and outputs {.tabset}
#### Input files

* Field descriptions for **'MSTS_addPar.txt'**
  a. HELCOM_subbasin: /Character variable/ name of a subbasin (assessment unit) of the Baltic Sea as defined in HELCOM (https://gis.ices.dk/geonetwork/srv/api/records/225df9db-bfdf-4388-8ccb-fa4b99053a36)
  
  b. RefFISH: /Character variable/ reference period defined based on expert opinion and Weight-At-Age or other herring/sprat feeding condition indices
  
  c. RefCHL: /Character variable/ reference period defined based on expert opinion and chlorophyll-a, Secchi depth or other eutrophication proxy
  
  d. Hist_data_source: /Character variable/ data source to be used (note: 'input' or 'COMBINE') for calculation of the threshold values (and supplementary parameters: lambda, ref-mean, ref-sd). input - if yearly means are to be used for reference period; COMBINE - if data extract from ICES is to be used for reference period 
  
  e. defGES_TZB[mg_m-3]: /Numeric variable/ defined threshold value for Total Stock (total zooplankton biomass) parameter
  
  f. defGES_MS[$\mu$g_ind-1]: /Numeric variable/ defined threshold value for Mean Size parameter

  g. Per_start[month]: /Numeric variable/ the starting month of the period considered as relevant for specific subbasin (period used for threshold value calculation and to be used for the assessment)
  
  h. Per_end[month]: /Numeric variable/ the ending month of the period considered as relevant for specific subbasin (period used for threshold value calculation and to be used for the assessment)
  
  i. Depth_minLimit[m]: /Numeric variable/ data covering layer to the defined depth or deeper are included in the calculation
  
  j. Depth_maxLimit[m]: /Numeric variable/ data covering layer not deeper than defined are included in the calculations

* Field descriptions for **'MSTS_hist_data.txt'**. No missing data are allowed in timeseries. If NA's are present an estimate must be provided before running the script!

  a. year: /Numeric variable/ year of the data
  
  b. HELCOM_subbasin: /Character variable/ name of a subbasin (assessment unit) of the Baltic Sea as defined in HELCOM (https://gis.ices.dk/geonetwork/srv/api/records/225df9db-bfdf-4388-8ccb-fa4b99053a36)
  
  c. ICES_area: /Character variable/ name of a subbasin of the Baltic Sea as defined in ICES (https://gis.ices.dk/geonetwork/srv/eng/catalog.search#/metadata/c784a0a3-752f-4b50-b02f-f225f6c815eb)
  
  d. RefFISH: /Logical variable/ reference period defined based on expert opinion and Weight-At-Age or other herring/sprat feeding condition indices (1 - used as a reference year; 0 - not used as a reference year)
  
  e. RefCHL: /Logical variable/ reference period defined based on expert opinion and chlorophyll-a, Secchi depth or other eutrophication proxy (1 - used as a reference year; 0 - not used as a reference year)
  
  f. TZB: /Numeric variable/ Total zooplankton biomass (Total Stock proxy), annual mean for period relevant in MSTS-assessment, see 'MSTS_addPar.txt' Per_start[month], Per_end[month]
  
  g. MS: /Numeric variable/ Mean Size of zooplankton, annual mean for period relevant in MSTS-asessment, see 'MSTS_addPar.txt' Per_start[month], Per_end[month]
  
  h. Notes: /Character variable/ Notes - if relevant
  
  i. Data from: /Character variable/ Raw data holder - if relevant
  
  j. Reference: /Character variable/ Reference to data aggregator/s - if relevant

* Field descriptions for **'MSTS_ZPspecies.txt'** (redundant if only MSTS-relevant species are included in data, i.e., without predatory cladocerans and meroplankton). Includes all variations of zp taxa reported to ICES.

  a. SPECI: /Character variable/

  b. SPECI_name: /Character variable/

  c. AphiaID: /Numeric variable/

  d. WoRMS_name: /Character variable/
  
  e. AphiaID_accepted: /Character variable/
  
  f. WoRMS_accepted_name: /Character variable/
  
  g. MSTS_use: /Logical variable/ describing whether Species is used for MSTS-assessment, only herbivorous holoplankton are included (1 - yes; 0- no)

* Field descriptions for **'MSTS_statn.txt'**. Information relevant for data aggregation if reference data provided in unaggregated way. It defines stations (STATN) to be included in the reference threshold value setting.

  a. HELCOM_subbasin: /Character variable/ name of a subbasin of the Baltic Sea as defined in HELCOM (https://gis.ices.dk/geonetwork/srv/api/records/225df9db-bfdf-4388-8ccb-fa4b99053a36)
  
  b. STATN: /Character variable/ station name as reported to ICES DOME (see ICES Station Vocabulary)
MSTS_calc: /Logical variable/ describing whether Station is used for calculations of MSTS threshold values (1 - yes; 0- no)

  c. Layer_considered: /Character variable/ describing the type of water layer considered in MSTS for the specific Station data (C - whole water column, from top to bottom; TXX - top XX meters of the water column)
  
  d. BotDepth_considered[m]: /Numeric variable/ bottom depth of the water layer considered for MSTS in meters

* **ICES example file 'ICES_HELCOM_zp.txt'**. 
An example data from https://data.ices.dk/view-map describing zooplankton biological community. Example data are limited to represent four HELCOM subbasins: Bothnian Bay, Bothnian Sea, Gulf of Riga, Western Gotland Basin.

#### Output files

* Tabular outputs

  1. **Aggregated data**: abundance and wet weight per sample covering the defined depth layer
  
  2. **Fin Assessment Data**: 
   
   * MSTS assessment data 'MSTS_assessment_df_[date].txt'. Holds yearly data per HELCOM subbasin (assessment unit) on untransformed and transformed TZB and MS, z-scores, lowerCusum, Alarm column (1 for sub-GES; 0 for GES)
   
   * MSTS assessment confidence values 'MSTS_confidence_[date].txt'. n_ObPerYear - number of datapoins (samples) per year included in the assessment (temporal confidence); n_STATNobserved - number of different stations included in the assessment (spatial confidence).
   
   * GES values calculated/used in the assessment 'MSTS_GES_[date].txt'.
  
  3. **Graphical outputs**
  
  * Trend plot with Kendall's statistics for TZB and MS
  
  * Trend plot for lower CuSum
  
  * MSTS plot
  





