# ZooplanktonMeanSizeTotalStock
### Brief intro to the indicator

The HELCOM core indicator Zooplankton mean size and total stock (MSTS) is used to assess the status of pelagic habitat in the HELCOM Holistic Assessment of the Baltic Sea. The indicator report with information on indicator concept and results can be found on the HELCOM website, [here](https://helcom.fi/baltic-sea-trends/indicators/).

___

### Indicator data

The indicator use data collected as part of the [HELCOM COMBINE monitoring program](https://helcom.fi/action-areas/monitoring-and-assessment/monitoring-guidelines/combine-manual/). Data are reported to the COMBINE database maintained by [ICES](https://data.ices.dk/). The indicator calculation script here use data extracts from the COMBINE database, as well as it allows to use aggregated data (yearly means) as input. 

___

### Input files
##### **1. 'MSTS_addPar.txt'**

| Column name | Variable type | Description |
|  :-----  |  :---  |  :-------- |
| HELCOM_subbasin | Character | name of a subbasin (assessment unit) of the Baltic Sea [as defined in HELCOM](https://gis.ices.dk/geonetwork/srv/api/records/225df9db-bfdf-4388-8ccb-fa4b99053a36) |
| RefFISH | Character | reference period defined based on expert opinion and Weight-At-Age or other herring/sprat feeding condition indices. year_from:year_to |
| RefCHL | Character | reference period defined based on expert opinion and chlorophyll-a, Secchi depth or other eutrophication proxy. year_from:year_to |
| Hist_data_source | Character | data source to be used (note: 'input' or 'COMBINE') for calculation of the threshold values (and supplementary parameters: lambda, ref-mean, ref-sd). input - if yearly means are to be used for reference period; COMBINE - if data extract from ICES is to be used for reference period |
| defGES_TZB[mg_m-3] | Numeric | defined threshold value for Total Stock (total zooplankton biomass) parameter |
| defGES_MS[ug_ind-1] | Numeric | defined threshold value for MEan Size parameter |
| Per_start[month] | Numeric | reference period defined based on expert opinion and chlorophyll-a, Secchi depth or other eutrophication proxy | the starting month of the period considered as relevant for specific subbasin (period used for threshold value calculation and to be used for the assessment) |
| Per_end[month] | Numeric | reference period defined based on expert opinion and chlorophyll-a, Secchi depth or other eutrophication proxy | the ending month of the period considered as relevant for specific subbasin (period used for threshold value calculation and to be used for the assessment) |
| Depth_minLimit[m] | Numeric | data covering layer to the defined depth or deeper are included in the calculation |
| Depth_maxLimit[m] | Numeric | data covering layer not deeper than defined are included in the calculations |

##### **2. 'MSTS_hist_data.txt'**

Provides information on assessment units idetified as 'input' in Hist_data_source column in 'MSTS_addPar.txt'. No missing data are allowed in timeseries. If NA's are present an estimate must be provided before running the script!

| Column name | Variable type | Description |
|  :-----  |  :---  |  :-------- |
| year | Numeric | year |
| HELCOM_subbasin | Character | name of a subbasin (assessment unit) of the Baltic Sea [as defined in HELCOM](https://gis.ices.dk/geonetwork/srv/api/records/225df9db-bfdf-4388-8ccb-fa4b99053a36) |
| ICES_area | Character | name of a subbasin of the Baltic Sea [as defined in ICES](https://gis.ices.dk/geonetwork/srv/eng/catalog.search#/metadata/c784a0a3-752f-4b50-b02f-f225f6c815eb) |
| RefFISH | Logical | reference period defined based on expert opinion and Weight-At-Age or other herring/sprat feeding condition indices (1 - used as a reference year; 0 - not used as a reference year) |
| RefCHL | Logical | reference period defined based on expert opinion and chlorophyll-a, Secchi depth or other eutrophication proxy (1 - used as a reference year; 0 - not used as a reference year) |
| TZB | Numeric | Total zooplankton biomass (Total Stock proxy), annual mean for period relevant in MSTS-assessment, see 'MSTS_addPar.txt' Per_start[month], Per_end[month] |
| MS | Numeric | Mean Size of zooplankton, annual mean for period relevant in MSTS-asessment, see 'MSTS_addPar.txt' Per_start[month], Per_end[month] |
| Notes | Character | Notes - if relevant |
| Data from | Character | Raw data holder - if relevant |
| Reference | Character | Reference to data aggregator/s - if relevant |

##### **3. 'MSTS_ZPspecies.txt'**

Redundant if only MSTS-relevant species are included in data, i.e., without predatory cladocerans and meroplankton). Includes all variations of zp taxa reported to ICES.

| Column name | Variable type | Description |
|  :-----  |  :---  |  :-------- |
| SPECI | Character | Taxon code as reported to ICES |
| SPECI_name | Character | Taxon name as reported to ICES |
| AphiaID | Numeric | AphiaID as defined in [WoRMS](https://www.marinespecies.org/) |
| WoRMS_name | Character | Taxon name as defined in [WoRMS](https://www.marinespecies.org/) |
| AphiaID_accepted | Numeric | AphiaID for the accepted taxon name as defined in [WoRMS](https://www.marinespecies.org/) |
| WoRMS_accepted_name | Numeric | Accepted taxon name as defined in [WoRMS](https://www.marinespecies.org/) |
| MSTS_use | Logical | describing whether taxon is used for MSTS-assessment, only herbivorous holoplankton are included (1 - yes; 0- no) |

##### **4. 'MSTS_statn.txt'**

Information relevant for data aggregation (if necesaray). It defines stations (STATN) to be included in the MSTS calculations.

| Column name | Variable type | Description |
|  :-----  |  :---  |  :-------- |
| HELCOM_subbasin | Character | name of a subbasin (assessment unit) of the Baltic Sea [as defined in HELCOM](https://gis.ices.dk/geonetwork/srv/api/records/225df9db-bfdf-4388-8ccb-fa4b99053a36) |
| STATN | Character | station name as reported to ICES DOME |
| MSTS_calc | Logical | describing whether Station is used for calculations of MSTS threshold values (1 - yes; 0- no) |
| Layer_considered | Character | describing the type of water layer considered in MSTS for the specific Station data (C - whole water column, from top to bottom; TXX - top XX meters of the water column |
| BotDepth_considered[m] | Numeric | bottom depth of the water layer considered for MSTS in meters |

##### **5. ICES example file 'ICES_HELCOM_zp.txt'**
An example data from https://data.ices.dk/view-map describing zooplankton biological community. Example data are limited to represent four HELCOM subbasins: Bothnian Bay, Bothnian Sea, Gulf of Riga, Western Gotland Basin.

___

### Output files

* Tabular outputs

  1. **Aggregated data**: abundance and wet weight per sample covering the defined depth layer
  
  2. **Fin Assessment Data**: 
   
    * MSTS assessment data 'MSTS_assessment_df_[date].txt'. Holds yearly data per HELCOM subbasin (assessment unit) on untransformed and transformed TZB and MS, z-scores, lowerCusum, Alarm column (1 for sub-GES; 0 for GES)
   
    * MSTS assessment confidence values 'MSTS_confidence_[date].txt'. n_ObPerYear - number of datapoins (samples) per year included in the assessment (temporal confidence); n_STATNobserved - number of different stations included in the assessment (spatial confidence).
   
    * GES values calculated/used in the assessment 'MSTS_GES_[date].txt'.
  
* **Graphical outputs**
  
  1. Trend plot with Kendall's statistics for TZB and MS
  
  2. Trend plot for lower CuSum
  
  3. MSTS plot
  
___





