# Estimating-epidemic-growth-rate-from-wastwater-data
Codes associated with the paper "The impact of signal variability on epidemic growth rate estimation from wastewater surveillance data"

The data used for this project is protected. This repo includes a script to create dummy data of the same format.

## Contents
### create_dummy_data.py
Creates fake hospital data and fake wastewater data in the same format as the data made available to us from PHS and SEPA

### aggregate_raw_hospital_data_to_DZs.py
Aggregates hospital data to the datazone (DZ) level producing a dictionary of lists  (a time series for each DZ)  

### aggregate_DZ_hospital_data_to_catchements.py
Aggregates admissions from datazone level to the site catchments for each wastwateter treatment site. Uses files kept in the data folder to make the mapping. Produces a number of csv files that populate the "Site level admissions" folder.

### WWlibrary.py
A set of functions which perform the main computations: likelihood calculation, hill-climb optimisation of parameters, choosing whether to add a new inflection point

### WW_all_thresholds.py
Loops over some of the sites and a range of threshold values (threshold is the main hyperparametr in the analysis). Hard coded variable data_source determines whether to analyse the hospital data or the wastewater data. Ouputs are pickles stored in the relevant subdirectory of the pickles folder

### WW_all_sites.py
Loops over all the sites. Hard coded variable data_source determines whether to analyse the hospital data or the wastewater data. Ouputs are pickles stored in the relevant subdirectory of the pickles folder

### figure_population_vs_likelihood.py
Produces figure S1 from the paper.

### figure_threshold_vs_changes.py
Produces figure 3 from the paper.

### figure_threshold_vs_waves.py
Produces figure 4 from the paper.

### figure_growth_rate_lag.py
Produces figure 5 from the paper.

### residuals_vs_flow_top10.py
Produces figure S2 from the paper and outputs some results in the terminal that are mentioned in the last paragraph of the results section.



