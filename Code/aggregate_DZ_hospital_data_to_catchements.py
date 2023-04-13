
import pandas as pd
import matplotlib.pyplot as plt
from datetime import datetime
from datetime import timedelta
import numpy as np
import pickle as pk


data_date='dummy'
##################### Sites #################
# this file contains the datazones of each catchment
DZ_df=pd.read_csv('../Data/dz_ssa_summary_20210216.csv')

# this file contains the alternative names of each site
# for the dumy data read only part of the file
look_up_df=pd.read_csv('../Data/catchmentsitelookup.csv',dtype=str,nrows=10)

# make a list of the two rows
# catchmentname is how it is listed in the SEPA file
catchmentname=look_up_df['catchmentname'].tolist()
# WWRNA_name is the name in the DZ-to-catchment file
WWRNA_name=look_up_df['WWRNA_simplename'].tolist()
# make a conversion dictionary DZ name to SEPA name
other_name_of={}
while catchmentname:
    other_name_of[WWRNA_name.pop()]=catchmentname.pop()

# list the names
WWRNA_names=sorted(list(set(look_up_df['WWRNA_simplename'].tolist())))

# include this to get the estmated numbers of each type
admissions=pk.load(open('../pickles/DZ_level_admissions/admissions_'+data_date+'.p','rb'))


#measure days from this day
time_zero=datetime.strptime('2021-05-01', '%Y-%m-%d')
# how many days do we have? use a random DZ
#random DZ
days=len(list(admissions.values())[0])

# get dates
date_list=[datetime.date((time_zero+timedelta(day))) for day in range(days)]

# recored population of each site (optional)
population={}

for site in WWRNA_names:    
    # select rows that correspond to the site
    site_df=DZ_df[DZ_df['Sampled Sewer Area']==other_name_of[site]]
    
    # the sum the populations of all the DZs
    population[site]=sum(site_df['Pop2011'].tolist())
    
    # initialise admission counts at 0
    adms=[0 for j in range(days)]
    
    # loop over all the DZs in the site
    for i,row in site_df.iterrows():
        
        DZ=row['DataZone2011']
        if DZ in admissions:
        
            prop=row['prop_dz']
            #print(DZ)
            for j in range(days):
                # add the appropriate proportion to the total
                adms[j]=adms[j]+admissions[DZ][j]*prop
            
            
    output={'Date':date_list,
            'admissions':adms
            }

    # store the output under the SEPA name
    pd.DataFrame(output).to_csv('../Site_level_admissions/'+site+'_admissions.csv',index=False)
    #print('|'+site+'|',len(site))
    #print('|'+site.strip(' ')+'|',len(site))
    #print(pd.DataFrame(output).head())
    print(site,len(date_list))
    

#pk.dump(population,open('../pickles/site_populations.p','wb'))
output={'Site':[],'Population':[]}

#print()
sorted_pops=sorted([(site,population[site]) for site in population], key=lambda item: item[1],reverse=True)
for site,population in sorted_pops:
    print(site,population)
    output['Site'].append(site)
    output['Population'].append(population)

pd.DataFrame(output).to_csv('../Data/population_count.csv',index=False)


