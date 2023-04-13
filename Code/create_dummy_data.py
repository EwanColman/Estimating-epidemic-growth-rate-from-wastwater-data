# Creates a dunmmy hospital data file
# uses ONS and site look up files

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from datetime import datetime
from datetime import timedelta
from random import random
from random import choice


# approximate population of made up area covered by catchments
population=4000000

# time frame of data
duration=100

# to generate dummy data we need a list some prevalence from the ONS and a list of datazones

################ Catchment information ###################
# this file contains the datazones of each catchment
DZ_df=pd.read_csv('../Data/dz_ssa_summary_20210216.csv')

datazone_list=DZ_df['DataZone2011'].tolist()
print(DZ_df.head())

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

# measure days from this day
time_zero=datetime.strptime('2021-01-01', '%Y-%m-%d')
################### ONS data ##################
df=pd.read_excel('../Data/20230113covid19infectionsurveydatasetsscotland1.xlsx',sheet_name='1a',usecols='A:B',names=['Time','Rate'],skiprows=6,nrows=115).replace('-',0)


prevalence=[]
for i,row in df.iterrows():
    t=row['Time']
    start=(datetime.strptime(t[:t.find('to ')-1].strip(),'%d %B %Y')-time_zero).days
    end=(datetime.strptime(t[t.find('to ')+3:].strip(),'%d %B %Y')-time_zero).days
    # get middle of interval

    # only adds times after time_zero
    prevalence=prevalence+[round(population*row['Rate']/100) for i in range(max(0,start),max(0,end))]

###################### Generate synthetic data ###################################

# simulate the bernoulli trial hospitalisation
# alpha is the probability of a non-zero
delta_hfr=0.007
omicron_hfr=0.001

omicron_time=(datetime.strptime('2021-11-01', '%Y-%m-%d')-time_zero).days

synthetic_admissions=[]
for t in range(duration):
    if t<omicron_time:
        alpha=delta_hfr
    else:
        alpha=omicron_hfr
        

    # add it up one shedder at a time
    admissions=0
    for i in range(prevalence[t]):
        if random()<alpha:
            admissions=admissions+1

    synthetic_admissions.append(admissions)
# alternative: use the normal distribution directly

# remove 0s
sample_value=[]
sample_time=[]
for t in range(duration):
    if synthetic_admissions[t]>0:
    
        sample_value.append(synthetic_admissions[t])
        sample_time.append(t)

################## Plot #####################
plt.scatter(sample_time,sample_value)
plt.xlabel('Time (days)')
plt.ylabel('Hospitlisation')
#############################################

# assign admissions to datazones

output={'id':[],
        'admission_date':[],
        'datazone_2011':[],
        'main_condition':[],
        'other_condition_1':[],
        'other_condition_2':[],
        'other_condition_3':[],
        'other_condition_4':[],
        'other_condition_5':[]
        }

datazone_admissions={}
for DZ in datazone_list:
    datazone_admissions[DZ]=[0 for t in range(duration)]

n=0

for t in range(duration):
   
    date=(time_zero+timedelta(t)).strftime('%Y%m%d')
    
    for i in range(synthetic_admissions[t]):
        # pick a random DZ
        DZ=choice(datazone_list)
        
        #for each admission update the lists
        output['id'].append('id'+str(n))
        n=n+1
        output['admission_date'].append(date)
        output['datazone_2011'].append(DZ)
        
        

        if random()<0.1:
            output['main_condition'].append('U071')
            for i in range(1,6): 

                output['other_condition_'+str(i)].append('noncovid')


        else:
            output['main_condition'].append('noncovid')
   
            # include 6 and 7 to get some fully noncovid entries
            m=choice([1,2,3,4,5,6,7])
            for i in range(1,6): 
                if i==m:
                    output['other_condition_'+str(i)].append('U071')
                else:
                    output['other_condition_'+str(i)].append('noncovid')

pd.DataFrame(output).to_csv('../Data/smr01_dummy.csv',index=False)


######### dummy wastewater data #############

output={'Date':[],'Site':[],'N1 Reported Value':[],'Flow (l/day)':[],'Ammonia (mg/l)':[]}

# get list of sites (use only the first 10)
WWRNA_name=look_up_df['WWRNA_simplename'].tolist()[0:10]

n=0
for site in WWRNA_name:
    print(n,site)
    n=n+1
    # samples every 2 or 3 days
    sample_time=[0]
    for i in range(int(duration/2)):
        if random()<0.5:
            sample_time.append(sample_time[-1]+2)
        else:
            sample_time.append(sample_time[-1]+3)
    
    sample_time=[t for t in sample_time if t<duration]
    # format like SEPA data
    sample_date=[(time_zero+timedelta(t)).strftime('%d/%m/%Y') for t in sample_time]
    
    
    # choose how much dispersion there will be
    D=0.05 # based on what we estimated
    
    max_RNA_count=1
    # generate the number of people actively shedding
    # divide by 10 to rescale 
    infections=[int(prevalence[t]/10) for t in sample_time]
    # params of the rna contrbution distribution
    
    # simulate the sampling contribution process
    # zero-inflated distribution
    # alpha is the probability of a non-zero
    delta_shedding_proportion=0.003
    omicron_shedding_proportion=0.001

    # simulate the sampling    
    sample_value=[]
    for i in range(len(infections)):
        if sample_time[i]<omicron_time:
            alpha=delta_shedding_proportion
        else:
            alpha=omicron_shedding_proportion
        
        mu=0.006
        sigma=0.03
        # add it up one shedder at a time
        gene_copies=0
        for i in range(infections[i]):
            if random()<alpha:
                #gene_copies=gene_copies+np.random.gamma((mu/sigma)**2,((sigma**2)/mu))
    
                gene_copies=gene_copies+np.random.lognormal(np.log(mu),1.25)
    
        sample_value.append((10**6)*gene_copies)

    # simulate the flow and ammonia
    
    flow=[max(0,(10**6)*(np.sin(2*np.pi*t/365)+np.random.normal(1,0.5))) for t in sample_time]
    ammonia=[max(0,20*(np.cos(2*np.pi*t/365)+np.random.normal(0,0.5))) for t in sample_time]
    
    
    output['Date']=output['Date']+sample_date
    output['Site']=output['Site']+[site for t in sample_time]
    output['N1 Reported Value']=output['N1 Reported Value']+sample_value
    output['Flow (l/day)']=output['Flow (l/day)']+flow
    output['Ammonia (mg/l)']=output['Ammonia (mg/l)']+ammonia


pd.DataFrame(output).to_csv('../Data/wastewater_dummy.csv',index=False,sep='\t',encoding='utf-16')
