import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import pickle as pk
from datetime import datetime
from random import shuffle

from WWlibrary import model
from WWlibrary import fit_the_function

# uncomment for using actiual data
#ww_file='../Data/RNA Monitoring Project - SW Sample Number, N1 Description, N1 Reported Value, N1 Sample 1, N1 Sa.csv'

# dummy files
ww_file='../Data/wastewater_dummy.csv'

# choose the data you want to look at
data_source='Admissions'#'RNA count' #'Admissions'#

# measure days from this day
time_zero=datetime.strptime('2021-02-01', '%Y-%m-%d')

# choose the site
df=pd.read_csv('../Data/population_count.csv')
sites=df['Site'].tolist()
# choose the top n sites (ordered by population)
for site in sites:

    if data_source=='RNA count':
        #print('Wastewater')
        ##################### Wastewater ###################################
        df=pd.read_csv(ww_file,sep='\t',encoding='utf-16').fillna(0)
        
        # add it to the dataframe
        df['days_since_march1']=[(datetime.strptime(str(d), '%d/%m/%Y')-time_zero).days for d in df['Date'].tolist()]
        df=df[df['days_since_march1']>=0]
        df=df[(df['N1 Reported Value']>0)&(df['N1 Reported Value']<2000000)]
        
        last_day=max(df['days_since_march1'])
        sorted_df=df.sort_values('days_since_march1',ascending=True)
        
        site_df=sorted_df[sorted_df['Site']==site]
        
        sample_time=site_df['days_since_march1'].tolist()
        # divide by a million
        sample_value=[i/1000000 for i in site_df['N1 Reported Value'].tolist()]
    
        plt.scatter(sample_time,sample_value)
            
    
    elif data_source=='Admissions':
    
        ######### Admissions #################################################
        df=pd.read_csv('../Site_level_admissions/'+site+'_admissions.csv')
        
        df['days_since_march1']=[(datetime.strptime(str(d), '%Y-%m-%d')-time_zero).days for d in df['Date'].tolist()]
        df=df[df['days_since_march1']>=0]
        df=df[(df['admissions']>0)]
        
        sample_value=[i for i in df['admissions'].tolist()]
        sample_time=df['days_since_march1'].tolist()
        ####################################################################
    
    
    number_of_data_points=len(sample_value)
    
    #############################################
    
    for improvement_threshold in range(3):
        print(data_source,site,'threshold=',improvement_threshold)
        
        record=fit_the_function(sample_time,sample_value,improvement_threshold)
        # get the last recorded param values
        params=record[-1]
        ###################
        # update the params list of variables
        length=int((len(params)-1)/2)
        inflections=list(zip(params[1+length:],params[1:1+length]))
        dispersion=params[0]
    
        # create plot axis    
        plt.figure()
        plt.title(site+', smoothness='+str(improvement_threshold))
        # plot wastewater data
        plt.scatter(sample_time,sample_value,s=10)
        # get the mean of the model
        
        mean=[np.exp(model(t,inflections)) for t in range(min(sample_time),max(sample_time)+10)]
         
        upper=[I+(dispersion*I)**(1/2) for I in mean]
        lower=[I-(dispersion*I)**(1/2) for I in mean]
        
        
        
        plt.plot(range(min(sample_time),max(sample_time)+10),mean)
        plt.fill_between(range(min(sample_time),max(sample_time)+10),upper,lower,alpha=0.3)
        #plt.ylim([0,0.7])
        #plt.yscale('log')
        #plt.xlim([0,150])
        
        plt.xlabel('Time (days)')
        plt.ylabel('N1 Reported Value')
        
        # store the rolling results
        pk.dump(record,open('../pickles/different_thresholds/'+data_source+'/'+site+'_record_'+str(improvement_threshold)+'.p','wb'))

