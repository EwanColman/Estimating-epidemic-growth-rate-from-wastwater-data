import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import pickle as pk
from datetime import datetime
from random import shuffle
from random import random

from WWlibrary import model
from WWlibrary import fit_the_function
from WWlibrary import likelihood
from scipy import stats



# dummy files
ww_file='../Data/wastewater_dummy.csv'

improvement_threshold=8

# measure days from this day
time_zero=datetime.strptime('2021-02-01', '%Y-%m-%d')


df=pd.read_csv('../Data/population_count.csv')

population={}
for i,row in df.iterrows():
    population[row['Site']]=row['Population']






mean_likelihoods=[]
log_population=[]
number_of_points=[]

n=0
for site in population:
    
    df=pd.read_csv(ww_file,sep='\t',encoding='utf-16').fillna(0)
    
    # add it to the dataframe
    df['days_since_march1']=[(datetime.strptime(str(d), '%d/%m/%Y')-time_zero).days for d in df['Date'].tolist()]
    df=df[df['days_since_march1']>=0]
    df=df[(df['N1 Reported Value']>0)&(df['N1 Reported Value']<2000000)]
    
    last_day=max(df['days_since_march1'])
    sorted_df=df.sort_values('days_since_march1',ascending=True)
    
    site_df=sorted_df[sorted_df['Site']==site]
    
    # restrict analysis to those with 100+ data points
    if len(site_df)>5:
    
        sample_time=site_df['days_since_march1'].tolist()
        # divide by a million
        sample_value=[i/1000000 for i in site_df['N1 Reported Value'].tolist()]

        samples=len(sample_value)
    
        #############################################
    
        record=pk.load(open('../pickles/different_sites/RNA count/'+site+'_record_'+str(improvement_threshold)+'.p','rb'))
            
        params=record.pop()
        # update the params list of variables
        length=int((len(params)-1)/2)
        inflections=list(zip(params[1+length:],params[1:1+length]))
        dispersion=params[0]
        
        
        l=likelihood(params,sample_time,sample_value)
        
        mean_likelihoods.append(l/len(sample_time))
        
        # pops as a list
        log_population.append(np.log(population[site]))
        n=n+1
        print(n,site,'population:',str(population[site])+',',len(sample_time),'data points, likelihood=',l/len(sample_time))

        number_of_points.append(len(sample_time))

        # plot to check convergence
        # create plot axis    
        plt.figure()
        if site=='Banff':
            plt.title(site+', smoothness='+str(improvement_threshold))
            # plot wastewater data
            plt.scatter(sample_time,sample_value,s=10)
            # get the mean of the model        
            mean=[np.exp(model(t,inflections)) for t in range(min(sample_time),max(sample_time)+10)]
             
            upper=[I+(dispersion*I)**(1/2) for I in mean]
            lower=[I-(dispersion*I)**(1/2) for I in mean]
    
            plt.plot(range(min(sample_time),max(sample_time)+10),mean)
            plt.fill_between(range(min(sample_time),max(sample_time)+10),upper,lower,alpha=0.3)
            plt.xlabel('Time (days)')
            plt.ylabel('N1 Reported Value')
    else:
        print('REMOVED',site,len(site_df),'data points')

pearson, p=stats.pearsonr(log_population,mean_likelihoods)
print(pearson,p)

fig=plt.figure(figsize=(10,3))

fig.add_subplot(121)

plt.scatter(log_population,mean_likelihoods,s=1)
plt.title('Mean log likelihood, correlation='+str(round(pearson,2))+', p='+str(round(p,2)))
plt.xlabel('Population')
plt.xticks([i*np.log(10) for i in range(3,6)],['$10^{'+str(i)+'}$' for i in range(3,6)])

pearson, p=stats.pearsonr(number_of_points,mean_likelihoods)
print(pearson,p)

fig.add_subplot(122)
plt.scatter(number_of_points,mean_likelihoods,s=1)
plt.title('Mean log likelihood, correlation='+str(round(pearson,2))+', p='+str(round(p,2)))
plt.xlabel('Number of data points')
plt.savefig('../figures/population_vs_model_fit.png',format='png',dpi=300,bbox_inches='tight')



