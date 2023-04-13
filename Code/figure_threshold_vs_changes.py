import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import pickle as pk
from datetime import datetime
from random import shuffle

from WWlibrary import model

df=pd.read_csv('../Data/population_count.csv')
sites=df['Site'].tolist()[:10]
#sites.remove('Grantown on Spey')

# measure days from this day
time_zero=datetime.strptime('2021-02-01', '%Y-%m-%d')
# Get list of dates for axes 
dates=['01/0'+str(i)+'/2021' for i in range(3,10)]\
    +['01/'+str(i)+'/2021' for i in range(10,13)]\
    +['01/0'+str(i)+'/2022' for i in range(1,10)]\
    +['01/'+str(i)+'/2022' for i in range(10,13)]\
    +['01/'+str(i)+'/2023' for i in range(1,3)]
    
dates_words=[datetime.strptime(str(d), '%d/%m/%Y').strftime('%b') for d in dates]
dates_numerical=[(datetime.strptime(str(d), '%d/%m/%Y')-time_zero).days for d in dates]


data_source='RNA count'#'RNA count' #'Admissions'#

# dummy files
ww_file='../Data/wastewater_dummy.csv'



number_of_changes={}
dispersion_list={}

for site in sites:



    if data_source=='RNA count':
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
    
            
    
    elif data_source=='Admissions':
    
        ######### Admissions #################################################
        df=pd.read_csv('../Site_level_admissions/'+site+'_admissions.csv')
        
        df['days_since_march1']=[(datetime.strptime(str(d), '%Y-%m-%d')-time_zero).days for d in df['Date'].tolist()]
        df=df[df['days_since_march1']>=0]
        df=df[(df['admissions']>0)]
        
        sample_value=[i/25 for i in df['admissions'].tolist()]
        sample_time=df['days_since_march1'].tolist()
        ####################################################################

    print(site)
    #choose wastewater only (no hospital for this one)
        
    # record the number of changes for each threshold value
    number_of_changes[site]=[]
    dispersion_list[site]=[]


    if site=='Seafield':
        fig = plt.figure(figsize=(10,4.5))
        gs = fig.add_gridspec(16, 2,width_ratios=[2,1])
        plt.subplots_adjust(hspace=0,wspace=0.2)
        
        
        plot=0
        ax=fig.add_subplot(gs[plot:plot+5,0]) 
        plt.scatter(sample_time,sample_value,c='m',s=2)
        ax.text(10,1.2,'Sample RNA copies (millions)')#,size=fs)
        ax.set_xlim([0,770])
        ax.set_ylim([0,1.5])
        ax.set_xticks([])
        ax.set_yticks([0,0.5,1])
        ax.text(-60,1.4,'A',size=20)

    for improvement_threshold in range(3):
        
        # create plot axis
        
        # load data
        record=pk.load(open('../pickles/different_thresholds/'+data_source+'/'+site+'_record_'+str(improvement_threshold)+'.p','rb'))    
        # convert to plotable model output
        params=record.pop()
        length=int((len(params)-1)/2)
        inflections=list(zip(params[1+length:],params[1:1+length]))
        start=min(sample_time)
        end=max(sample_time)
        mean=[np.exp(model(t,inflections)) for t in range(start,end+10)]
        dispersion=params[0]
        # confidence intervals
        upper=[I+(dispersion*I)**(1/2) for I in mean]
        lower=[I-(dispersion*I)**(1/2) for I in mean]
        
        
        
        #plt.figure()
        #plt.title(site+', smoothness='+str(improvement_threshold))
        # plot wastewater data
        
        if site=='Seafield' and (improvement_threshold==1 or improvement_threshold==2):
            plot=plot+5
            ax=fig.add_subplot(gs[plot:plot+5,0]) 
            plt.scatter(sample_time,sample_value,c='m',s=2)
        
            plt.plot(range(start,end+10),mean,c='k')
            plt.fill_between(range(start,end+10),upper,lower,color='k',alpha=0.2)
            ax.text(10,1.2,'Model with threshold='+str(improvement_threshold))#,size=fs)
            ax.set_xlim([0,770])
            ax.set_ylim([0,1.5])
            ax.set_xticks([])
            ax.set_yticks([0,0.5,1])
        #plt.savefig('../oct_figures/'+data_source+'_series_'+str(improvement_threshold)+'.png',format='png',dpi=300,bbox_inches='tight')
    
        

        number_of_changes[site].append(len(inflections))

        dispersion_list[site].append(dispersion)
        ##############################################
    
    if site=='Seafield':
        plt.xticks(dates_numerical,dates_words,rotation=60)

# years on axes

lh=-0.65
th=-0.86
year1=(datetime.strptime(str('01/01/2022'), '%d/%m/%Y')-time_zero).days
year2=(datetime.strptime(str('01/01/2023'), '%d/%m/%Y')-time_zero).days
ax.plot([15,year1-15],[lh,lh],c='k',linewidth=1,clip_on=False)
ax.text(year1/2-20,th,'2021')
ax.plot([year1+15,year2-15],[lh,lh],c='k',linewidth=1,clip_on=False)
ax.text((year1+year2)/2-20,th,'2022')
ax.plot([year2+15,770-15],[lh,lh],c='k',linewidth=1,clip_on=False)
ax.text((year2+770)/2-20,th,'2023')

#label

# remove the one from Seafield that didn't converge
number_of_changes['Seafield'][1]=None



colors=plt.cm.Paired(np.linspace(0, 1, 12))

ax = fig.add_subplot(gs[0:6,1])
n=0
for site in sites:
    plt.plot(range(1,3),number_of_changes[site][1:3],linewidth=3,c=colors[n],linestyle='-',marker='o',markersize=4,alpha=0.5)
    n=n+1
#plt.ylabel('Number of changes')
plt.xlabel('Threshold')
plt.ylim([0,80])
plt.xlim([0.5,14.5])
plt.text(0.5,84,'Number of change points')
ax.text(-2.6,77,'B',size=20)

ax = fig.add_subplot(gs[10:16,1])
n=0
for site in sites:
    plt.plot(number_of_changes[site][1:15],dispersion_list[site][1:15],c=colors[n],linewidth=3,linestyle='-',marker='o',markersize=4,alpha=0.5,label=site)
    n=n+1
plt.xlabel('Number of change points')
#plt.ylabel('Dispersion')
plt.ylim([0,0.1])
plt.xlim([2,75])
plt.yticks([0,0.05,0.1])
#plt.legend()
plt.text(2,0.105,'Dispersion, $D$')
ax.text(-13,0.11,'C',size=20)
plt.savefig('../figures/threshold_vs_changes.png',format='png',dpi=300,bbox_inches='tight')



