import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import pickle as pk
from datetime import datetime
from random import shuffle

from WWlibrary import model

# dummy files
ww_file='../Data/wastewater_dummy.csv'

df=pd.read_csv('../Data/population_count.csv')
sites=df['Site'].tolist()[:10]

data_source='RNA count'
# measure days from this day
time_zero=datetime.strptime('2021-02-01', '%Y-%m-%d')
# Get list of dates for axes 
dates=['01/0'+str(i)+'/2021' for i in range(2,10)]\
    +['01/'+str(i)+'/2021' for i in range(10,13)]\
    +['01/0'+str(i)+'/2022' for i in range(1,10)]\
    +['01/'+str(i)+'/2022' for i in range(10,13)]\
    +['01/'+str(i)+'/2023' for i in range(1,3)]
    

dates_words=[datetime.strptime(str(d), '%d/%m/%Y').strftime('%b') for d in dates]
dates_numerical=[(datetime.strptime(str(d), '%d/%m/%Y')-time_zero).days for d in dates]


fig = plt.figure(figsize=(10,3))
gs = fig.add_gridspec(16, 2, width_ratios=[2,1])
plt.subplots_adjust(hspace=0,wspace=0.2)


plot=-5


top=0.3
bottom=-0.3

earliness_list={}
growth_interval={}
growth_period_list={}

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

    print(len(sample_value))
        
    earliness_list[site]=[]
    growth_period_list[site]=[]
    
    for improvement_threshold in range(3): # change to 1,15
        good_threshold=True

        record=pk.load(open('../pickles/different_thresholds/'+data_source+'/'+site+'_record_'+str(improvement_threshold)+'.p','rb'))
        
        rate_series=[0]
        time_series=[0]
        
        growth_starts=[]
        decline_periods=[]
        growth_interval=[]
        
        # start of red/green background
        start=0
        for n in range(5,len(sample_value)-1):
            x=sample_time[:n]
            y=sample_value[:n]
            
            # end is the time of the nth data point
            end=max(x)
            params=record.pop(0)
           
            var=params[0]**2
            length=int((len(params)-1)/2)
            inflections=list(zip(params[1+length:],params[1:1+length]))
        
            # calculate the rate of increase
            # times
            t1,t2=inflections[-2][0],inflections[-1][0]
            #print(t1,t2)
            # heights
            h1,h2=inflections[-2][1],inflections[-1][1]
            rate=np.log(h1/h2)/(t1-t2)
        
            rate_series.append(rate)
            time_series.append(end)
            #increase since last time point?
            if rate>0 and rate_series[-2]<=0: 
                # store the info for later
                #growth_starts.append(end)
                start=end
            if rate<=0 and rate_series[-2]>0:
                # store the info for later
                growth_interval.append((start,end))
                start=end
        
        # if still in growth when reach the end add that period
        if rate>0:
            growth_interval.append((start,end))
        
        # add to list for plotting
        growth_period_list[site].append(len(growth_interval))
        
        if site=='Seafield' and improvement_threshold in [0,1,2]:          

            plot=plot+5
            ax=fig.add_subplot(gs[plot:plot+5,0])
            ax.set_xlim([0,770])
            ax.set_ylim([bottom,top])
            ax.set_xticks([])
            #ax.set_xticklabels(dates_words,rotation=60)
            ax.set_yticks([-0.2,0,0.2])
            ax.set_yticklabels([-0.2,0,0.2])
            ax.text(10,0.17,'Threshold='+str(improvement_threshold))#,size=fs)
    
            # rate should always be flat but with discontinuities
            new_time_series=[]
            new_rate_series=[0]
            for i in range(len(time_series)):
                new_time_series=new_time_series+[time_series[i],time_series[i]]
                new_rate_series=new_rate_series+[rate_series[i],rate_series[i]]
            new_time_series.append(time_series[i])
        
        
            ax.plot(new_time_series,new_rate_series,c='k',linewidth=1)
            plt.axhline(0,linestyle=':',color='k')

            # plot the growth periods as shades
            for (start,end) in growth_interval:
                plt.fill_between([start,end],[bottom,bottom],[top,top], linewidth=0, color='r',alpha=0.2)

            #for j in range(len(rate_series)):
            #    print(time_series[j],rate_series[j])
            if improvement_threshold==4:
                #ax.set_ylabel('Growth rate')
                ax.text(0,0.35,'Real-time growth rate estimate')

plt.xticks(dates_numerical,dates_words,rotation=60)

lh=-0.66
th=-0.8
# years on axes
year1=(datetime.strptime(str('01/01/2022'), '%d/%m/%Y')-time_zero).days
year2=(datetime.strptime(str('01/01/2023'), '%d/%m/%Y')-time_zero).days
ax.plot([15,year1-15],[lh,lh],c='k',linewidth=1,clip_on=False)
ax.text(year1/2-20,th,'2021')
ax.plot([year1+15,year2-15],[lh,lh],c='k',linewidth=1,clip_on=False)
ax.text((year1+year2)/2-20,th,'2022')
ax.plot([year2+15,770-15],[lh,lh],c='k',linewidth=1,clip_on=False)
ax.text((year2+770)/2-20,th,'2023')

# label
ax.text(-100,1.45,'A',size=20)

# remove the one from Seafield that didn't converge
growth_period_list['Seafield'][0]=None

# Paired is the name of the cmap
colors=plt.cm.Paired(np.linspace(0, 1, 12))

ax = fig.add_subplot(gs[:,1])
n=0
for site in sites:
    plt.plot(range(3),growth_period_list[site],linewidth=3,linestyle='-',c=colors[n],marker='o',markersize=4,alpha=0.5)
    n=n+1
plt.ylim([0,40])
plt.xlim([0.5,14.5])
plt.xlabel('Threshold')

plt.text(0.5,41,'Number of growth periods')
ax.text(-2.6,39,'B',size=20)
plt.savefig('../figures/threshold_vs_waves.png',format='png',dpi=300,bbox_inches='tight')



