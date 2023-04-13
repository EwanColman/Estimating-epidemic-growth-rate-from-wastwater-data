import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import pickle as pk
from datetime import datetime
from random import shuffle

from WWlibrary import model
from scipy import stats

# dummy files
ww_file='../Data/wastewater_dummy.csv'

data_source='Admissions'

threshold={'Admissions':2,#8 in paper
           'RNA count':2}#8 in paper

df=pd.read_csv('../Data/population_count.csv')
sites=df['Site'].tolist()[:10]

# measure days from this day
time_zero=datetime.strptime('2021-02-01', '%Y-%m-%d')
# dates used in paper
#hospital_data_start_time=(datetime.strptime('2021-05-01', '%Y-%m-%d')-time_zero).days
#end_time=(datetime.strptime('2022-07-29', '%Y-%m-%d')-time_zero).days

# dates used for dummy data
hospital_data_start_time=(datetime.strptime('2021-03-01', '%Y-%m-%d')-time_zero).days
end_time=(datetime.strptime('2021-03-20', '%Y-%m-%d')-time_zero).days

# Get list of dates for axes 
dates=['01/0'+str(i)+'/2021' for i in range(5,10)]\
    +['01/'+str(i)+'/2021' for i in range(10,13)]\
    +['01/0'+str(i)+'/2022' for i in range(1,10)]
    #+['01/'+str(i)+'/2022' for i in range(10,10)]

dates_words=[datetime.strptime(str(d), '%d/%m/%Y').strftime('%b') for d in dates]
dates_numerical=[(datetime.strptime(str(d), '%d/%m/%Y')-time_zero).days for d in dates]


fig = plt.figure(figsize=(10,3))
gs = fig.add_gridspec(4, 2,height_ratios=[4,2,3,3],width_ratios=[2,1])
plt.subplots_adjust(hspace=0,wspace=0.3)



correlation={}
p_val={}

for site in sites:
    #print()
    #print(site)
    daily_rate_series={}
    for data_source in ['RNA count','Admissions']:

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
        
        
        elif data_source=='Admissions':
            #print('Admissions')
            ######### Admissions #################################################
            df=pd.read_csv('../Site_level_admissions/'+site+'_admissions.csv')

            df['days_since_march1']=[(datetime.strptime(str(d), '%Y-%m-%d')-time_zero).days for d in df['Date'].tolist()]
            df=df[df['days_since_march1']>=0]
            df=df[(df['admissions']>0)]
            
            sample_value=df['admissions'].tolist()
            sample_time=df['days_since_march1'].tolist()
            
            print(site,max(df['Date'].tolist()),max(sample_time))
            ####################################################################
        
        record=pk.load(open('../pickles/different_thresholds/'+data_source+'/'+site+'_record_'+str(threshold[data_source])+'.p','rb'))
        
        rate_series=[0]
        time_series=[0]
        
        start=0
        for n in range(5,len(sample_value)):
            x=sample_time[:n]
            y=sample_value[:n]
        
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
        
        
        
        # convert rates into a daily series
        daily_rate_series[data_source]=[]
       
        # interpolate
        for i in range(len(time_series)-1):
            for d in range(time_series[i],time_series[i+1]):
                daily_rate_series[data_source].append(rate_series[i])
        
        #if its seafield then plot the time series
        if site=='Seafield':
            if data_source=='RNA count':
                ax=fig.add_subplot(gs[2,0])
                ax.set_xticks([])
                
                ax.text(hospital_data_start_time+10,0.15,'From waste water')
                ax.text(hospital_data_start_time,0.36,'Real-time growth rate estimate')
                ax.text(hospital_data_start_time-40,0.4,'B',size=20)
                
            else:
                # plot raw data and best fit
                ax=fig.add_subplot(gs[0,0]) 
                plt.scatter(sample_time,sample_value,c='m',s=2)
                plt.xlim([hospital_data_start_time,end_time+20])
        
                #params=record[-1]
                length=int((len(params)-1)/2)
                inflections=list(zip(params[1+length:],params[1:1+length]))
                start=min(sample_time)
                end=max(sample_time)
                mean=[np.exp(model(t,inflections)) for t in range(start,end+10)]
                dispersion=params[0]
                upper=[I+(dispersion*I)**(1/2) for I in mean]
                lower=[I-(dispersion*I)**(1/2) for I in mean]
       
                plt.plot(range(start,end+10),mean,c='k')
                plt.fill_between(range(start,end+10),upper,lower,color='k',alpha=0.2)
                
                # fake legend
                plt.plot([hospital_data_start_time+10,hospital_data_start_time+20],[27,27],c='k')
                plt.fill_between([hospital_data_start_time+10,hospital_data_start_time+20],[24,24],[30,30],color='k',alpha=0.2)
                ax.text(hospital_data_start_time+25,24,'Model with threshold=8')
                plt.scatter([hospital_data_start_time+17],[35],c='m',s=5)
                ax.text(hospital_data_start_time+25,32,'Hospital admissions')
                ax.text(hospital_data_start_time-40,40,'A',size=20)
                
                ax.set_xlim([hospital_data_start_time,end_time+20])
                ax.set_ylim([0,40])
                ax.set_xticks([])
                ax.set_yticks([10,20,30])
                #ax.set_ylabel('Admissions')
                
                ax=fig.add_subplot(gs[3,0])
                ax.text(hospital_data_start_time+10,0.15,'From hospital admissions')
                ax.set_xticks(dates_numerical)
                ax.set_xticklabels(dates_words,rotation=60)
                #ax.set_ylabel()
            
            ax.set_xlim([hospital_data_start_time,end_time+20])
            ax.set_ylim([-0.3,0.3])
            
            #ax.set_xticklabels(dates_words,rotation=60)
            ax.set_yticks([-0.2,0,0.2])
            ax.set_yticklabels([-0.2,0,0.2])
    

            ax.plot(range(hospital_data_start_time,end_time),daily_rate_series[data_source][hospital_data_start_time:end_time],c='k',linewidth=1)
            plt.axhline(0,linestyle=':',color='k')
    

    
    correlation[site]=[]
    p_val[site]=[]
    lag=[]
    for i in range(10):
        l=5-i
        # trim series and include lag
        #if l>0:
        trimmed_Admissions=daily_rate_series['Admissions'][hospital_data_start_time:end_time]
        trimmed_RNA=daily_rate_series['RNA count'][hospital_data_start_time-l:end_time-l]
        # else:
        #     trimmed_Admissions=daily_rate_series['Admissions'][0:400]
        #     trimmed_RNA=daily_rate_series['RNA count'][-l:400-l]
        pearson, p=stats.pearsonr(trimmed_RNA,trimmed_Admissions)
        #print(lag,p)
        lag.append(l)
        correlation[site].append(pearson)
        p_val[site].append(p)
        
 
# years on axes
lh=-0.78
th=-0.95
year=(datetime.strptime(str('01/01/2022'), '%d/%m/%Y')-time_zero).days
ax.plot([hospital_data_start_time+10,year-10],[lh,lh],c='k',linewidth=1,clip_on=False)
ax.text((year+hospital_data_start_time)/2-10,th,'2021')
ax.plot([year+10,end_time+10],[lh,lh],c='k',linewidth=1,clip_on=False)
ax.text((year+end_time)/2-10,th,'2022')

 
    
ax=fig.add_subplot(gs[0:4,1])
plt.axvline(0,linewidth=1,c='k',zorder=0)
plt.axhline(0,linewidth=1,c='k',zorder=0)

# get colours
colors=plt.cm.Paired(np.linspace(0, 1, 12))

lags=[]
n=0
for site in sites:
    
    ax.plot(lag,correlation[site],label=site,linewidth=3,linestyle='-',c=colors[n],markersize=2,alpha=0.5,zorder=1)
    
    best_lag=0
    # find max 
    biggest=0
    for i in range(10):
        l=5-i
        if correlation[site][i]>biggest:
            best_lag=l
            biggest=correlation[site][i]
            
        
    plt.scatter([best_lag],[-0.5],s=8,c=colors[n])
    n=n+1
    print(site,'lag=',best_lag)
    
    lags.append(best_lag)
    
print(lags)
print(np.mean(lags))
    
ax.set_ylim([-0.52,0.68])
ax.set_xlim([-50,50])
#ax.set_ylabel('Correlation coefficient')
plt.text(-50,0.71,'Correlation coefficient')
ax.set_xlabel('Lag from wastewater to admissions')
ax.text(-75,0.68,'C',size=20)
                
        

mean_correlation=[sum([correlation[site][i]/10 for site in sites]) for i in range(10)]
# for i in range(100):
#     l=50-i
#     if total_correlation[i]>biggest:
#         best_lag=l
#         biggest=total_correlation[i]

plt.plot([5-i for i in range(10)],mean_correlation,c='k',linestyle='--',linewidth=4)

plt.savefig('../figures/hospital_vs_wastewater.png',format='png',dpi=300,bbox_inches='tight')

print('Admissions: '+str(threshold['Admissions'])+', RNA count: '+str(threshold['RNA count']), ', total r:'+str(biggest)+', lag:'+str(best_lag))