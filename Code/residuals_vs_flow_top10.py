import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import pickle as pk
from datetime import datetime
from random import shuffle

from WWlibrary import model
from WWlibrary import fit_the_function
from scipy import stats


# dummy files
ww_file='../Data/wastewater_dummy.csv'
# measure days from this day
time_zero=datetime.strptime('2021-02-01', '%Y-%m-%d')

# choose the site
df=pd.read_csv('../Data/population_count.csv')
# get top ten
sites=['Banff']+df['Site'].tolist()


fig = plt.figure(figsize=(10,12))
num_rows=5
gs = fig.add_gridspec(5, 2)
plt.subplots_adjust(hspace=1,wspace=0.2)

combined_x=[]
combined_y=[]

improvement_threshold=8
ax_num=0
count=0
sites_with_flow_data=0
while sites:
    
    site=sites.pop(0)

    #print('Wastewater')
    ##################### Wastewater ###################################
    df=pd.read_csv(ww_file,sep='\t',encoding='utf-16').fillna(0)
    
    # add it to the dataframe
    df['days_since_march1']=[(datetime.strptime(str(d), '%d/%m/%Y')-time_zero).days for d in df['Date'].tolist()]
    df=df[df['days_since_march1']>=0]
    df=df[(df['N1 Reported Value']>0)&(df['N1 Reported Value']<2000000)]
    df=df[df['Flow (l/day)']>0]
    
    
    last_day=max(df['days_since_march1'])
    sorted_df=df.sort_values('days_since_march1',ascending=True)
    site_df=sorted_df[sorted_df['Site']==site]
    
    sample_time=site_df['days_since_march1'].tolist()
    # divide by a million
    sample_value=[i/1000000 for i in site_df['N1 Reported Value'].tolist()]
    
    flow_value=[i/1000000 for i in site_df['Flow (l/day)'].tolist()]
    ammonia_value=[i for i in site_df['Ammonia (mg/l)'].tolist()]
    
    # enter here the minimum size to be considered worth including (in the paper its 100)
    if len(sample_time)>10:
        sites_with_flow_data=sites_with_flow_data+1
        # load data
        record=pk.load(open('../pickles/different_sites/RNA count/'+site+'_record_'+str(improvement_threshold)+'.p','rb'))    
        # convert to plotable model output
        params=record.pop()
        length=int((len(params)-1)/2)
        inflections=list(zip(params[1+length:],params[1:1+length]))
        start=min(sample_time)
        end=max(sample_time)
        mean=[np.exp(model(t,inflections)) for t in range(0,end+10)]
        dispersion=params[0]
        
        
        
        
        residuals=[]
        z_scores=[]
        for i in range(len(sample_value)):
            value=sample_value[i]
            time=sample_time[i]
            
            z_score=(value-mean[time])/((dispersion*mean[time])**(1/2))
            
            residuals.append(value-mean[time])
            z_scores.append(z_score)
            
        # choose flow or ammonia
        x=flow_value
        #x=ammonia_value
        y=z_scores
            
        pearson, p=stats.pearsonr(x,y)
        
        slope, intercept, r, p, se=stats.linregress(x,y)
        
        if ax_num<10:
            #print(ax_num%num_rows,int(ax_num/num_rows))
            ax=fig.add_subplot(gs[ax_num%num_rows,int(ax_num/num_rows)])
            
            ax.set_title(site+', p='+str('%.1g' % p)+', slope='+str('%.1g' % slope))
            plt.xlabel('Flow (ML/day)')
            if ax_num==1:
                plt.ylabel('z-score')
            plt.axhline(0,linewidth=1,c='k')
            plt.scatter(x,y,s=1)
            
            #line of regression
            x_vals=[min(x),max(x)]
            y_vals=[intercept+slope*i for i in x_vals]
            plt.plot(x_vals,y_vals)
            
            ax_num=ax_num+1
    
        if p<0.05:
            print(count,site+', 100 ML of flow reduces mean by '+str(slope*100)+' standard deviations')    
            count=count+1
        


print(sites_with_flow_data,'sites with flow data')
plt.savefig('../figures/flow_vs_residuals_top_10.png',format='png',dpi=300,bbox_inches='tight')

