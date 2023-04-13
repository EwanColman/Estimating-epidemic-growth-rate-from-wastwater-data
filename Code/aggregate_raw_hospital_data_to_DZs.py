import pandas as pd
import time
from datetime import datetime
import matplotlib.pyplot as plt
import pickle as pk

data_date='dummy'

start=time.time()
df = pd.read_csv('../Data/smr01_'+data_date+'.csv').fillna('NAN')

print('time to read csv:',int(time.time()-start))




df=df[df['admission_date']>=20210201]

for c in df.columns:
    print(c)
print()

#filter down to just the covids
df=df[(df['main_condition'].isin(['U072']))
      |(df['other_condition_1'].isin(['U071','U072']))
      |(df['other_condition_2'].isin(['U071','U072']))
      |(df['other_condition_3'].isin(['U071','U072']))
      |(df['other_condition_4'].isin(['U071','U072']))
      |(df['other_condition_5'].isin(['U071','U072']))
    ]

# only as main condition
#df=df[df['main_condition']=='U071']

# time in days
time_zero=datetime.strptime('2020-03-01', '%Y-%m-%d')
df['time_in_days']=[(datetime.strptime(str(d), '%Y%m%d')-time_zero).days for d in df['admission_date']] 


id_list=df['id'].tolist()
print('Before:',len(id_list),'IDs',len(set(id_list)),'unique')

start=time.time()

# get a dataframe only including the first time the id appears in the data
first_admission_df=df[['admission_date','datazone_2011','main_condition','id','time_in_days']].groupby(['id'],as_index=False).min()

# store the first admission time of all the id's
first_admission={}
date=first_admission_df['time_in_days'].tolist()
ID=first_admission_df['id'].tolist()
while ID:
    first_admission[ID.pop()]=date.pop()

# add a column to the full dataframe giving the time of first admission of the id in the row
df['time_of_first_admission']=[first_admission[i] for i in df['id'].tolist()]
# get a df of admissions that happen at least 60 later than the first
second_admission_df=df[df['time_in_days']-df['time_of_first_admission']>60]
# take the first appearance of these second admissions
second_admission_df=second_admission_df[['admission_date','datazone_2011','main_condition','id','time_in_days']].groupby(['id'],as_index=False).min()
print(len(second_admission_df))

# should probably remove third, fourth admissions too but not coded that yet 
# Some checks suggest that it won;t make much difference

# put the first and second admissions back together
df=pd.concat([first_admission_df,second_admission_df])

print('time to group by date:',int(time.time()-start))

print(df.head())
print(df.columns)

id_list=df['id'].tolist()
print('After:',len(id_list),'IDs',len(set(id_list)),'unique')



date_list=set(df['admission_date'].tolist())

print('earliest:',min(date_list),'latest:',max(date_list))

# aggregate to datazone
#covid_df = pd.pivot_table(df[df['main_condition'].isin(['U071'])], index=['admission_date','datazone_2011'], aggfunc=len).reset_index()
covid_df = pd.pivot_table(df, index=['admission_date','datazone_2011'], aggfunc=len).reset_index()


# get DZs
DZ_list=list(set(covid_df['datazone_2011'].tolist()))

# time stuff
time_zero=datetime.strptime('2021-02-01', '%Y-%m-%d')
max_day=(datetime.strptime(str(max(date_list)), '%Y%m%d')-time_zero).days
print('max day',max_day,max(date_list))


admissions={}
for DZ in DZ_list:
    admissions[DZ]=[0 for i in range(max_day+1)]


date=covid_df['admission_date'].tolist()
number=covid_df['main_condition'].tolist()
datazone=covid_df['datazone_2011'].tolist()
i=0
while date:
    i=i+1
    if i%100000==0:
        print(i)
    d=date.pop()
    n=number.pop()
    DZ=datazone.pop()
    day_numerical=(datetime.strptime(str(d), '%Y%m%d')-time_zero).days
    admissions[DZ][day_numerical]=n
    
total=[sum([admissions[DZ][i] for DZ in DZ_list]) for i in range(max_day+1)]
    

plt.plot(total)
    
pk.dump(admissions,open('../pickles/DZ_level_admissions/admissions_'+data_date+'.p','wb'))
