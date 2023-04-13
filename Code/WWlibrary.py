import numpy as np
import pickle as pk
from datetime import datetime
from random import shuffle

# takes the inflection points and returns the log of the mean at time t
def model(t,inflection):
    # find which interval t is in
    # if its out of the range
    if t>inflection[-1][0]:
        k=len(inflection)-1
    else:
        k=1
        while t>inflection[k][0]:
            #print(t,k,inflection[k])
            k=k+1
        
    x0=inflection[k-1][0]
    y0=np.log(inflection[k-1][1])
    
    x1=inflection[k][0]
    y1=np.log(inflection[k][1])
    
    m=(y0-y1)/(x0-x1)
    c=(y1*x0-y0*x1)/(x0-x1)
        
    return m*t+c

# calculates the likelihood of the model given the data arrays x and y
def likelihood(params,time,rna):
    dispersion=params[0]
    length=int((len(params)-1)/2)
    inflections=list(zip(params[1+length:],params[1:1+length]))
    #inflection=[60]+[end]
    #rate=[0.08,-0.06]
    #time=x
    #rna=y
        
    # loop over all te data points and add up theor likelihoods
    log_likelihood=0
    for i in range(len(time)):
        # time of sampling
        t=time[i]
        # amont sampled
        v=rna[i]
        # modeeled number of infected
        mu=np.exp(model(t,inflections))
        var=mu*dispersion
        log_normal=-(1/2)*np.log(np.pi*2*var)-(1/2)*((v-mu)**2)/var
        # add to total
        log_likelihood=log_likelihood+log_normal
    return log_likelihood

# creates the ranges of parameter values for the algorithm to explore
def convert_to_ranges(inflections,dispersion,start,end,upper_limit):
    
    # PART 1: dispersion
    # 0 represents the index of dispersion
    # create the variable ranges around the initial values given
    variable_range={0:[i*upper_limit/200 for i in range(1,20)]}
    # create an index for each variable 
    index={0:int(dispersion*200/upper_limit)-1}

    ####################################
    # only give a range to the variables we want to optimize over
    # i.e. only the most recent 5(?) time points 
    ####################################

    # PART 2: model values at the inflection points
    M=1.5*upper_limit # max value
    m=0.01 # min value
    N=100 # number of values
    
    n=0
    # list of heights
    heights=[i[1] for i in inflections]
    
    for h in heights:
        n=n+1
        
        # if its in the last 5 time points
        if n>len(heights)-5: # change to something more general

            variable_range[n]=[m+i*(M-m)/N for i in range(N)]
            index[n]=min(N-1,int(N*(h-m)/(M-m)))
            
        else:
            # otherwise only include the established value
            variable_range[n]=[h]
            index[n]=0

    # PART 3: time points of inflecions

    # get list of time points for the inflections
    inflection_times=[i[0] for i in inflections]

    for i in range(len(inflection_times)):
        n=n+1
        # start and end have no range
        if i==0:
            lower_bound=start
            upper_bound=start+1
        elif i==len(inflection_times)-1:
            lower_bound=end
            upper_bound=end+1
        # the first adjustable inflection point 
        elif i==1:
            lower_bound=start+1
            # upper bound is mid-point of the first and second current times
            upper_bound=int((inflection_times[i]+inflection_times[i+1])/2)
        elif i==len(inflection_times)-2:
            # get the lower bound that doen't intersect with the previous
            lower_bound=int((inflection_times[i-1]+inflection_times[i])/2)
            upper_bound=end
            
        else:
            # get the lower bound that doen't intersect with the previous
            lower_bound=int((inflection_times[i-1]+inflection_times[i])/2)
            upper_bound=int((inflection_times[i]+inflection_times[i+1])/2)
            
        variable_range[n]=[j for j in range(lower_bound,upper_bound)]  
        
        # where does the current inflection point appear in the list?
        #print(inflection_times[i],lower_bound,upper_bound)
        index[n]=min(inflection_times[i]-lower_bound,len(variable_range[n])-1)
    
        # if its old then overwrite
        if n<len(heights)+len(inflection_times)-5:
            # only include the established value in the list of possibilities
            variable_range[n]=[inflection_times[i]]
            index[n]=0
    
    return variable_range,index

# finds the local optima
def hill_climb(variable_range,index,x,y):
    # create these lists
    params=[] # the parameter values
    free_variables=[] # the variables that can be adjusted
    # add the current value for each variable
    for variable in variable_range:  
        # want to keep the old result as the starting point for the new
        # find the index that corresponds to the variable
        params.append(variable_range[variable][index[variable]])
        # add the free ones to the list
        if len(variable_range[variable])>1:# 
            free_variables.append(variable)
    
    # to start the best oarams are the old params
    E_old=likelihood(params,x,y)
    
    E_best=E_old.copy()
    best_params=params.copy()

    # now look for something better
    maxima_found=False

    ###
    while not maxima_found:

        var_list=free_variables.copy()

        #loop over every free variable, find the best perturbation
        for variable in var_list:
            #print(variable,'of',len(var_list))
            # copy the dictionary of indexes so we can perturb it
            new_index=index.copy()
            # first try the down direction
            new_index[variable]=max(0,new_index[variable]-1)
            # new params
            perturbed_params=[variable_range[var][new_index[var]] for var in variable_range]
            
            # one last test, if it fails then do not mark as an improvement
            # inflections have to be 7+ days apart
            # get inflection times from params
            valid=True
            length=int((len(perturbed_params)-1)/2)
            inflection_times=perturbed_params[1+length:]
            for i in range(1,len(inflection_times)):
                # if any of them are too close then it becomes invalid
                if inflection_times[i]-inflection_times[i-1]<7:
                    valid=False
            
            
            # new max likelihood
            E_new=likelihood(perturbed_params,x,y)
            
            # is it an improvement
            if E_new>E_best and valid:
                E_best=E_new
                best_params=perturbed_params.copy()
                best_index=new_index.copy()
        
            else:
                # if not then try the other direction
                new_index=index.copy()
                # try the up direction
                new_index[variable]=min(new_index[variable]+1,len(variable_range[variable])-1)
            
                perturbed_params=[variable_range[var][new_index[var]] for var in variable_range]
                
                # one last test, if it fails then do not mark as an improvement
                # inflections have to be 7+ days apart
                # get inflection times from params
                valid=True
                length=int((len(perturbed_params)-1)/2)
                inflection_times=perturbed_params[1+length:]
                for i in range(1,len(inflection_times)):
                    # if any of them are too close then it becomes invalid
                    if inflection_times[i]-inflection_times[i-1]<7:
                        valid=False

                #
                E_new=likelihood(perturbed_params,x,y)
   
                if E_new>E_best and valid:
                    E_best=E_new
                    best_params=perturbed_params.copy()
                    best_index=new_index.copy()
            
                
            
        # after all that, is the best improvement an improvement?
        if E_old<E_best and valid:
            params=best_params.copy()
            index=best_index.copy()
            E_old=E_best
            
        # otherwise there is no way to improve
        else:
            maxima_found=True

    return params

# run the procedure to sequentially add data, re-fit, and test candidate inflections
def fit_the_function(sample_time,sample_value,improvement_threshold):
    # inital conditions for very start (first 23 points)
    

    m=5
    x=sample_time[:m]
    y=sample_value[:m]
    
    # find the scale of the data
    upper_limit=max(sample_value)
    
    # choose an initial dispersion in the middle of the range
    dispersion=upper_limit/20
    
    end=max(x)+1
    start=min(x)
    
    inflections=[(start,y[0]),
                  (end,y[-1])]
    
    # 'record' will be the record of model outputs at each time point
    # in future would be better to get a range of probabilities of each outcome
    # so we can attach a measure of confidence 
    record=[]
    #################
    # up to 158
    for n in range(m,len(sample_time)):
        
        #print('Threshold='+str(improvement_threshold)+', First '+str(n)+' samples')
        x=sample_time[:n]
        y=sample_value[:n]
        
        # new sample so need to extend the range
        # get new end
        end=max(x)+1
        new_end=(end,np.exp(model(end,inflections)))
        # remove the last point
        inflections.pop()
        # and repace with the later one
        inflections.append(new_end)
        
        # use the current inflections array to create variabl ranges
        variable_range,index=convert_to_ranges(inflections,dispersion,start,end,upper_limit)
    
        ## run the fitting function
        #params_without=simulated_annealing(variable_range,index)
        params_without=hill_climb(variable_range,index,x,y)
        # get the likelihood
        likelihood_without_inflection=likelihood(params_without,x,y)
        #print('without new inflection:',likelihood_without_inflection)
        # get new inflections and sd
        length=int((len(params_without)-1)/2)
        inflections=list(zip(params_without[1+length:],params_without[1:1+length]))

        # add the new initial value for inflection point
        # the last one
        last_inflection=max(i[0] for i in inflections if i[0]!=end) 
        # initial guess: insert at the mid point between the last inflection and the end
        new_inflection=int((end+last_inflection)/2)
        # add this point 
        inflections.append((new_inflection,np.exp(model(new_inflection,inflections))))
        # probably don;t need to sort it, but do it anyway (delete later?)
        inflections =sorted(inflections, key = lambda item: item[0])
        # process the inflections array
        variable_range,index=convert_to_ranges(inflections,dispersion,start,end,upper_limit)
        
        params_with=hill_climb(variable_range,index,x,y)
        # get the likelihood
        likelihood_with_inflection=likelihood(params_with,x,y)
        
        #print('with new inflection:',likelihood_with_inflection)
        
        # choose which one to go with
        if likelihood_with_inflection-likelihood_without_inflection>improvement_threshold:
            params=params_with
            #print('Inflection added')
        else:
            params=params_without
    
        # update the params list of variables
        length=int((len(params)-1)/2)
        inflections=list(zip(params[1+length:],params[1:1+length]))
        dispersion=params[0]
    
        # store the model so that it can be used for prediction
        record.append(params)

    return record


##################

# calculates the likelihood of the model given the data arrays x and y
def likelihood2(params,time,rna):
    dispersion=params[0]
    length=int((len(params)-1)/2)
    inflections=list(zip(params[1+length:],params[1:1+length]))
    #inflection=[60]+[end]
    #rate=[0.08,-0.06]
    #time=x
    #rna=y
        
    # loop over all te data points and add up theor likelihoods
    log_likelihood=0
    for i in range(len(time)):
        # time of sampling
        t=time[i]
        # amont sampled
        v=rna[i]
        # modeeled number of infected
        mu=np.exp(model(t,inflections))
        var=mu*dispersion
        log_normal=-(1/2)*np.log(np.pi*2*var)-(1/2)*((v-mu)**2)/var
        # add to total
        log_likelihood=log_likelihood+log_normal
    return log_likelihood

