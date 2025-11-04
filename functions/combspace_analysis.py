import numpy as np
import pandas as pd
from scipy.integrate import simpson
from scipy.stats import entropy,dirichlet




#remove anything that is not included more than 2 times. We have a nested dictionary structure
def filter_states(integration_sites,threshold):
    """Due to many double integrations, there are some values that only are included once. This function filters these out given a cutoff"""
    filtered_integration_sites={}
    for key,promoters in integration_sites.items():
        temp_dict={}
        for key2,values in promoters.items():
            if values>= threshold: 
                temp_dict[key2]=values
        filtered_integration_sites[key]=temp_dict
    return filtered_integration_sites



def sample_combinatorial_designs(filtered_integration_sites,N):
    """Due to impossibility of brute-forcing, we sample a large set that has each possible promoter choice equally represented.
    Input: filtered integration sites: consists of a nested dictionary where each key is a gene and the values are the promoter strengths"""
    dataframe=np.zeros((N,len(list(filtered_integration_sites.keys()))))
    max_differences=0
    for k,i in enumerate(filtered_integration_sites.keys()):
        choices=filtered_integration_sites[i]

        choices=np.array(list(choices.keys()))
        sampling=np.random.choice(choices,N)

        dataframe[:,k]=sampling

    dataframe=pd.DataFrame(dataframe,columns=filtered_integration_sites.keys())
    return dataframe,max_differences


def generate_frequency_Gene(X_i,pred_y,threshold):
    X_i=np.array(X_i)
    values=np.where(pred_y>threshold)[0]
    subset=X_i[values]
    promoter, counts=np.unique(subset,return_counts=True)
    freq=counts/len(subset)
    freq_dict=dict(zip(promoter,freq))
    return freq_dict

def get_remaining_designs(X_i,pred_y,threshold):
    X_i=np.array(X_i)
    values=np.where(pred_y>threshold)[0]
    subset=X_i[values]
    return len(subset)



def scan_combinatorial_space(combinatorial_design_space,init, n_evals=400):
    """"Scans over the combinatorial space and calculates the frequency matrix. 
    INPUT:
    combinatorial_design_space: consists of all combinatorial designs (or a sampling). The last column is the predicted value
    init: part from where to start the thresholding. Mostly useful when only interested in producers above a certain value
    n_evals: number of evaluations on the interval [0-max production]"""
    max_production=np.max(combinatorial_design_space.iloc[:,-1])-0.00001
    threshold_range=np.linspace(0,max_production,n_evals)

    enz_names=list(combinatorial_design_space.columns[:-1])
    pred_y_name=combinatorial_design_space.columns[-1]

    combinatorial_design_space=combinatorial_design_space.sort_values(by=pred_y_name)
    combinatorial_design_space=np.array(combinatorial_design_space)

    #this part of the function maps for each feature (gene) generate_frequency_Gene over the threshold_range
    X=combinatorial_design_space[:,:-1]
    pred_y=combinatorial_design_space[:,-1]

    nrows,ncols=np.shape(X)
    frequency_dict={}

    #do this ones to retrieve remaining designs given thresholding
    
    remaining_designs=list(map(lambda t: get_remaining_designs(X[:,0],pred_y,t),threshold_range))
    
    for i in range(ncols):
        print(i)
        frequency_threshold_matrix=list(map(lambda t: generate_frequency_Gene(X[:,i],pred_y,t),threshold_range))
        frequency_dict[enz_names[i]]=frequency_threshold_matrix
    return frequency_dict,enz_names,threshold_range,remaining_designs


def get_dirichlet_deviation(states,N_samples,int_sigma):
    """Deviation from baseline when each promoter is equilly represented
    int_sigma: number of standard deviations from the baseline"""
    alpha=(np.ones(states)/states)*(N_samples)

    var=dirichlet.var(alpha)[0]
    deviation=1/states+(np.sqrt(var)*int_sigma)
    return deviation
    
    
#get combinations
def construct_paired_combinatorial_space(combinatorial_design_space):
    """Makes a paired feature combinatorial design space"""
    shape=len(combinatorial_design_space.columns)-1
    a=np.ones((shape,shape))
    idx=np.triu_indices(a.shape[0],1)
    combis=[(idx[0][i],idx[1][i]) for i in range(len(idx[0]))]


    my_combination_dict={}
    for i in combis:
        row_i=combinatorial_design_space.iloc[:,i[0]].values
        row_j=combinatorial_design_space.iloc[:,i[1]].values
        combination=list(zip(row_i,row_j))
        string=combinatorial_design_space.columns[i[0]]+"_"+combinatorial_design_space.columns[i[1]]
        my_combination_dict[string]=combination
        


    c2_combinatorial_design_space=pd.DataFrame(my_combination_dict)
    c2_combinatorial_design_space['p-CA']=combinatorial_design_space['p-CA']
    return c2_combinatorial_design_space



