
import numpy as np
import itertools
import pandas as pd
import xgboost as xgb
from xgboost import XGBClassifier,XGBRegressor
import time
import matplotlib.pyplot as plt
import os
import sys




def get_possibilities(unique_values,dim,pos,ref_strain):
    """Returns the list of possible combinations for the algorithm
    Input: 
    unique_values: number of unique values to test for a certain gene
    dim: number of dimensions
    pos: number of genes to consider for optimization. 
    """
    ref_strain_instances=np.where(ref_strain!=0)[0]
    # print(ref_strain_instances)

    features=np.arange(0,dim,1)
    combinations=list(itertools.product(features,unique_values))

    #this makes sure that the parent strain genes are having other perturbation values
    for m,i in enumerate(combinations):
        for k in ref_strain_instances:
            if i[0]==k:
                # print(i,k)
                i=(k,np.log10(10**i[1]+10**ref_strain[k]))
                # print(i,k)
        combinations[m]=i
    
    combinations=list(itertools.combinations(combinations,pos))

    return combinations

# i=(k,np.log10(10**i[1]+10**ref_strain[k]))
def predict_combinations(model,combinations,pos,ref_strain):
    """Predict the combinatorial possibilities.
      Number of positions is retrieved from combinations list itself, and should be similar to pos in get_possibilities"""

    data_matrix=np.zeros((len(combinations),len(ref_strain)))

    ref_strain_instances=np.where(ref_strain!=0)[0]
    for i in ref_strain_instances:
        fill_data_with_ref=[ref_strain[i]]*len(combinations)
        data_matrix[:,i]=fill_data_with_ref

    #we need to change some aspects to ensure 

    for k,combi in enumerate(combinations):
        for m in range(pos):
            data_matrix[k,combi[m][0]]=combi[m][1]

    #     for k in range(pos):
    #         vector[combi[k][0]]=combi[k][1]
        
    x=xgb.DMatrix(data_matrix)
    prediction=model.predict(x)


    return prediction

# def change_combinations_for_reference_strain(ref_strain,combinations,unique_values):
#     """Given the reference strains, change the promoter strengths according to what is already in there"""
#     parent_strain_instances=np.where(ref_strain!=0)
#     for i in parent_strain_instances:
#         for k in combinations:
#             for m in k: