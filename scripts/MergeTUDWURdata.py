"""This scripts integrates data from previous rounds together"""

import pandas as pd
import numpy as np
import os



TUD_count_matrix_filename="data/processed/CycleTUD/strain_count_matrix_tud.csv"

WUR_count_matrix_filename=("data/raw/CycleWUR/"
                           "Copy of Shiki-factory-15-flowcells-_assembly-based_cp_nr-id97.5-cov95_shk0046.xlsx")
WUR_count_matrix=pd.read_excel(WUR_count_matrix_filename,index_col=0).T
TUD_count_matrix=pd.read_csv(TUD_count_matrix_filename,index_col=0)


#processing steps in WUR count matrix
POT_names=[]
for i in WUR_count_matrix.columns.to_list():
    POT_names.append(i.replace("__","-")) #some change in structure
WUR_count_matrix.columns=POT_names

PAL_colonies=[]
for i in WUR_count_matrix.index.to_list():
    if "PAL" in i: #only get PAL_colonies
        PAL_colonies.append(i)

WUR_count_matrix=WUR_count_matrix.filter(PAL_colonies,axis=0)
print(np.shape(WUR_count_matrix))
WUR_count_matrix=WUR_count_matrix.iloc[:,list(np.where(WUR_count_matrix.sum(axis=0)!=0)[0])]
# WUR_count_matrix
WUR_TUD_count_matrix=pd.concat([TUD_count_matrix,WUR_count_matrix])
print(np.shape(WUR_TUD_count_matrix))

WUR_TUD_count_matrix=WUR_TUD_count_matrix.fillna(value=0)



# ## this count matrix can now be saved and converted to a dense numeric matrix
WUR_TUD_count_matrix.to_csv("data/processed/IntegratedData_WURTUD/WUR_batch_corrected_screening.csv")

##
WUR_screening_filename="data/processed/CycleWUR/WUR_batch_corrected_screening.csv"
TUD_screening_filename="data/processed/CycleTUD/TUD_batch_corrected_screening.csv"

##
WUR_screening_results=pd.read_csv(WUR_screening_filename,index_col=0)['batch_corrected_Coumaric_acid']
TUD_screening_results=pd.read_csv(TUD_screening_filename,index_col=1)['batch_corrected_Coumaric_acid']

WUR_TUD_screening_results=pd.concat([TUD_screening_results,WUR_screening_results])
WUR_TUD_screening_results.to_csv("data/processed/IntegratedData_WURTUD/"
                                 "batch_corrected_screening_results_integrated.csv")