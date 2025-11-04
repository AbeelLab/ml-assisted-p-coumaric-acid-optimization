"""Batch correction of wur screennig results
based on SHK0046"""

import pandas as pd
import numpy as np



screening_filename="data/raw/CycleWUR/Screening of Shiki Project_complete data file.xlsx"
WUR_count_matrix_filename=("data/raw/CycleWUR/"
                           "Copy of Shiki-factory-15-flowcells-_assembly-based_cp_nr-id97.5-cov95_shk0046.xlsx")

screening=pd.read_excel(screening_filename,header=0)

#get the transformation factors from controls
SHK0046=screening[screening['SampleName']=="SHK0046_"]
SHK0046=SHK0046[SHK0046['Coumaric Acid']!=-1]
transformation_factors=SHK0046.groupby('ANA plate')['Coumaric Acid'].mean()

#normalize against SHK0046
## we use a reference point estimate from our current data to normalize this data for proper integration. This means that we calculate [SHK0046]/[SHK0066] ratio
## and multiply the SHK0046 corrected column
screening_dataset=screening.join(transformation_factors,on="ANA plate",how="left",rsuffix="_normalization")
screening_dataset['batch_corrected_Coumaric_acid']=(screening_dataset['Coumaric Acid']/screening_dataset['Coumaric Acid_normalization'])*0.6549170264045384

PAL_screening_dataset=screening_dataset[screening_dataset['SampleName'].str.contains("PAL")]


#we now want to filter the ones where we had sequencing data
WUR_count_matrix=pd.read_excel(WUR_count_matrix_filename,index_col=0)
PAL_screening_dataset["SampleName"]=PAL_screening_dataset['SampleName'].str.replace(" ","_").values
intersection=np.intersect1d(PAL_screening_dataset['SampleName'],WUR_count_matrix.columns.to_list())

PAL_screening_dataset=PAL_screening_dataset.set_index(keys="SampleName",drop=False)
PAL_screening_dataset=PAL_screening_dataset.filter(items=intersection,axis=0)

PAL_screening_dataset.to_csv("data/processed/CycleWUR/WUR_batch_corrected_screening.csv")