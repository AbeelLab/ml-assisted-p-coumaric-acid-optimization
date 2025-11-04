"""Constructs the numeric matrix from pro-orf-ter matrices"""

import pandas as pd
import re
import numpy as np


filename1 = "data/processed/IntegratedData_WURTUD/WUR_TUD_count_matrix.csv"

strain_count_matrix = pd.read_csv(filename1, index_col=0)
strain_count_matrix = strain_count_matrix.rename({'Agos_TEF1.pro_0003-Sc_GND2.orf': "Agos_lox_TEF1.pro-Sc_GND2.orf"},
                                                 axis=1)
print(np.shape(strain_count_matrix))

# filter out some pots that are accidently annotated by biobricks
filter_out = pd.read_csv("data/raw/POTS_tofilter_for_numeric.csv", index_col=0)
filter_out = pd.Series(filter_out['0']).to_list()

# ugly but its 6pm and I cannot think of anything better
for name in filter_out:
    if name in strain_count_matrix.columns.to_list():
        strain_count_matrix = strain_count_matrix.drop(name, axis=1)

print(np.shape(strain_count_matrix))
strain_count_matrix=strain_count_matrix.drop(labels='Unnamed: 1',axis=1)

# screening_results=pd.read_csv("data/merged_dataset_batchcorrected.csv",index_col=0)

promoter_strengths = pd.read_csv("results/GFP_PromoterStrengths/quantified_promoter_strength_rewritten.csv", index_col=0)
promoter_strengths = promoter_strengths['FI'].to_dict()


# make a numeric vector for each column in strain_count_matrix that has the proper numeric value attached to that particular
# promoter.

def get_numeric_matrix(strain_count_matrix, promoter_strengths):
    """include promoter strength in the count matrix"""
    pro_orfs = strain_count_matrix.columns.to_list()
    colonies = strain_count_matrix.index.to_list()
    numeric_values_pro = []
    for i in pro_orfs:
        promoter = i.split("-")[0]

        promoter = re.sub('.pro*', "", promoter)
        promoter = re.sub('.Pro*', "", promoter)
        promoter = re.sub('_0001', "", promoter)

        # some exception due to mismatch in biobrick annotation
        # if promoter=="Agos_lox_TEF1":
        #     promoter="Agos_lox_TEF1.pro"

        numeric_values_pro.append(promoter_strengths[promoter])

    # reverse log scale for summing (10**log(A)+10**log(B))
    numeric_values_pro = 10 ** np.array(numeric_values_pro)
    strain_numeric_matrix = np.array(strain_count_matrix) * numeric_values_pro
    strain_numeric_matrix = pd.DataFrame(strain_numeric_matrix, index=colonies, columns=pro_orfs)

    last_parts = [col.split('-')[-1] for col in strain_numeric_matrix.columns]
    print(last_parts)

    #some formatting steps that are necessary for proper grouping (was not necessary before)
    # last_parts=[col.replace(".orf_0001",".orf") for col in last_parts]
    # last_parts = [col.replace(".orf_0002", ".orf") for col in last_parts]
    # print(last_parts)

    strain_numeric_matrix = strain_numeric_matrix.groupby(last_parts, axis=1).sum()
    genes = strain_numeric_matrix.columns.to_list()
    # make dense
    strain_numeric_matrix = np.array(strain_numeric_matrix)
    strain_numeric_matrix = np.log10(strain_numeric_matrix, out=np.zeros_like(strain_numeric_matrix),
                                     where=(strain_numeric_matrix != 0))
    strain_numeric_matrix = pd.DataFrame(strain_numeric_matrix, index=colonies, columns=genes)
    return strain_numeric_matrix


strain_numeric_matrix = get_numeric_matrix(strain_count_matrix, promoter_strengths)

strain_numeric_matrix=strain_numeric_matrix.drop(labels=['Shyg_KanMX.orf_0003'],axis=1)
strain_numeric_matrix.to_csv("data/processed/IntegratedData_WURTUD/strain_numeric_matrix.csv")
