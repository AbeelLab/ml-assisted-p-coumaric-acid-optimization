"""This script processes the GFP-fluorescence to a promoter strength.csv file. Also plots the mean promoter strengths"""

import pandas as pd
import numpy as np
import regex as re
import FlowCal
import os
from matplotlib import pyplot as plt



def promoter_strength_dictionary(files):
    """Get the strengths of the promoters for the biological replicates, return a dictionary"""
    keys = []
    values = []
    for i in range(len(files)):

        splitted = re.split("_", files[i])  # split files based on _
        for i in splitted:
            if re.search("-", i) is not None:  # if there is an element with a stripe, this is the promoter-term pair
                i = i.replace("(2)", "")  # replace with 2
                keys.append(i)
                mean_fluorescence, median_fluorescence = get_mean_fluorescence(i,
                                                                               files)  # promoter i, get the mean fluorescence
                values.append(mean_fluorescence)
    promoter_strengths = dict(zip(keys, values))
    promoter_strengths = pd.Series(promoter_strengths)
    return promoter_strengths


def get_mean_fluorescence(promoter_name, files):

    files_to_consider = []
    for k, i in enumerate(files):
        if re.search(promoter_name, i):
            files_to_consider.append(i)
    mean_list = []
    median_list = []
    # fig,ax=plt.subplots(figsize=(4,4))
    for k, i in enumerate(files_to_consider):
        p = FlowCal.io.FCSData("data/raw/PromoterScreen_FCS_Elif/" + str(i))
        p_trans = FlowCal.transform.to_rfi(p, channels='GFP-A')

        p = np.mean(p_trans[:, 'GFP-A'])
        # p=np.log10(p)
        mean_list.append(np.mean(p))
        median_list.append(np.median(p))

    mean_fluorescence = np.log10(np.mean(mean_list))
    median_fluorescence = np.log10(np.mean(median_list))
    return mean_fluorescence, median_fluorescence


files = os.listdir("data/raw/PromoterScreen_FCS_Elif/")
promoter_strengths = promoter_strength_dictionary(files)

promoter_strengths=pd.DataFrame(promoter_strengths.sort_values(),columns=["FI"])

promoter_strengths.to_csv("results/GFP_PromoterStrengths/quantified_promoter_strength.csv")



fig,ax=plt.subplots(figsize=(4,6))
ax.barh(promoter_strengths.index.to_list(), promoter_strengths["FI"],color="#00468BFF")
ax.set_xlabel("Log GFP-fluorescence",fontsize=16)
ax.set_ylabel("Promoter-terminator pair",fontsize = 16)
ax.set_title("Promoter strengths",fontsize=20)
# plt.tight_layout()
plt.show()
fig.savefig("figures/GFP_promoter_strength/GFP_promoter_strength.svg",bbox_inches='tight')
fig.savefig("figures/GFP_promoter_strength/GFP_promoter_strength.png",bbox_inches='tight')
