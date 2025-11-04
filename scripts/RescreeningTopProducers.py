"""Analysis of rescreened values of the top X producers"""

import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt
from scipy.stats import linregress, spearmanr

directory = "data/raw/CycleTUD/Rescreening_Top/"
sample_list = pd.read_excel(directory + "SampleList_AI4Bio Top86 rescreen_Feb2024.xlsx")
screening = pd.read_excel(directory + "Overview_f31982-f31985.xlsx", header=4)

# screening=screening.iloc[:,0:6]
screening = screening.iloc[:, 0:6]
positions = screening['Well'].values
positions = [s[0] + str(int(s[1:])) for s in positions]
screening['Well'] = positions

screening_dataset = pd.merge(sample_list, screening, how="left", left_on=['ANA plate', 'Well'],
                             right_on=["Title", "Well"])

#get controls
controls = ["SHK0058", "SHK0066"]

#only take the controls that are similar
controls = screening_dataset[screening_dataset['Original SampleName'].isin(controls)]

transformation_factors = controls.groupby(by="ANA plate")['Coumaric acid'].median()
screening_dataset = screening_dataset.join(transformation_factors, on="ANA plate", how="left", rsuffix="_normalization")
screening_dataset['batch_corrected_Coumaric_acid'] = screening_dataset['Coumaric acid'] / screening_dataset[
    'Coumaric acid_normalization']

# this doesnt look amazing for plate correction. Might want to discuss with Joep an Rianne
controls[controls['Original SampleName'] == "SHK0066"]['Coumaric acid'].plot(kind="kde", label="SHK0066")
# controls[controls['Original SampleName']=="SHK0058"]['Coumaric acid'].plot(kind="kde",label='SHK0058')
# controls[controls['Original SampleName']=="SHK0046"]['Coumaric acid'].plot(kind="kde",label='SHK0046')
plt.legend()
plt.show()
#
sorted = controls.sort_values(by="Coumaric acid")

from matplotlib.colors import ListedColormap

color_dict = dict(zip(sorted['ANA plate'].unique(), ['#00468BFF', '#ED0000FF', '#42B540FF', '#925E9FFF']))
# Create a list of numeric values corresponding to each category
categories = sorted['ANA plate']
category_numeric = [list(color_dict.keys()).index(cat) for cat in categories]
# Create a colormap based on the color dictionary
cmap = ListedColormap([color_dict[cat] for cat in color_dict])
# Plot the scatter plot with specified colors
plt.scatter(np.arange(0, len(sorted['Coumaric acid'])), sorted['Coumaric acid'], c=category_numeric, cmap=cmap)
plt.title("sorted control strain production")
# Show the plot
plt.colorbar(label='ANA plate')  # Add colorbar for reference
plt.ylabel("[p-CA]")
plt.show()


previous_measurements=pd.read_excel(directory+"Plate lay-out AI4Bio sequencing Coumaric acid hits 14Feb2024.xlsx")
# previous_measurements=previous_measurements[previous_measurements['ANA plate']!="f31883"]
previous_measurements=previous_measurements.filter(items=['SampleName','batch_corrected_Coumaric_acid'])

#remove this plate 24



newdataset=screening_dataset.filter(['Original SampleName','batch_corrected_Coumaric_acid','ANA plate'])
previous_measurements=previous_measurements.rename({"batch_corrected_Coumaric_acid":"m1_pCA"},axis=1)
newdataset=newdataset.rename({'batch_corrected_Coumaric_acid':'m2_pCA'},axis=1)

## Do some filtering of measurements that failed
final_ds=pd.merge(newdataset,previous_measurements,how="left",left_on="Original SampleName",right_on='SampleName')
final_ds=final_ds.dropna()
# final_ds=final_ds.sort_values('m1_pCA')

final_ds=final_ds[final_ds['m2_pCA']>0]

corrected_remeasured_screening_results_top=final_ds[['SampleName','m2_pCA']]

corrected_remeasured_screening_results_top=corrected_remeasured_screening_results_top.rename({"m2_pCA":"batch_corrected_Coumaric_acid"},axis=1)
corrected_remeasured_screening_results_top.to_csv("data/processed/"
                                                  "IntegratedData_WURTUD/"
                                                  "2504_topX_corrected_screening_results.csv")






plate_test=final_ds.sort_values(['m2_pCA'])

color_dict = dict(zip(plate_test['ANA plate'].unique(), ['#00468BFF', '#ED0000FF', '#42B540FF', '#925E9FFF']))
# Create a list of numeric values corresponding to each category
categories = plate_test['ANA plate']
category_numeric = [list(color_dict.keys()).index(cat) for cat in categories]
# Create a colormap based on the color dictionary
cmap = ListedColormap([color_dict[cat] for cat in color_dict])
# Plot the scatter plot with specified colors
plt.scatter(np.arange(np.shape(plate_test)[0]),plate_test['m2_pCA'],c=category_numeric, cmap=cmap)
plt.ylabel("p-CA")
plt.title("ANA Plate coloring per strain measurement")
plt.show()



mean_m2=final_ds.groupby(['Original SampleName'])['m2_pCA'].mean()
std_m2=final_ds.groupby(['Original SampleName'])['m2_pCA'].std()
ANA_plate=final_ds.groupby(['Original SampleName'])['m2_pCA'].std()

measurement_2=pd.concat([mean_m2,std_m2],axis=1)
measurement_2.columns=['mean','std']
measurement_2=measurement_2.sort_values('mean')

plt.bar(np.arange(0,np.shape(measurement_2)[0]),measurement_2['mean'],yerr=measurement_2['std'])
plt.xlabel("strains (sorted)")
plt.ylabel("plate corrected [p-CA]")
plt.show()

#%%

final_ds = final_ds.groupby(by="SampleName").median()

fig, ax = plt.subplots(figsize=(4, 4))


ax.scatter(final_ds['m1_pCA'], final_ds['m2_pCA'], c="#00468BFF")
ax.set_ylabel("New Measurement (4 replicates)")
ax.set_xlabel("Old Measurement (1 replicate)")
ax.set_title("Rescreened top 86 strains")

fig.tight_layout()
plt.show()

fig.savefig("figures/library_visualization/rescreened_top86.png",bbox_inches='tight')
fig.savefig("figures/library_visualization/rescreened_top86.svg",bbox_inches='tight')



# res= linregress(final_ds['m2_pCA'],final_ds['m1_pCA'])