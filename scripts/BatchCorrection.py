"""Batch correction and some visualization of the screening results.
Notebook  221223_screening_results.ipynb for more plots
"""

# In[1]:


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt



import seaborn as sns



# ## Loading data and preprocessing
# Preprocessing steps:
# 1. Remove all failed measurements (-1 values)
# 2. Remove all strains with glucose>1 (these did not grow)
# 3. Remove batch plate effects by using controls per plate as reference (SHK0058 and SHK0066)
#

# In[2]:


screening=pd.read_excel("data/raw/CycleTUD/Screening/Overview_f31860-f31891.xlsx",header=4)

#remove 0 from A01, A02--> A1,A2
screening=screening.iloc[:,0:6]
positions=screening['Well'].values
positions = [s[0] + str(int(s[1:])) for s in positions]
screening['Well']=positions

sample_list=pd.read_excel("data/raw/CycleTUD/Screening/SampleList_ANA_CoumlibrariesForARF.xlsx")


screening_dataset=pd.merge(sample_list,screening,how="left",left_on=['ANA plate','position'],right_on=["Title","Well"])
screening_dataset['Design Id']=screening_dataset['Design Id'].astype(str)
# screening_dataset.to_csv("merged_screening.csv")
# np.shape(screening_dataset)
screening_dataset.drop(list(np.where(screening_dataset['Coumaric acid']==-1)[0]),axis=0,inplace=True)
screening_dataset=screening_dataset[screening_dataset['Glucose']<=1]
controls=screening_dataset[screening_dataset['Design Id']=="nan"]




# In[3]: an uncorrected pCA screening plot

#plot of all samples production values
fig,ax=plt.subplots(figsize=(4,4))
colors=screening_dataset['Trehalose']
p=ax.scatter(screening_dataset['Glucose'],screening_dataset['Coumaric acid'],c=colors,cmap="Reds")
ax.set_xlabel("[Glucose]")
ax.set_ylabel("[p-CA]")
ax.set_ylim(0,2)
ax.set_title("p-CA production (all samples)")
cbar = fig.colorbar(p, ax=ax)
cbar.set_label('[Trehalose]')
fig.tight_layout()
fig.savefig("figures/library_visualization/uncorrected_screening_TUD.png")
fig.savefig("figures/library_visualization/uncorrected_screening_TUD.svg")
fig.show()





# In[4]:


subset_library1=np.where(screening_dataset['Design Id']=="Coum_Library1")[0][-1]

library1=screening_dataset.iloc[0:subset_library1,:]
controls_l1=library1[library1['Design Id']=="nan"]

subset_library2=np.where(screening_dataset['Design Id']=="Coum_Library2")[0][-1]
library2=screening_dataset.iloc[subset_library1+1:subset_library2,:]
controls_l2=library2[library2['Design Id']=="nan"]


subset_library3=np.where(screening_dataset['Design Id']=="Coum_Library3")[0][-1]
library3=screening_dataset.iloc[subset_library2+1:subset_library3,:]
controls_l3=library3[library3['Design Id']=="nan"]


subset_library4=np.where(screening_dataset['Design Id']=="Coum_Library4")[0][-1]
library4=screening_dataset.iloc[subset_library3+1:subset_library4,:]
controls_l4=library4[library4['Design Id']=="nan"]




## Controls within library 1
control_mean_l1=controls_l1.groupby("SampleName")['Coumaric acid'].mean()
control_std_l1=controls_l1.groupby("SampleName")['Coumaric acid'].std()

control_mean_l2=controls_l2.groupby("SampleName")['Coumaric acid'].mean()
control_std_l2=controls_l2.groupby("SampleName")['Coumaric acid'].std()

control_mean_l3=controls_l3.groupby("SampleName")['Coumaric acid'].mean()
control_std_l3=controls_l3.groupby("SampleName")['Coumaric acid'].std()

control_mean_l4=controls_l4.groupby("SampleName")['Coumaric acid'].mean()
control_std_l4=controls_l4.groupby("SampleName")['Coumaric acid'].std()

means=[control_mean_l1,control_mean_l2,control_mean_l3,control_mean_l4]
std=[control_std_l1,control_std_l2,control_std_l3,control_std_l4]

## plot of the controls per library
fig,ax=plt.subplots(figsize=(6,4))
x=np.arange(len(means[0]))
ax.bar(x,height=control_mean_l1,width=0.20,yerr=control_std_l1,label="lib1")
ax.bar(x+0.2,height=control_mean_l2,width=0.2,yerr=control_std_l2,label="lib2")
ax.bar(x+0.4,height=control_mean_l3,width=0.2,yerr=control_std_l3,label="lib3")
ax.bar(x+0.6,height=control_mean_l4,width=0.2,yerr=control_std_l4,label="lib4")
ax.set_xticks(x, control_mean_l1.keys())
ax.set_ylabel("[p-CA]")
ax.set_title("Control strains per library")
ax.legend(loc="upper left",bbox_to_anchor=(1.04,1))
fig.tight_layout()
fig.savefig("figures/library_visualization/controls_per_library_design_TUD.png")
fig.savefig("figures/library_visualization/controls_per_library_design_TUD.svg")
fig.show()




# In[6]:



controls=controls[controls['SampleName']!="Blanc"]
controls=controls[controls['SampleName']!="SHK0046"]
print(np.shape(controls))
production_per_plate_controls_mean=controls.groupby(by='GEN plate name')[['Coumaric acid']].mean()
production_per_plate_controls_median=controls.groupby(by='GEN plate name')[['Coumaric acid']].median()
production_per_plate_controls_std=controls.groupby(by='GEN plate name')[['Coumaric acid']].std()



## plot of the controls per plate
mean_controls=production_per_plate_controls_mean['Coumaric acid'].values
median_controls=production_per_plate_controls_median['Coumaric acid'].values
x=np.arange(len(mean_controls))


fig,ax=plt.subplots(figsize=(12,4))
ax.bar(x-0.2,height=mean_controls,width=0.4,label="mean controls",yerr=production_per_plate_controls_std['Coumaric acid'].values)
ax.bar(x+0.2,height=median_controls,width=0.4,label="median controls")
ax.set_xticks(x, list(production_per_plate_controls_median.index.get_level_values(0)),rotation=90)
ax.set_ylabel("[p-CA]")
ax.set_title("Mean/Median control p-CA production per plate")
ax.legend(loc="upper left",bbox_to_anchor=(1.0,1))
ax.set_xlabel("Microtiter plates")
fig.tight_layout()
fig.savefig("figures/library_visualization/controls_per_MTP_TUD.png")
fig.savefig("figures/library_visualization/controls_per_MTP_TUD.svg")
fig.show()


# In[7]:


controls=screening_dataset[screening_dataset['Design Id']=="nan"]
fig,ax=plt.subplots(figsize=(6,4))
#check with a test whether SHK0058 and SHK0066 are significantly different
# sns.kdeplot(controls[controls['SampleName']=="SHK0046"]['Coumaric acid'],alpha=0.8,label="SHK0046")
sns.kdeplot(controls[controls['SampleName']=="SHK0046"]['Coumaric acid'],
            alpha=0.8,label="SHK0046",
            linewidth=3.0,c="#ED0000FF")
sns.kdeplot(controls[controls['SampleName']=="SHK0058"]['Coumaric acid'],
            alpha=0.8,label="Parent +KanMX (SHK0058)",linewidth=3.0,c="#00468BFF")
sns.kdeplot(controls[controls['SampleName']=="SHK0066"]['Coumaric acid'],
            alpha=0.8,label="Parent (SHK0066)",linewidth=3.0,c="#ADB6B6FF")

ax.axvline(0.275967677494958,linestyle="--",c="black")
ax.axvline(0.42099603783548517,linestyle="--",c="black")
ax.set_xlabel("[pCA concentration]")

ax.legend(loc="upper left",bbox_to_anchor=(1.0,1))
ax.set_title("Density for control strains")
fig.tight_layout()
fig.savefig("figures/library_visualization/controls_kde_pCA_production_TUD.png")
fig.savefig("figures/library_visualization/controls_kde_pCA_production_TUD.svg")

fig.show()

# In[8]:


### Calculate the ratio between SHK0046 and SHK0066

#remove zero counts, these likely failed
controls=controls[controls['Coumaric acid']>0.01]
med_SHK0046=np.median(controls[controls['SampleName']=="SHK0046"]['Coumaric acid'])

med_SHK0066=np.median(pd.concat([controls[controls['SampleName']=="SHK0066"],controls[controls['SampleName']=="SHK0058"]])['Coumaric acid'])

print(med_SHK0046/med_SHK0066)


# ## Batch effects
# We are going to correct batch effects between plates by taking the median strain of SHK0058 and SHK0066. We thus exclude the SHK0046 strain.

# In[20]:


#remove SHK0046
transformation_factors=controls[controls['SampleName']!="SHK0046"]
transformation_factors=transformation_factors.groupby(by="GEN plate name")['Coumaric acid'].median()
screening_dataset=screening_dataset.join(transformation_factors,on="GEN plate name",how="left",rsuffix="_normalization")
# screening_dataset['batch_cor_Coumaric acid']=screening_dataset['Coumaric acid']/screening_dataset['Coumaric acid_normalization']

# transformation_factors=controls[controls['SampleName']!="SHK0046"]
# trehalose_transform=transformation_factors.groupby(by="GEN plate name")['Trehalose'].median()
# screening_dataset=screening_dataset.join(trehalose_transform,on="GEN plate name",how="left",rsuffix="_normalization")


# In[21]: The one line that corrected batch effects

screening_dataset['batch_corrected_Coumaric_acid']=screening_dataset['Coumaric acid']/screening_dataset['Coumaric acid_normalization']
# screening_dataset['batch_corrected_treha']=screening_dataset['Trehalose']/screening_dataset['Trehalose_normalization']

screening_dataset.to_csv("data/processed/CycleTUD/TUD_batch_corrected_screening.csv")

# In[22]:


#plot of all samples production values



# Batch corrected
fig,ax=plt.subplots(figsize=(4,4))
colors=screening_dataset['Trehalose']
p=ax.scatter(screening_dataset['Glucose'],screening_dataset['batch_corrected_Coumaric_acid'],c=colors,cmap="Reds")
ax.set_xlabel("[Glucose]")
ax.set_ylabel("Plate corrected [p-CA]")
ax.set_ylim(0,4)
plt.axhline(1,linestyle="--",c="black")

ax.set_title("Plate corrected [p-CA]")
cbar = fig.colorbar(p, ax=ax)
cbar.set_label('[Trehalose]')
# plt.colorbar(label='[Trehalose]')
fig.tight_layout()
fig.savefig("figures/library_visualization/batch_corrected_screening_TUD.png")
fig.savefig("figures/library_visualization/batch_corrected_screening_TUD.svg")
fig.show()


# In[23]:


## Colored by plate
plates=screening_dataset['GEN plate name'].values
unique_values = np.unique(plates)
colors = plt.cm.tab20c(np.linspace(0, 1, len(unique_values)))
# Map each unique value to its corresponding color
color_mapping = dict(zip(unique_values, colors))
# Create a list of colors for each data point based on the mapping

color_list = [color_mapping[value] for value in plates]

# Batch corrected
plt.scatter(screening_dataset['Glucose'],screening_dataset['batch_corrected_Coumaric_acid'],c=color_list,label="plates")
plt.legend()
plt.xlabel("[Glucose]")
plt.ylabel("[p-CA]")
plt.ylim(0,4)
plt.title("p-CA production all strains (batch corrected)")
plt.show()



# In[25]:


sorted=screening_dataset[screening_dataset['Sample - sampleFunction']!="Control"]
sorted=sorted.sort_values(by="batch_corrected_Coumaric_acid",ascending=False)



pCA_sorted=sorted['batch_corrected_Coumaric_acid']
plates_sorted=sorted['GEN plate name']
libraries=sorted['Design Id']

unique_libraries=np.unique(libraries)
m=np.arange(len(unique_libraries))
unique_colors_dict_lib=dict(zip(unique_libraries,m))
coloring_lib=[unique_colors_dict_lib[i] for i in libraries]

unique_plates=np.unique(plates_sorted)
m=np.arange(len(unique_plates))
unique_colors_dict=dict(zip(unique_plates,m))
coloring=[unique_colors_dict[i] for i in plates_sorted]

fig,ax=plt.subplots(figsize=(20,4))
plt.scatter(np.arange(0,len(plates_sorted)),pCA_sorted,c=coloring_lib,s=30,cmap="Set1",label="Strains")
plt.axhline(1,c="black",linestyle="--",label="Top 188 strains")
# plt.axhline(1.81,c="black",linestyle="--",label="Control production")
# plt.axvline(1690,c="black",linestyle="--")
# plt.axvline(186,c="black",linestyle="--")
plt.xlabel("sorted strains")
plt.ylabel("batch corrected coumaric acid value")
plt.legend()
plt.title("Screening resuls")
plt.show()





# In[27]:

## density plot of the
fig,ax=plt.subplots(figsize=(4,4))
ax.hist(screening_dataset['batch_corrected_Coumaric_acid'],
        bins=40, color="#00468BFF")
ax.set_xlabel("normalized [p-CA]")
ax.set_ylabel("# of strains")
ax.axvline(np.median(screening_dataset['batch_corrected_Coumaric_acid']),linestyle="--",c="black",label="median")
ax.text(np.median(screening_dataset['batch_corrected_Coumaric_acid']) + 0.04, ax.get_ylim()[1] * 0.9, "SHK0066", color="black", verticalalignment='center')
ax.set_title("[p-CA] screening")
fig.tight_layout()
fig.savefig("figures/library_visualization/screening_histogram_TUD.png",dpi=1200)
fig.savefig("figures/library_visualization/screening_histogram_TUD.svg")
fig.show()

