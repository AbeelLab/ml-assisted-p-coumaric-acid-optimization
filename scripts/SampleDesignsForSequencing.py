"""Script that choose the designs after screening for further sequencing.
Takes the 3 MTPs for the top strains, 2 MTPs for the mid-range producers, and 1 MTP
for strains that produce less than the parent strain"""
import os
# In[1]:


from collections import Counter
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import ListedColormap
from scipy import stats

# ## Preprocessing
# 1. Remove nonmeasured strains (-1 values in CA)
# 2. Remove strains with high glucose values (<=1)
# 3. Remove batch plate effects
#

# In[2]:


screening = pd.read_excel("data/raw/CycleTUD/Screening/Overview_f31860-f31891.xlsx", header=4)
#remove 0 from A01, A02--> A1,A2
screening = screening.iloc[:, 0:6]
positions = screening['Well'].values
positions = [s[0] + str(int(s[1:])) for s in positions]
screening['Well'] = positions

sample_list = pd.read_excel("data/raw/CycleTUD/Screening/SampleList_ANA_CoumlibrariesForARF.xlsx")

screening_dataset = pd.merge(sample_list, screening, how="left", left_on=['ANA plate', 'position'],
                             right_on=["Title", "Well"])
screening_dataset['Design Id'] = screening_dataset['Design Id'].astype(str)
# screening_dataset.to_csv("merged_screening.csv")
# np.shape(screening_dataset)
screening_dataset.drop(list(np.where(screening_dataset['Coumaric acid'] == -1)[0]), axis=0, inplace=True)
screening_dataset = screening_dataset[screening_dataset['Glucose'] <= 1]
screening_dataset



# In[3]:


## Remove plate batch effects
controls = screening_dataset[screening_dataset['Design Id'] == "nan"]
controls = controls[controls['SampleName'] != "Blanc"]
transformation_factors = controls[controls['SampleName'] != "SHK0046"]
transformation_factors = transformation_factors.groupby(by="GEN plate name")['Coumaric acid'].median()
screening_dataset = screening_dataset.join(transformation_factors, on="GEN plate name", how="left",
                                           rsuffix="_normalization")
screening_dataset['batch_corrected_Coumaric_acid'] = screening_dataset['Coumaric acid'] / screening_dataset[
    'Coumaric acid_normalization']

#remove controls and sort
sorted_screening_dataset = screening_dataset[screening_dataset['Sample - sampleFunction'] != "Control"]
sorted_screening_dataset = sorted_screening_dataset.sort_values(by="batch_corrected_Coumaric_acid", ascending=False)

# In[4]:


pCA_sorted = sorted_screening_dataset['batch_corrected_Coumaric_acid']
plates_sorted = sorted_screening_dataset['GEN plate name']
libraries = sorted_screening_dataset['Design Id']

unique_libraries = np.unique(libraries)
m = np.arange(len(unique_libraries))
unique_colors_dict_lib = dict(zip(unique_libraries, m))
coloring_lib = [unique_colors_dict_lib[i] for i in libraries]

unique_plates = np.unique(plates_sorted)
m = np.arange(len(unique_plates))
unique_colors_dict = dict(zip(unique_plates, m))
coloring = [unique_colors_dict[i] for i in plates_sorted]

fig, ax = plt.subplots(figsize=(20, 4))
plt.scatter(np.arange(0, len(plates_sorted)), pCA_sorted, s=30, cmap="Set1", label="Strains")
plt.axhline(1, c="black", linestyle="--", label="Control production")
plt.axhline(1.654523, c="black", linestyle="--", label="Top 285 strains")
plt.axvline(1690, c="black", linestyle="--")
plt.axvline(285, c="black", linestyle="--")
plt.xlabel("sorted strains")
plt.ylabel("batch corrected coumaric acid value")
plt.legend()
plt.title("Screening resuls (no controls)")
plt.text(x=30, y=3, s="3*96-wells plates\n +10 that were left")
plt.text(x=900, y=3, s="2*96-wells plates")
plt.text(x=1800, y=3, s="1*96-wells plate")
plt.show()


# In[6]:


#sampling the top 285 designs (3*96-3 blanco)
topdesigns_3_plates = sorted_screening_dataset.iloc[0:295]

#sampling 190 designs (2*96-2 blanco) from everything between 1 and 1.65
# - Bins based on p-coumaric acid value (take a bin size that 190 is divisible for simplicity)
# - Group by library and bin label,
temp = sorted_screening_dataset.iloc[295:, ]
temp = temp[temp['batch_corrected_Coumaric_acid'] >= 1]
range_to_check = np.linspace(np.min(temp['batch_corrected_Coumaric_acid']),
                             np.max(temp['batch_corrected_Coumaric_acid']), 6)
bin_labels = np.zeros(len(temp['batch_corrected_Coumaric_acid']))
for i in range(len(range_to_check) - 1):
    x = np.where(temp['batch_corrected_Coumaric_acid'] >= range_to_check[i])
    bin_labels[list(x)] = i + 1
temp['bin_labels'] = bin_labels
middesigns_2_plates = temp.groupby(by=["Design Id", "bin_labels"]).sample(9,
                                                                          random_state=1)  #180 strains now (we can add 10 to the top performers)

#For the worse strains we sample 95 strains (1*96 wells plate -1)
temp = sorted_screening_dataset[sorted_screening_dataset['batch_corrected_Coumaric_acid'] <= 1.0]
baddesigns_1_plate = temp.groupby(by="Design Id").sample(23, random_state=1)

#After this I still have 10 left from middle designs and 3 from the bad performers
# 13 strains: 1 for the parent strain.
# 12 strains


sampled_designs = pd.concat([topdesigns_3_plates, middesigns_2_plates, baddesigns_1_plate])
sampled_designs = sampled_designs.drop("bin_labels", axis=1)

sampled = np.zeros(len(sorted_screening_dataset['SampleName']))
a, ind1, ind2 = np.intersect1d(sorted_screening_dataset['SampleName'], sampled_designs['SampleName'],
                               return_indices=True)
sampled[list(ind1)] = 1
sorted_screening_dataset['sampled'] = sampled



# In[21]:


pCA_sorted = sorted_screening_dataset['batch_corrected_Coumaric_acid']
plates_sorted = sorted_screening_dataset['GEN plate name']
sampled = sorted_screening_dataset['sampled']

unique_libraries = np.unique(libraries)
m = np.arange(len(unique_libraries))
unique_colors_dict_lib = dict(zip(unique_libraries, m))
coloring_lib = [unique_colors_dict_lib[i] for i in libraries]

# unique_plates=np.unique(plates_sorted)
# m=np.arange(len(unique_plates))
# unique_colors_dict=dict(zip(unique_plates,m))
# coloring=[unique_colors_dict[i] for i in plates_sorted]
cMap = ListedColormap(['#00468BFF', '#ED0000FF'])  #,'#42B540FF','#925E9FFF'])
fig, ax = plt.subplots(figsize=(20, 4))
plt.scatter(np.arange(0, len(plates_sorted)), pCA_sorted,
            c=sampled,
            s=30,
            cmap=cMap, label="Strains", alpha=0.9)
plt.axhline(1, c="black", linestyle="--", label="Control production")
plt.axhline(1.654523, c="black", linestyle="--", label="Top 285 strains")
plt.axvline(1690, c="black", linestyle="--")
plt.axvline(285, c="black", linestyle="--")
plt.xlabel("sorted strains")
plt.ylabel("batch corrected coumaric acid value")
plt.legend()
plt.colorbar()
# plt.xlim(500,600)
# plt.ylim(1,1.6)
plt.title("Screening resuls (no controls)")
plt.show()

# In[16]:


plt.scatter(np.arange(0, len(sampled_designs['batch_corrected_Coumaric_acid'])),
            sampled_designs['batch_corrected_Coumaric_acid'].sort_values(), s=1)
plt.xlabel("Strains")
plt.ylabel("[batch-normalized p-CA]")
plt.show()

plt.scatter(sampled_designs['Glucose'], sampled_designs['batch_corrected_Coumaric_acid'].sort_values(), s=1)
# sampled_designs.to_excel("100124_designs_for_sequencing.xlsx")
plt.show()

# In[152]:


range_to_check = np.linspace(np.min(temp['batch_corrected_Coumaric_acid']),
                             np.max(temp['batch_corrected_Coumaric_acid']), 100)
cumulative = []
for i in range_to_check:
    # print(i)
    frequency = np.sum(temp['batch_corrected_Coumaric_acid'] < i) / len(temp['batch_corrected_Coumaric_acid'])
    cumulative.append(frequency)

plt.scatter(range_to_check, cumulative)
# plt.ylim(0,1)
plt.xlabel("[p-CA] value")
plt.ylabel("Frequency")
plt.show()

top5percent = temp
top5percent_libraries = Counter(top5percent['Design Id'])
keys = ['Coum_Library1', 'Coum_Library2', 'Coum_Library3', 'Coum_Library4']
top5percent_libraries = {i: top5percent_libraries[i] for i in keys}

plt.bar(top5percent_libraries.keys(), top5percent_libraries.values())

plt.ylabel("# of strains")
plt.show()


#%%

## density plot of the

pCA_sorted = sorted_screening_dataset['batch_corrected_Coumaric_acid']
sampled = sorted_screening_dataset['sampled']

cMap = ListedColormap(['#00468BFF', '#ED0000FF'])  #,'#42B540FF','#925E9FFF'])


cdf = np.arange(1, len(pCA_sorted) + 1) #/ len(pCA_sorted)

# Plot the CDF
fig, ax =plt.subplots(figsize=(8, 4))
ax.scatter(cdf, pCA_sorted,
         c=list(sampled.values),
         cmap=cMap,
         label="Empirical CDF")
ax.set_ylabel("Normalized [p-CA]")
ax.set_xlabel("Strains")
ax.set_xlim(-50,len(pCA_sorted)+50)
plt.title("Strain selection for sequencing")

# Define zoom-in range (adjust these values)
x_min, x_max = 0.23 *len(pCA_sorted), 0.24* len(pCA_sorted)  # Adjust for your dataset
y_min, y_max = 1.35, 1.45

# Inset scatter plot (zoomed-in region)
axins = inset_axes(ax, width="40%", height="50%", loc="upper right")  # Adjust size and position
axins.scatter(cdf, pCA_sorted,
              c=sampled, cmap=cMap, s=30, label="Sampled")
axins.set_xlim(x_min, x_max)
axins.set_ylim(y_min, y_max)

# Remove inset tick labels for cleaner look
axins.set_xticklabels([])
axins.set_yticklabels([])

# Add a rectangle to indicate the zoomed-in region
ax.indicate_inset_zoom(axins, edgecolor="black")

box_x = [x_min, x_max, x_max, x_min, x_min]
box_y = [y_min, y_min, y_max, y_max, y_min]


# axins_pos = axins.get_position()
ax.plot((x_max, 1500), (y_max, 3.5), linewidth=1,linestyle="--", color="black")
ax.plot((x_min, 1500), (y_min, 1.6), linewidth=1,linestyle="--", color="black")
ax.axhline(1, c="black", linestyle="--")
ax.text(x=0, y=1, s="SHK0066", fontsize=10, color='black', ha='left', va='bottom')

# plt.tight_layout()
plt.show()
fig.savefig("figures/library_visualization/selected_strains_for_sequencing_cdf.png",dpi=300)
fig.savefig("figures/library_visualization/selected_strains_for_sequencing_cdf.svg")


#%%
sorted_screening_dataset