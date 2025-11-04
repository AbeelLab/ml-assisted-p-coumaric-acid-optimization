"""This script contains the analysis of the design of experiment rounds
# ### Questions to be answered in this notebook
# 1. What was the effect of alpha on the distribution? Which strategy worked best?
# 2. Did this strategy outperform the best design build?
# 3. How did the gene ladder compare to the feature importance method findings?
"""

# In[3]:

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import spearmanr, linregress
import scipy

# In[4]:


#load screening data
screening_data = pd.read_excel("data/raw/CycleTUDValidation/Overview_f32451_f32457.xlsx", index_col=0, header=0,
                               sheet_name="NMR")

# we first like to inspect all the controls and their distribution per plate

#blank, control strain (wur), parent strain (tud), best strain (library), best strain (rebuild),best strain 2 (library), best strain (rebuild), ood strain (top 4 predicted)
controls_all_plates = ["SHK0046", "SHK0066", "SHK0073"]
controls_screening = screening_data[screening_data['Strain'].isin(controls_all_plates)]

# rebuilds=["SHK0068","SHK0081","SHK0070","SHK0082"]

plate_names = controls_screening['file'].unique()

#plot controls for all plates
fig, ax = plt.subplots(figsize=(4, 4))
X_axis = np.arange(len(plate_names))

pointer = [-0.15, 0, 0.15]
for i, control in enumerate(controls_all_plates):
    SHK_mean = controls_screening[controls_screening['Strain'] == control].groupby("file")['Coumaric acid'].mean()
    SHK_std = controls_screening[controls_screening['Strain'] == control].groupby("file")['Coumaric acid'].std()

    ax.bar(X_axis + pointer[i], SHK_mean, width=0.15, label=control, yerr=SHK_std)
ax.legend()
ax.set_xticks(X_axis, plate_names)
ax.set_ylabel("[p-CA] (NMR)")
fig.savefig("figures/ValidationRound/0512_allplatecontrols_screening.png")
plt.show()
### Plate f32457 (all rebuilds and validation)


rebuilds = ["SHK0068", "SHK0081", "SHK0070", "SHK0082"]  #validation
screening_data_f32457 = screening_data[screening_data['file'] == "f32457"]
rebuilds_data = screening_data_f32457[screening_data_f32457['Strain'].isin(rebuilds)]

# Group by 'Strain' and calculate mean and standard deviation
control_mean = screening_data_f32457[screening_data_f32457['Strain']=="SHK0066"]['Coumaric acid'].mean()

grouped = rebuilds_data.groupby("Strain")['Coumaric acid']
means = grouped.mean()/control_mean  # This is a pandas Series
means = means.reindex(index=rebuilds)
stds = grouped.std()/control_mean  # This is also a pandas Series
stds = stds.reindex(index=rebuilds)



# Plotting with error bars
full_names =["SHK0068\n original\n 10 genes", "SHK0081 \n rebuild \n 10 genes",
                       "SHK0070 \n original \n 20 genes", "SHK0082 \n rebuild, \n 20 genes"]
fig, ax = plt.subplots(figsize=(4, 4))
means.plot(kind="bar", yerr=stds, ax=ax, capsize=4, color="#00468BFF" ,edgecolor='black')
# Customize the plot
ax.set_ylabel('Normalized [p-CA]')
ax.axhline(1, linestyle="--", color="black")
ax.set_title('Normalized [p-CA] of best strains')
ax.set_xticks(range(len(full_names)))  # Correct tick positions
ax.set_xticklabels(full_names, rotation=0)  # Set labels with rotation

plt.tight_layout()
fig.savefig("figures/ValidationRound/0512_topbuilds.png")
fig.savefig("figures/ValidationRound/0512_topbuilds.svg")
plt.show()

# ### Plot the other strains
# plot the rest of the files

# In[5]:


#first we merge the two files on position, but now we use the REPORT file
screening_data = pd.read_excel("data/raw/CycleTUDValidation/Overview_f32451_f32457.xlsx", header=4,
                               sheet_name="REPORT").iloc[:, :-7]
screening_data = screening_data.rename(columns={"Title": "Analyse plate", "Well": "Position destination"})

sample_list = pd.read_excel("data/raw/CycleTUDValidation/SampleList_Screen20241104_AI4BioDoE.xlsx", index_col=0,
                            header=0, sheet_name="Sheet1")

screening_data = pd.merge(screening_data, sample_list, on=["Analyse plate", "Position destination"])

#some data cleaning, we first remove all -1 values from the coumaric acid, because these are failed measurements


# ### Sara normalized the data already
# I have changed the normalization from mean to median for consistency with previous analysis.
# 



# In[7]:


data = pd.read_csv("data/raw/CycleTUDValidation/"
                   "DBTL2_rescreening_results_ann_norm_median.csv", index_col=0)
annotations_to_include = data['annotation'].unique()

annotations_to_include = ['alpha=0', 'alpha=2.5', 'alpha=5', "out-of-distribution"]  #visualization of alpha experiments
#
data = data[data['annotation'].isin(annotations_to_include)]
data[data['annotation'] == "alpha=0"]['norm_pCA'].plot(kind="kde", label="$\\beta=0$")
data[data['annotation'] == "alpha=2.5"]['norm_pCA'].plot(kind="kde", label="$\\beta=2.5$")
data[data['annotation'] == "alpha=5"]['norm_pCA'].plot(kind="kde", label="$\\beta=5$")

plt.legend()
plt.show()

# Sort data
sorted_pCA = data.sort_values('norm_pCA')
indices = np.arange(len(data['norm_pCA']))

# Map each annotation to a color
unique_labels = annotations_to_include
color_map = {label: plt.cm.tab20(i / len(unique_labels)) for i, label in enumerate(unique_labels)}

# Get colors for the annotations
colors = [color_map[label] for label in data['annotation']]






# In[9]:

data = pd.read_csv("data/raw/CycleTUDValidation/"
                   "DBTL2_rescreening_results_ann_norm_median.csv", index_col=0)

# Filter the data
mean_var = data[data['Sequence info'] != "multi gene integration"]

# Compute stats
grouped_by_experiment_mean = mean_var.groupby("annotation")['norm_pCA'].mean()
grouped_by_experiment_std = mean_var.groupby("annotation")['norm_pCA'].var()

# Labels and order
xlabels = ['$\\beta=0 (parent)$', '$\\beta=2.5 (parent)$',
           '$\\beta=5$ (parent)', '$\\beta=2.5$ (top strain)']
keys = ['alpha=0', 'alpha=2.5', 'alpha=5', 'out-of-distribution']

# Create subplots with correct layout
fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, figsize=(10, 4))

# Plot means
ax1.set_title("Mean [p-CA] production")
ax1.bar(range(len(keys)),
        grouped_by_experiment_mean[keys],
        color="#00468BFF")
ax1.set_xticks(range(len(keys)))
ax1.set_xticklabels(xlabels, rotation=45)
ax1.set_ylabel("Mean normalized [p-CA]")

# Plot variances
ax2.set_title("Variance [p-CA] production")
ax2.bar(range(len(keys)),
        grouped_by_experiment_std[keys],
        color="#00468BFF")
ax2.set_xticks(range(len(keys)))
ax2.set_xticklabels(xlabels, rotation=45)
ax2.set_ylabel("Variance normalized [p-CA]")

# Fix layout and save
plt.tight_layout()
plt.show()

# Save figures AFTER layout is fixed
fig.savefig("figures/ValidationRound/mean_var_production.png")
fig.savefig("figures/ValidationRound/mean_var_production.svg")


#%%
data = pd.read_csv("data/raw/CycleTUDValidation/"
                   "DBTL2_rescreening_results_ann_norm_median.csv", index_col=0)
SHK0073_median = data[data['Strain']=="SHK0073"]['norm_pCA'].median()


#%% plot of the contribution of individiual genes (perturbations)
data = pd.read_csv("data/raw/CycleTUDValidation/"
                   "DBTL2_rescreening_results_ann_norm_median.csv", index_col=0)
# Filter for rows where annotation is "individual contribution"
individual_contribution = data[data['annotation'] == "individual contribution"]
individual_contribution.to_csv("results/ValidationRound/individual_gene_contribution.csv")
# Calculate mean and std

average = pd.DataFrame(individual_contribution.groupby('design')['norm_pCA'].mean()).rename(columns={'norm_pCA': 'mean_norm_pCA'})
std = pd.DataFrame(individual_contribution.groupby('design')['norm_pCA'].std()).rename(columns={'norm_pCA': 'std_norm_pCA'})

# Combine mean and std
average = average.merge(std, on='design')

# Most frequent strings per design
most_frequent_strings = (
    individual_contribution[['P1', 'O1', 'T1', 'design']].groupby('design')
    .agg(lambda x: x.mode().iloc[0] if not x.mode().empty else None)
)

# Merge most frequent strings into average dataframe
average = average.merge(most_frequent_strings, on='design', how='left')

# Sort by mean_norm_pCA
average = average.sort_values('mean_norm_pCA', ascending=False)

# Prepare labels
P1 = average['P1'].astype(str)
O1 = average['O1'].astype(str)

P1 = [i.replace("Sc_", "") for i in P1]
O1 = [i.replace("Sc_", "") for i in O1]
O1 = [i.replace("_OFP.orf_0001", "") for i in O1]
O1 = [i.replace("_0002", "") for i in O1]
P1 = pd.Series([i.replace(".pro", "") for i in P1])
O1 = pd.Series([i.replace(".orf", "") for i in O1])

names = P1 + "-" + O1

# Plot with error bars
fig, ax1 = plt.subplots(figsize=(8, 4),)
ax1.bar(names, average['mean_norm_pCA'], yerr=average['std_norm_pCA'], color="#00468BFF", capsize=3)
ax1.set_ylim(0,2)
ax1.set_xticklabels(names, rotation=45)
ax1.set_ylabel("Normalized [p-CA]")
ax1.set_xlabel("Library element")
ax1.axhline(1.01193953456951,linestyle="--", color="black")
# plt.tight_layout()
plt.show()
fig.savefig("figures/ValidationRound/Individual_contribution_production.png", bbox_inches="tight")
fig.savefig("figures/ValidationRound/Individual_contribution_production.svg", bbox_inches="tight")

#%%
from scipy.stats import ttest_ind
controls = data[data['annotation'] =="plate_control_SHK0066"]

results =[]

return_significant_genes = []
for o1_group, group_df in individual_contribution.groupby('O1'):
    test_group = group_df['norm_pCA'].dropna()
    # Perform independent two-sample t-test (assuming unequal variance)

    t_stat, p_val = ttest_ind(test_group,controls['norm_pCA'])
    if t_stat > 0 and p_val < 0.05:
        return_significant_genes.append(o1_group)

print(return_significant_genes)
#%%
# It has to be noted that the ARO4 (first one) has a different promoter
# This was because we did not get the other ARO4 promoter working (may be because
# of too large protein load?

data = pd.read_csv("data/raw/CycleTUDValidation/DBTL2_rescreening_results_ann_norm_median.csv", index_col=0)
data_gene_ladder = data[data['annotation'] == "Gene Ladder"]
individual_contribution = data[data['annotation'] == "individual contribution"]
aro4 = individual_contribution[individual_contribution['O1']=="Sc_PHA2.orf"]
aro4 =  aro4[aro4['P1']=="Sc_ENO2.pro"]

data_gene_ladder = pd.concat([aro4,data_gene_ladder])
#Sc_PGK1.pro_0001
# data_gene_ladder[['O1',"O2","O3","O4","O5","O6","norm_pCA"]]
data_gene_ladder_mean = data_gene_ladder.groupby(['P1', 'O1', "P2", "O2", "P3", "O3", "P4", "O4", "P5", "O5", "P6", "O6"],
                                            as_index=False).agg({"norm_pCA": "mean"})
std = data_gene_ladder.groupby(['P1', 'O1', "P2", "O2", "P3", "O3", "P4", "O4", "P5", "O5", "P6", "O6"],
                                            as_index=False).agg({"norm_pCA": "std"})



data_gene_ladder_mean.to_csv("results/ValidationRound/gene_ladder.csv")
#


data_gene_ladder_mean['number_of_genes'] = ["PHA2","+ARO4","+ARO2",
                                            "+RKI1","+YHM2","+C4H"]

data_gene_ladder_mean.loc[0,'norm_pCA'] = 1.37 #this is the average of the two promoters
# that came out in the individual contributions ##FIX

fig,ax2 = plt.subplots(figsize=(4, 4))
ax2.bar(data_gene_ladder_mean['number_of_genes'],
        data_gene_ladder_mean['norm_pCA'],
        yerr=2*std['norm_pCA'],
       color="#00468BFF",)


ax2.set_xlabel("Genes")
ax2.set_ylim(0,2)
# ax2.yaxis.set_ticks([])           # removes ticks
# ax2.set_yticklabels([])           # removes tick labels

# ax2.set_ylabel("Normalized pCA")
plt.tight_layout()
plt.axhline(1.01193953456951,linestyle="--", color="black")
plt.show()
fig.savefig("figures/ValidationRound/gene_ladder.png")
fig.savefig("figures/ValidationRound/gene_ladder.svg")
# individual_contribution = data[data['annotation'] == "individual contribution"]
# individual_contribution[individual_contribution['O1']=="Sc_ARO4_OFP.orf_0001"].groupby(["P1","O1","T1"],as_index=False).agg({"norm_pCA":"mean"})
#



# In[67]:

## we try to establish how well the model performs in the external validation set

strain_studio_data = pd.read_excel("data/raw/CycleTUDValidation/Strainstudio_AI4bioDoE.xlsx",
                                   sheet_name="strain studio export")
strain_studio_data = strain_studio_data[['name', 'Optional_column_2']]
strain_studio_data = strain_studio_data.rename({"name": "design", "Optional_column_2": "predicted_pCA"}, axis=1)
data_merged = pd.merge(data, strain_studio_data, on=['design'], how='left')

#for this analysis, we only include the truly confirmed strains
data_merged = data_merged[data_merged['Sequence info'] == "confirmed"]
# #alpha=
alpha5 = data_merged[data_merged['annotation'] == "alpha=5"]
alpha25 = data_merged[data_merged['annotation'] == "alpha=2.5"]
alpha0 = data_merged[data_merged['annotation'] == "alpha=0"]
alpha_ood = data_merged[data_merged['annotation'] == "out-of-distribution"]



# plt.violinplot([alpha0,alpha5])
fig, (ax1,ax2) = plt.subplots(nrows= 1,
                       ncols=2,
                       figsize=(12, 4),
                       gridspec_kw={'width_ratios': [2, 1]})  # ratio of widths)
violins = ax1.violinplot(dataset=[alpha0['norm_pCA'],
                                  alpha25['norm_pCA'],
                                  alpha5['norm_pCA'],
                                  alpha_ood['norm_pCA']],
               vert=False,
               showmedians=True,
               showextrema=False, )

ax1.set_yticks([1, 2, 3, 4])
ax1.set_yticklabels(["$\\beta=0$ (Parent) \n Exploration",
                        "$\\beta=2.5$ (Parent) \n Balanced",
                        "$\\beta=5$ (Parent) \n Exploitation",
                     "$\\beta=2.5$ (Top) \n Balanced",
                     ], rotation=0)
ax1.set_xlabel("Normalized [p-CA]")

# Define matching colors
colors = ["#00468BFF", "#ED0000FF","#925E9FFF","#ADB6B6FF", ]

# Apply colors to each violin body
for i, body in enumerate(violins['bodies']):
    body.set_facecolor(colors[i])
    body.set_edgecolor("black")
    body.set_alpha(1)

# Color median line black
if 'cmedians' in violins:
    violins['cmedians'].set_color("black")
    violins['cmedians'].set_linewidth(1.5)

if 'cbars' in violins:
    violins['cbars'].set_color("black")
# Color min and max lines black (extrema)
if 'cmins' in violins:
    violins['cmins'].set_color("black")
    violins['cmins'].set_linewidth(1.2)
if 'cmaxes' in violins:
    violins['cmaxes'].set_color("black")
    violins['cmaxes'].set_linewidth(1.2)

# Define jitter amount to avoid overlapping points
jitter = 0.01
# Dataset list must match the order of the violinplot
datasets = [alpha0['norm_pCA'],
            alpha25['norm_pCA'],
            alpha5['norm_pCA'],
            alpha_ood['norm_pCA']]

# Add individual data points
for i, data in enumerate(datasets):
    y = np.random.normal(loc=i + 1, scale=jitter, size=len(data))  # vertical jitter
    ax1.scatter(data, y, color='black', edgecolor='black', s=10, alpha=0.7, zorder=3)


# Add vertical lines
ax1.axvline(1.776583, color='black', linestyle='--', linewidth=1.5)
ax1.axvline(1.0, color='black', linestyle='--', linewidth=1.5)

# Add text annotations for each line
ax1.text(1.776583 + 0.02, 1, "Greedy\n(best)", rotation=0,
         verticalalignment='center', color='black', fontsize=10)

ax1.text(1.0 + 0.02, 3.5, "SHK0066", rotation=0,
         verticalalignment='center', color='black', fontsize=10)
## adding significance statements



ax2.scatter(alpha0['norm_pCA'], alpha0['predicted_pCA'], label="$\\beta=0$ (parent)", c="#00468BFF")
ax2.scatter(alpha25['norm_pCA'], alpha25['predicted_pCA'], label="$\\beta=2.5$ (parent)",c="#ED0000FF")
ax2.scatter(alpha5['norm_pCA'], alpha5['predicted_pCA'], label="$\\beta=5$ (parent)", c="#925E9FFF")
ax2.scatter(alpha_ood['norm_pCA'], alpha_ood['predicted_pCA'], label="$\\beta=2.5$ (top strain)",
           c="#ADB6B6FF")


ax2.plot([0, 1], [0, 1], transform=ax2.transAxes, linestyle="--", c="black")

ax2.set_xlabel("True normalized [p-CA]")
ax2.set_ylabel("Predicted [p-CA]")
# ax2.legend(loc="upper left", bbox_to_anchor=(1.0, 1))

plt.show()
# fig.savefig("figures/ValidationRound/exploration_exploitation_scenarios_plot.png", bbox_inches="tight")
# fig.savefig("figures/ValidationRound/exploration_exploitation_scenarios_plot.svg", bbox_inches="tight")
#




fig.savefig("figures/ValidationRound/exploration_exploitation_violins_plot.png", bbox_inches="tight")
fig.savefig("figures/ValidationRound/exploration_exploitation_violins_plot.svg", bbox_inches="tight")

#%%


_, _, r, p, _ = scipy.stats.linregress(list(alpha0['norm_pCA']), list(alpha0['predicted_pCA']))
print(r ** 2)
_, _, r, p, _ = scipy.stats.linregress(list(alpha25['norm_pCA']), list(alpha25['predicted_pCA']))
print(r ** 2)
_, _, r, p, _ = scipy.stats.linregress(list(alpha5['norm_pCA']), list(alpha5['predicted_pCA']))
print(r ** 2)
_, _, r, p, _ = scipy.stats.linregress(list(alpha_ood['norm_pCA']), list(alpha_ood['predicted_pCA']))
print(r ** 2)

spr_0 = spearmanr(list(alpha0['norm_pCA']), list(alpha0['predicted_pCA'])).statistic
spr_25 = spearmanr(list(alpha25['norm_pCA']), list(alpha25['predicted_pCA'])).statistic
spr_5 = spearmanr(list(alpha5['norm_pCA']), list(alpha5['predicted_pCA'])).statistic
spr_ood = spearmanr(list(alpha_ood['norm_pCA']), list(alpha_ood['predicted_pCA'])).statistic

print(spr_0, spr_25, spr_5, spr_ood)

all_scenarios = pd.concat([alpha5,alpha25, alpha0, alpha_ood])
_, _, r, p, _ = scipy.stats.linregress(list(all_scenarios['norm_pCA']), list(all_scenarios['predicted_pCA']))
print(r**2)
spearman_all = spearmanr(list(all_scenarios['norm_pCA']), list(all_scenarios['predicted_pCA'])).statistic
print(spearman_all)