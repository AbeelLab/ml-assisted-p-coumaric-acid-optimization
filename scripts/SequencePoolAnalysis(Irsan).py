"""Sequence pool analysis (Irsan) reproduction"""
import pandas as pd
from anndata import AnnData
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.colors import ListedColormap
from scipy.stats import ranksums


file_path = "data/raw/SequencePool(Irsan)"
cassette_counts =  pd.read_csv(f'{file_path}/AI4BIO_Pools_and_Host_SHK0066_NoSecMappings_'
                           'NotAllData_genome_copy_cohort.txt', sep='\t',
                           index_col=0,
                           na_values='NA')

cassette_counts = cassette_counts[cassette_counts["feature"] == "cassette"]
cassette_counts["sample"] = cassette_counts["file_name"].str.replace("_merged_copynumber.txt", "")
sample_names = cassette_counts["sample"].unique().tolist()



# Pivot counts to wide format
cassette_coverage_wide = (
    cassette_counts[["name", "description", "mean_cov", "sample"]]
    .groupby(["name", "description", "sample"])["mean_cov"].mean()
    .reset_index()
    .pivot(index=["name", "description"], columns="sample", values="mean_cov")
    .reset_index()
)
cassette_coverage_matrix = cassette_coverage_wide[sample_names].to_numpy()
cassette_coverage_matrix_names = cassette_coverage_wide["name"]
cassette_coverage_matrix_df = pd.DataFrame(
    cassette_coverage_matrix,
    index=cassette_coverage_matrix_names,
    columns=sample_names
)
# Genome coverage
read_stats = pd.read_csv(f"{file_path}/AI4BIO_Pools_and_Host_SHK0066_NoSecMappings_NotAllData_read_stats.txt", sep="\t")
genome_size = 12e6
read_stats["sample"] = read_stats["file"].str.replace("_chopped.fastq.gz", "")
genome_coverage = (read_stats.set_index("sample")["sum_len"] / genome_size).round()

# Normalize counts
cassette_fractions = cassette_coverage_matrix_df.divide(genome_coverage, axis=1)
cassette_percentages = (cassette_fractions * 100).round(1)


# Create pseudo ExpressionSet
genes = cassette_coverage_wide[["name", "description"]].copy()
brick_matrix = genes["description"].str.split("-", expand=True)
brick_matrix.columns = ["lead_con", "prom", "orf", "term", "tail_con"]
genes = pd.concat([genes, brick_matrix], axis=1)
genes["sum_counts"] = cassette_percentages.sum(axis=1).values
genes["name"] = pd.Categorical(genes["name"], categories=genes.sort_values("sum_counts", ascending=False)["name"])

samples = pd.DataFrame(index=sample_names)
samples["group"] = ["before" if name == "SHK0066" else "after" for name in sample_names]

# Create AnnData object
adata = AnnData(X=cassette_percentages.T.values, obs=samples, var=genes)
adata.obs_names = cassette_percentages.columns
adata.var_names = cassette_percentages.index

# Annotate sample label and con_pair
adata.obs["sample_label"] = pd.Categorical(
    adata.obs_names, categories=["Lib1", "Lib2", "Lib3", "Lib4", "SHK0066"]
)

adata.var["con_pair"] = adata.var["lead_con"] + "-" + adata.var["tail_con"]

promoter_info = pd.read_csv(f"{file_path}/promoters_seqbase.txt", sep="\t")[["name", "strength"]]

# Cassette design info
clm_design = pd.read_csv(f"{file_path}/Library_compositions_AI4Bio_Coumaric_acid_Nov2023.csv")
cassette_name_matrix = clm_design["Composition Names"].str.split(" - ", expand=True).iloc[:, 1:8]
cassette_design_long = (
    pd.concat([clm_design[["library"]], cassette_name_matrix], axis=1)
    .melt(id_vars=["library"], var_name="slot", value_name="description")
    .replace({
        "description": {
            "Spar_TEF1.Pro": "Spar_TEF1.pro",

            "Xx_KanMX.orf": "Shyg_KanMX.orf"
        }
    })
    .merge(adata.var, on="description", how="left")
    .merge(promoter_info, left_on="prom", right_on="name", how="left"))

cassette_design_long["strength"] = pd.Categorical(
    cassette_design_long["strength"], categories=["WEAK", "MEDIUM", "STRONG"]
)
cassette_design_long["pot"] = (
    cassette_design_long["prom"] + " -> " + cassette_design_long["orf"] + " -> " + cassette_design_long["term"]
)
cassette_design_long["con_pair"] = cassette_design_long["lead_con"] + "-" + cassette_design_long["tail_con"]

# Filter low abundance cassettes
filtered_mask = cassette_percentages["SHK0066"] < 1
cassette_ids_filtered = cassette_percentages.index[filtered_mask]

#expression data
expr_df = adata.to_df().T
expr_long = expr_df.T.reset_index().melt(id_vars=["index"], var_name="sample_label", value_name="expression")
expr_long = expr_long.set_index('sample_label')
# expr_df.reset_index()

print(expr_long)
expr_long = pd.merge(expr_long,adata.var,left_index=True,right_index=True)
# expr_long= expr_long.reset_index("index", inplace=True)

expr_long = expr_long[
    expr_long["name"].isin(cassette_ids_filtered)
].groupby(["index", "con_pair", "orf"])["expression"].sum().reset_index()

expr_long= expr_long.rename(columns={"index":"library"})
print(expr_long)

modulated_genes = ['Atha_C4H.orf', 'Atha_CPR1.orf', 'Atha_PAL1.orf', 'Sc_ARO4_OFP.orf_0001']
nonmodulated = np.setdiff1d(expr_long['orf'].unique(),modulated_genes)
nonmodulated =list(nonmodulated)
print(nonmodulated)


expr_long= expr_long[expr_long['con_pair']!="F-G"]
#%%




libraries=['Lib1','Lib2','Lib3','Lib4']

values ={}
for library in libraries:
    temp = expr_long[expr_long["library"] == library]


    modulated_mean = temp[temp['orf'].isin(modulated_genes)].mean(numeric_only=True)
    modulated_std = temp[temp['orf'].isin(modulated_genes)].std(numeric_only=True)

    nonmodulated_mean = temp[temp['orf'].isin(nonmodulated)].mean(numeric_only=True)
    nonmodulated_std = temp[temp['orf'].isin(nonmodulated)].std(numeric_only=True)

    values[library] = {'modulated_mean':float(modulated_mean.values),
                       'nonmodulated_mean':float(nonmodulated_mean.values),
                       'modulated_std':float(modulated_std.values),
                       'nonmodulated_std':float(nonmodulated_std.values)}

    # print('a', modulated)
    # print('b', nonmodulated)
    # values[library] = [modulated, nonmodulated]

barplot_values =pd.DataFrame(values).T
barplot_values['design'] = ['300 ng/kb','200 ng/kb', '100 ng/kb','50 ng/kb']

print(barplot_values)
# X locations for groups
x = np.arange(len(barplot_values.design))  # label locations
width = 0.35  # width of the bars

# Create the plot
fig, ax = plt.subplots(figsize=(8, 5))
ax.bar(x - width/2, barplot_values.modulated_mean/4, width,
       yerr=barplot_values.modulated_std/4,
       label='[ARO4, C4H, PAl1, CPR1] dosage change)')
ax.bar(x + width/2, barplot_values.nonmodulated_mean/14,width,
       yerr=barplot_values.modulated_std/14,
       label='constant dosage')

# Labels and formatting
ax.set_xlabel('Library')
ax.set_ylabel('Expression values')
ax.set_title('Gene dosage effect on sequence pool')
ax.set_xticks(x)
ax.set_xticklabels(barplot_values.design)
ax.legend()

plt.tight_layout()
plt.show()
#%%

plot_data = []

for library in libraries:
    temp = expr_long[expr_long["library"] == library]
    for gene in modulated_genes:
        val = temp[temp['orf'] == gene]['expression'].values
        val = val / len(val)
        for v in val:
            plot_data.append({'library': library,
                              'gene_type': 'Varied dosage',
                              'expression': v,})

    for gene in nonmodulated:
        val = temp[temp['orf'] == gene]['expression'].values
        val = val/len(val)
        for v in val:
            plot_data.append({'library': library,
                              'gene_type': 'Fixed dosage',
                              'expression': v,})

# Create DataFrame for plotting
violin_df = pd.DataFrame(plot_data)

for library in libraries:
    temp = violin_df[violin_df['library'] == library]
    a = temp[temp['gene_type'] == 'Modulated']['expression']
    b = temp[temp['gene_type'] == 'Nonmodulated']['expression']

    test = ranksums(a,b)
    print(library, test) #all are significant so I will not plot this






# Replace library names with more descriptive labels if needed
design_labels = ['300 ng/kb', '200 ng/kb', '100 ng/kb', '50 ng/kb']
library_mapping = dict(zip(libraries, design_labels))
violin_df['design'] = violin_df['library'].map(library_mapping)

# Plot
# plt.figure(figsize=(10, 6))
fig,ax = plt.subplots(figsize=(10, 5))

lancet_colors =['#00468BFF','#ED0000FF']
custom_cmap = ListedColormap(lancet_colors)

ax = sns.violinplot(data=violin_df,
               x='design',
               y='expression',
               hue='gene_type',
               split=True,
               inner='quart',  # adds a boxplot inside the violins
               palette=lancet_colors,
               cut=0,
               alpha=1)

# Get category order
# Formatting
ax.set_xlabel('Library')
ax.set_ylabel('Occurrence values')
ax.set_title('Dosage effect of cassettes in sequenced pool')
ax.legend(title='Gene set', loc="upper left")
plt.tight_layout()
plt.show()
fig.savefig('figures/SequencePoolAnalysis/SequencedPools.svg', bbox_inches='tight')
fig.savefig('figures/SequencePoolAnalysis/SequencedPools.png',bbox_inches="tight")