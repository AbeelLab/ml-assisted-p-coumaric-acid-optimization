"""Generates gene count plot"""


#%% gene inclusion fraction for the ai4bio library.

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Example: Tick names mapping for genes in pathway description
tick_names = {
    "Atha_C4H.orf": "C4H",
    "Atha_CPR1.orf": "CPR1",
    "Atha_PAL1.orf": "PAL1",
    # "Ec_AROL.orf": "AROL",
    "Sc_ARO1.orf": "ARO1",
    "Sc_ARO2.orf": "ARO2",
    "Sc_ARO4_OFP.orf": "ARO4",
    # "Sc_ARO7.orf": "ARO7",
    "Sc_ARO8.orf": "ARO8",
    "Sc_ENO2.orf": "ENO2",
    "Sc_FBA1.orf": "FBA1",
    "Sc_GND1": "GND1",
    "Sc_GND2.orf": "GND2",
    "Sc_PHA2.orf": "PHA2",
    "Sc_RKI1.orf": "RKI1",
    "Sc_SOL3.orf": "SOL3",
    "Sc_SOL4.orf": "SOL4",
    "Sc_TKL1.orf": "TKL1",
    "Sc_YHM2.orf": "YHM2",
    "Sc_ZWF1.orf": "ZWF1"
}

# Load the data
filepath=("data/processed/"
          "IntegratedData_WURTUD/strain_numeric_matrix.csv")
strain_numeric_matrix = pd.read_csv(filepath,
                                    index_col=0)
print(np.sum(strain_numeric_matrix == 0, axis=1).sum())

# Filter for Coum_Library strains and SHK0066 control
ai4bio_round = strain_numeric_matrix[strain_numeric_matrix.index.str.contains("Coum_Library", regex=True)]
control = strain_numeric_matrix[strain_numeric_matrix.index.str.contains("SHK0066", regex=True)]

# Create a control dictionary for subtraction
control_dict = dict(control.T)['SHK0066']
ai4bio_round = ai4bio_round.apply(
    lambda col: col - control_dict[col.name] if col.name in control_dict else col, axis=0
)

# Count genes that are nonzero
ai4bio_gene_count = np.sum(ai4bio_round != 0)

# Filter xs and ys based on pathway genes (keys in tick_names)
xs = [gene for gene in ai4bio_gene_count.index if gene in tick_names]  # Only genes in tick_names
ys = ai4bio_gene_count[xs].values / np.shape(ai4bio_round)[0]  # Normalize counts

# Create the plot
fig, ax = plt.subplots(figsize=(4, 4))
ax.bar(xs, ys,color="#00468BFF")

# Update x-tick labels based on tick_names
mapped_labels = [tick_names[gene] for gene in xs]  # Map gene names to pathway descriptions
ax.set_xticks(range(len(xs)))  # Set tick positions based on the number of bars
ax.set_xticklabels(mapped_labels, rotation=90)  # Update tick labels and rotate them
ax.set_ylabel("Fraction of genes")
ax.set_ylim(0,1)
# Add title and show the plot
ax.set_title("Gene abundance in strains")
plt.tight_layout()  # Adjust layout to fit rotated labels
fig.savefig("figures/library_visualization/gene_count_TUD_library.png",bbox_inches="tight",dpi=1200)
fig.savefig("figures/library_visualization/gene_count_TUD_library.svg",bbox_inches="tight")

plt.show()

