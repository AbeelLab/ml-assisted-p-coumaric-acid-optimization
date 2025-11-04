"""1. Performs parsimonious flux balance analysis on yeast8.
2.   """

 #%%
# Create a model
from cobra import Model, Reaction, Metabolite
from cobra.io import (load_json_model, save_json_model,
                      load_matlab_model, save_matlab_model,
                      read_sbml_model, write_sbml_model)
import regex as re
import numpy as np
import matplotlib.pyplot as plt
import cobra
import pandas as pd
# from brendapyrser import BRENDA



def find_metabolite(model,pattern):
    """Find metabolite in model"""
    met=""
    for metabolite in model.metabolites:
        name=model.metabolites.get_by_id(str(metabolite)).name
        if re.search(pattern,name)!=None:
            met=metabolite
    return met


def find_reaction(model,pattern):
    """Find reaction in model"""
    react=[]

    for reaction in model.reactions:
        name=reaction.name
        if re.search(pattern,name)!=None:
            react.append(reaction)
    return react

#load the model from Lu2019
model=read_sbml_model("data/raw/Yeast-GEM/yeast-GEM.xml")


print(len(model.species))

#%%
# add PAL route to model.
## We know that the PAL route is preferred over the TAL route, so we should add this to the stoichiometry

# PHE--> CIN+NH4: PAL1 (phenyl ammonia-lyase)
# CIN-->  p-coumarate:  C4H CINinnamate 4-hydroxylase
# CPR1: does not seem to have any effect and is not required?
PAL_reaction=Reaction("PAL1")
PAL_reaction.name="Phenylalanine ammonia-lyase 1"
#Find metabolite phenylalanine

#PAL reaction
phe=find_metabolite(model,"phenylalanine \[cytoplasm\]")
cin = Metabolite(
    'cin_c',
    formula='C9H7O2',
    name='cinnamate [cytoplasm]',
    compartment='c')
ammon=find_metabolite(model,"ammonium \[cytoplasm\]")
PAL_reaction.add_metabolites({
    phe:-1,
    cin:1,
    ammon:1})
model.add_reactions([PAL_reaction])


# C4H reaction
C4H_reaction=Reaction("C4H")
C4H_reaction.name="cinnamate 4-hydroxylase"
NADP_c=find_metabolite(model,"NADP\(\+\) \[cyt")
NADPH_c=find_metabolite(model,"NADPH \[cytoplasm]")
pCA = Metabolite(
    'p-CA_c',
    formula='C9H8O3',
    name='p-Coumaric acid [cytoplasm]',
    compartment='c')
C4H_reaction.add_metabolites({
    cin:-1,
    NADPH_c:-1,
    pCA:1,
    NADP_c:1})
model.add_reactions([C4H_reaction])

#add sink reaction
model.add_boundary(model.metabolites.get_by_id("p-CA_c"),type="sink")

#%%

### We then perform parsimonious Flux Balance Analysis
# with as an objective function the flux through c4H
# We add a minimization of flux constraint (the parsimonious part)
# Furthermore, we add the maintenance reaction to the objective function

model.objective={model.reactions.get_by_id("r_2111"):1}
pfba_solution_biomass = cobra.flux_analysis.pfba(model)

model.objective={model.reactions.get_by_id("C4H"):1, model.reactions.get_by_id("r_4046"):1}


# remove reactions that have no flux in optimized strain given yeast-GEM
x=model.optimize()
flux_reactions=np.where(np.abs(x.fluxes)>0.001)[0]
model.summary()
pfba_solution_pca = cobra.flux_analysis.pfba(model)

#%%


biomass_to_pca_diff = pfba_solution_pca.fluxes / pfba_solution_biomass.fluxes
# Change infinities to large quantities
biomass_to_pca_diff[np.isinf(biomass_to_pca_diff)] = 100000
interesting_targets = np.where(biomass_to_pca_diff > 0)[0]

# lists to be included in the .csv file
gene_names = []
dPCA = []
dBiomass = []
dPCAdBiomass = []
real_gene_name = []
EC_number = []

for i in interesting_targets:
    gene_names.append(model.reactions[i].name)
    dPCAdBiomass.append(biomass_to_pca_diff[i])
    dPCA.append(pfba_solution_pca.fluxes.values[i])
    dBiomass.append(pfba_solution_biomass.fluxes[i])

    real_gene = [k.name for k in list(model.reactions[i].genes)]
    real_gene_name.append(real_gene)
    try:

        EC_number.append(np.array(model.reactions[i].annotation['ec-code']))
    except:
        EC_number.append(np.array([]))

id_target_dict = dict(zip(gene_names, dPCAdBiomass))
names = ["Enzyme", "dPCA/dBiomass", "dPCA", "dBiomass", "symbolic_name", "EC_number"]
features = dict(zip(names,[gene_names, dPCAdBiomass, dPCA, dBiomass, real_gene_name, EC_number]))


id_target_df = pd.DataFrame(features)
id_target_df = id_target_df.sort_values(by="dPCA/dBiomass", ascending=False)
id_target_df['symbolic_name'] = [",".join(map(str, l)) for l in id_target_df['symbolic_name']]

# manually set gene name of the two genes that were not in the model already. Could also do that in the model itself probably, but this is easier.
id_target_df.loc[id_target_df['Enzyme'] == "cinnamate 4-hydroxylase", 'EC_number'] = "1.14.14.91"
id_target_df.loc[id_target_df['Enzyme'] == "cinnamate 4-hydroxylase", 'symbolic_name'] = "C4H"

id_target_df.loc[id_target_df['Enzyme'] == "Phenylalanine ammonia-lyase 1", 'EC_number'] = "4.3.1.24"
id_target_df.loc[id_target_df['Enzyme'] == "Phenylalanine ammonia-lyase 1", 'symbolic_name'] = "PAL"

#%%
# Plot the result
fig,ax=plt.subplots(figsize=(20,8))
ax.bar(id_target_df['Enzyme'],id_target_df['dPCA/dBiomass'],color='#00468BFF')
ax.set_ylim(0,3)

# Set xticks and rotate them
ax.set_xticks(range(len(id_target_df['Enzyme'])))
ax.set_xticklabels(id_target_df['Enzyme'], rotation=90)
ax.axhline(1,c="black",linestyle="--")
# ax.set_yscale("log",base=10)
ax.set_ylabel("Ratio dpCA/dBiomass")
# plt.show()
fig.savefig("figures/ChoosingGeneTargets/FBA_dpCAdBiomass_ration.png",bbox_inches='tight')
fig.savefig("figures/ChoosingGeneTargets/FBA_dpCAdBiomass_ration.svg",bbox_inches='tight')


#%%

#Relevant literature report
#Not every gene has an enzyme commission number.
# These are all sink reactions (modelling) or transport reactions (e.g. L-isoleucine transport).
# While BRENDA has interesting properties, I found the literature attached to
# be a bit too chemical and out of context of Metabolic Engineering.
# I therefore use STRING-db on a set of proteins to identify literature that is enriched in these targets.
#
# We perform the following procedure
#2. Choose enzymes where dPCA/dbiomass>1
#3. Perform string search on E.C. numbers . While you could use an API, I have found their GUI to be very useful.
# If a target has multiple E.C. numbers, only take the first one.
#4. Download .tsv file for relevant literature


# count at how many occurences the gene is found in the string database

Features_for_String=id_target_df[id_target_df['dPCA/dBiomass']>1]


#7 genes do not have an E.C. number attached (25 genes left).
#Using the GUI https://string-db.org, we find E.C. for 17 genes. One thing stringDB cannot do is find connection between multiple organisms, so two genes (PAL and C4H)
# need to be seperately checked for literature. So in total, for our targets, we find 19/25 features using stringDB.

# Now we use the literature that was defined. We separate PAL-TAL from the other targets, because these are from another organism
lit_data_other=pd.read_csv("data/GeneTargets/enrichmentOther.tsv",header=0,sep="\t")
lit_data_PALTAL=pd.read_csv("data/GeneTargets/enrichmentPAL_C4H.tsv",header=0,sep="\t")


#Now, we extract the number of mentions of each enzyme and divide by the total number of literature mentions. For PAL-C4H this will be one due to
#the way enrichment analysis works, but in this case that won't be a problem.
# for i,j in lit_data_other.iterrows():
#     string += j["matching proteins in your network (labels)"]
# string
mentions=','.join(lit_data_other["matching proteins in your network (labels)"])
# print(mentions)
mentions=mentions.split(",")

values, counts = np.unique(mentions, return_counts=True)
counts=counts/np.shape(lit_data_other)[0]
literature_dict=dict(zip(values,counts))

mentions=','.join(lit_data_PALTAL["matching proteins in your network (labels)"])
# print(mentions)
mentions=mentions.split(",")
values, counts = np.unique(mentions, return_counts=True)
counts=counts/np.shape(lit_data_PALTAL)[0] #both 1
literature_dict['C4H']=1
literature_dict['PAL']=1

lit_dict_series=pd.Series(literature_dict)

Features_for_String=Features_for_String.join(lit_dict_series.rename('Rel_Lit_References'),on="symbolic_name") #this doesnt work for double names so add manually
# Features_for_String
Features_for_String["Rel_Lit_References"][Features_for_String['symbolic_name']=='ARO3,ARO4']["Rel_Lit_References"]=0.363636
Features_for_String["Rel_Lit_References"][Features_for_String['symbolic_name']=='PFK2,PFK1']["Rel_Lit_References"]=0.345455
Features_for_String["Rel_Lit_References"][Features_for_String['symbolic_name']=='ERR2,ENO1,ERR3,ENO2,ERR1']["Rel_Lit_References"]=0.390909
Features_for_String["Rel_Lit_References"][Features_for_String['symbolic_name']=='SOL3,SOL4']["Rel_Lit_References"]=0.027273 # a homolog
Features_for_String["Rel_Lit_References"][Features_for_String['symbolic_name']=='TKL2,TKL1']["Rel_Lit_References"]=0.163636
Features_for_String["Rel_Lit_References"][Features_for_String['symbolic_name']=='GPM1,YOR283W']=0.054545 # a homolog



dataFile = 'data/GeneTargets/BRENDA/brenda_2023_1.txt'
brenda = BRENDA(dataFile)
EC_numbers=Features_for_String['EC_number']

Inhibitor_dictionary={}
for i in EC_number:
    if isinstance(i,str):
        try:
            r=brenda.reactions.filter_by_organism(species="Saccharomyces cerevisiae").get_by_id(i)
            Inhibitors=len(r.KIvalues.filter_by_organism("Saccharomyces cerevisiae").keys())
            # print(Inhibitors)
            Inhibitor_dictionary[i]=Inhibitors
        except:
            continue
    elif isinstance(i,list):
        for j in i:
            try:
                r=brenda.reactions.filter_by_organism(species="Saccharomyces cerevisiae").get_by_id(j)

                Inhibitors=len(r.KIvalues.filter_by_organism("Saccharomyces cerevisiae").keys())

                Inhibitor_dictionary[i]=Inhibitors
            except:
                continue


Features_for_String['Inhibitors']=np.zeros(np.shape(Features_for_String)[0])

for i in Inhibitor_dictionary.keys():
    x=Features_for_String['EC_number'].str.contains(i)
    pos=np.where(x==True)
    Features_for_String['Inhibitors'].iloc[pos]=Inhibitor_dictionary[i]


Features_for_String.to_csv("results/GeneTargets/selected_gene_targets.csv")