import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
import sys
from functions.helper_functions import *

print(os.getcwd())

work_dir="../../AI4Bio_392clones_assembly_and_annotations/"
cov_dir=work_dir+"Coverage_per_contig_zip/" #change this back at some point 

#A yeast assembly should be around 12mb, with 20% difference at most
fasta_size_bounds=[12000000*0.8,12000000*1.2]

files=os.listdir(work_dir)

sequence_length=[]
counts_dict_list=[]
colony_names=[]

#construct the proper gff file
for file in files:
    if file.endswith(".gff"):
        colony_name=file.replace("_SeqBase_bricks_lifted.gff","")
        colony_name=colony_name.replace("_1_","")
        colony_name=colony_name.replace("_2_","")
        colony_name=colony_name.replace("_3_","")
        contig_coverage_file=cov_dir+colony_name+"_coverage_per_contig.txt"

        #2 find assembly size and make sure it is between bounds
        if os.path.isfile(work_dir+colony_name+"_assembly.fasta"):
            assembly_size=os.path.getsize(work_dir+colony_name+"_assembly.fasta")
        elif os.path.isfile(work_dir+"_1_"+colony_name+"_assembly.fasta"):
            assembly_size=os.path.getsize(work_dir+"_1_"+colony_name+"_assembly.fasta")
        elif os.path.isfile(work_dir+"_2_"+colony_name+"_assembly.fasta"):
            assembly_size=os.path.getsize(work_dir+"_2_"+colony_name+"_assembly.fasta")
        elif os.path.isfile(work_dir+"_3_"+colony_name+"_assembly.fasta"):
            assembly_size=os.path.getsize(work_dir+"_3_"+colony_name+"_assembly.fasta")

        if fasta_size_bounds[0]<assembly_size<fasta_size_bounds[1]:
            colony_names.append(colony_name)
            print(colony_name)

            gff_file=pd.read_csv(work_dir+file,sep="\t",header=None,comment="#")
            # print(gff_file.columns)
            # attributes are added to the dataframe
            gff_file.columns=['seqid','source','type','start','end','score','strand','phase','attributes']

            attributes_gff_file=get_attribute_dataframe(gff_file)  #get all attributes
            gff_file=pd.concat([gff_file,attributes_gff_file],axis=1)
            gff_file=gff_file.drop("attributes",axis=1)
            gff_file=gff_file.groupby(by=["seqid","strand"]).apply(lambda x: x.sort_values(by='start')).reset_index(drop=True)  #sort by contig and strand info

            # information for the coverage per contig such as mean depth and contig lengths
            # are added for later filtering.
            gff_file=add_contig_depth_assembly_info(gff_file,contig_coverage_file)

            #Finds pots and the integration site
            POTS,POTS_labels,multiple_copies=find_pots_and_integration_site(gff_file)


            gff_file['POT_label']=POTS_labels
            gff_file=gff_file[gff_file['POT_label'].notna()]


            # If we have properly called POTS we can group by POT label
            # Then we merge these and get the minimum and maximum value of the strand
            # With the exception of the integration site
            POT_matrix=merge_pots(gff_file)


            #add the integration site again
            int_site=gff_file[gff_file['POT_label']==0]

            int_site=int_site.filter(items=["seqid","strand","Name","start","end","contig_lengths","mean_depth","POT_label"],axis=1)
            POT_matrix=pd.concat([POT_matrix,int_site])
            POT_matrix=POT_matrix.groupby(by=["seqid"]).apply(lambda x: x.sort_values(by='start')).reset_index(drop=True)

            # add the gaplength
            POT_matrix['gaps']=get_gap_length(POT_matrix)

            # POT_matrix=get_int_site_labeling(POT_matrix,150)
            # print("beforena",np.shap
            # e(POT_matrix))

            # POT_matrix=POT_matrix[~POT_matrix['synteny_int_site'].isna()]
            # print(np.shape(POT_matrix))
            #move filtering completely to the  end together with summing similar terminator pairs
            POT_matrix=POT_matrix[POT_matrix['contig_lengths']>150000]
            POT_matrix=POT_matrix[POT_matrix['mean_depth']>20]

            #add a count column and sum
            counts=np.ones(np.shape(POT_matrix)[0])
            POT_matrix['counts']=counts
            POT_matrix = POT_matrix.groupby('Name').agg({"counts": 'sum', 'POT_length': 'mean'}).reset_index(level="Name")

            #sum similar pots (often different terminators)
            POT_matrix=sum_similar_POTs(POT_count_df=POT_matrix)
            POT_matrix=POT_matrix.rename({"Newname":"Name"},axis=1)


            counts_dict=dict(zip(POT_matrix['Name'].values,POT_matrix['counts']))

            counts_dict_list.append(counts_dict)



strain_count_matrix=pd.DataFrame(counts_dict_list,index=colony_names)
strain_count_matrix=strain_count_matrix.fillna(value=0)

#rename this KanMX
strain_count_matrix=strain_count_matrix.rename({'Agos_lox_TEF1.pro-Shyg_KanMX.orf-Agos_TEF1.ter':"Agos_lox_TEF1.pro-Shyg_KanMX.orf-Agos_TEF1.ter_0002"},axis=1)

#drop duplicates
strain_count_matrix = strain_count_matrix.loc[:,~strain_count_matrix.columns.duplicated()].copy()



#remove strains that have no integrations after QC
#rename this KanMX
strain_count_matrix=strain_count_matrix.rename({'Agos_lox_TEF1.pro-Shyg_KanMX.orf-Agos_TEF1.ter':"Agos_lox_TEF1.pro-Shyg_KanMX.orf-Agos_TEF1.ter_0002"},axis=1)
# strain_count_matrix=strain_count_matrix.rename({'Agos_lox_TEF1.pro-Shyg_KanMX.orf-Agos_TEF1.ter':"Agos_lox_TEF1.pro-Shyg_KanMX.orf-Agos_TEF1.ter_0002"},axis=1)

print(np.shape(strain_count_matrix))
#drop duplicates
strain_count_matrix = strain_count_matrix.loc[:,~strain_count_matrix.columns.duplicated()].copy()
print(np.shape(strain_count_matrix))
#remove strains that have no integrations after QC
empty_strains=strain_count_matrix.T.columns[np.where(strain_count_matrix.sum(axis=1)==0)[0]]
print(empty_strains)
strain_count_matrix=strain_count_matrix.drop(labels=empty_strains,axis=0)

# filter out some pots that need to be filtered out manually because they were called by biobricks,
# but are not actual pro-orf combinations that wer designed
int_site_filtering = [
    "Sc_ACT1.pro","Sc_FPS1.pro-Sc_FPS1.orf","Sc_GCN4.pro-Sc_GAL83.orf",
    "Sc_GCY1.pro-Sc_RKI1.orf","Sc_GSY1.pro", "Sc_HHF2.pro-Sc_NCE103.orf",
    "Sc_HXT7.pro-Sc_HXT7.orf_0001","Sc_LEU2.orf","Sc_MAD1.pro-Sc_MIG1.orf",
    "Sc_MAL.pro-Sc_MAL31.orf","Sc_TDH1.pro-Sc_RPE1.orf", "Sc_TIR2.pro-Sc_AUS1.orf",
    "Sc_TPO2.pro-Sc_LEU2.orf","Sc_TPS1.pro-Sc_TPS1.orf","Sc_TRP1.pro-Sc_TRP1.orf",
    "Sc_URA3.pro-Sc_URA3.orf","Sc_VPS68.pro","Sc_VPS68.pro-Sc_ADH1.orf",
    "Sc_YKT6.pro-Sc_RCN1.orf","Sc_ZPR1.pro-Sc_XKS1.orf","Sc_ACT1.pro_0001-Sc_ARO7.orf"
]
strain_count_matrix = strain_count_matrix.drop(int_site_filtering,axis=1)

strain_count_matrix.to_csv("data/processed/IntegratedData_WURTUD/WUR_TUD_count_matrix.csv")