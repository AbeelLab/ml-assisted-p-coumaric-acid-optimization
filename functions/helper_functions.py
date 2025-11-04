import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os



def get_attribute_dataframe(gff_file):
    # Create a DataFrame from the list of attribute dictionaries
    attributes_list = []
    
    # Iterate over each contig's attributes in the GFF file
    for contig_attributes in gff_file['attributes']:
        # print(contig_attributes)

        attributes = contig_attributes.split(";")
        attribute_dict = {}
        
        # Extract key-value pairs and store them in the attribute dictionary
        for attribute in attributes:
            key, value = attribute.split("=")
            attribute_dict[key] = value
        attributes_list.append(attribute_dict)
    
    # Create a DataFrame from the list of attribute dictionaries
    attributes_file = pd.DataFrame(attributes_list)
    return attributes_file


def add_contig_depth_assembly_info(gff_file,contig_coverage_file):
    """adds contig lengths and mean depth to gff file """
    cov_file=pd.read_csv(contig_coverage_file,sep="\t")
    lengths_contig=cov_file['endpos'].to_list()
    mean_depth=cov_file['meandepth'].to_list()
    names=cov_file['#rname']
    contig_lengths=dict(zip(names,lengths_contig))
    mean_depth=dict(zip(names,mean_depth))
    temp=pd.DataFrame([contig_lengths,mean_depth],index=["contig_lengths","mean_depth"]).T
    gff_file=gff_file.join(temp,how="left",on="seqid")

    return gff_file

# def get_descriptions(cassettes):
#     # Return the DataFrame containing descriptions
#     # Define keys for the elements in the description
#     descriptions = []
    
#     # Iterate over each description in the 'description' column of the cassettes DataFrame
#     for description in cassettes['description']:
#         description = description.split("-")
        
#         description_dict = {}
        
#         # Enumerate over elements in the description
#         for k, element in enumerate(description):
#             # if element=='5':
#                 # description_dict['flank']=element
#             if "pro" in element:
#                 description_dict['promoter']=element
#             elif "orf" in element:
#                 description_dict['CDS']=element
#             elif "ter" in element:
#                 description_dict['terminator']=element                
#             # else:
#                 # description_dict['connector']=element
#             else:
#                 pass
#         descriptions.append(description_dict)
#     descriptions = pd.DataFrame(descriptions)
#     return descriptions

# def find_pots(gff_file):
#     """Find pots and add to gff dataframe as a labelling """
#     POTS = []

#     for i in range(len(gff_file['type']) - 2):
#         # Check for terminator, CDS, and promoter in consecutive entries
#         if (gff_file['type'][i] == "terminator" and
#                 gff_file['type'][i + 1] == "CDS" and
#                 gff_file['type'][i + 2] == "promoter"):
#             POTS.append([i+2, i + 1, i])
#     for i in range(len(gff_file['type']) - 2):
#         # Check for terminator, CDS, and promoter in consecutive entries but in reverse order
#         if (gff_file['type'][i + 2] == "terminator" and
#                 gff_file['type'][i + 1] == "CDS" and
#                 gff_file['type'][i] == "promoter"):

#             POTS.append([i, i + 1, i+2])
#     #sanity check, if length of unique elements is not equal to len of the number of pots *3, then there is one counted twice
#     # if len(POTS)*3!=len(np.unique(POTS)):
#     #     print(len(POTS))
#     #     print(len(np.unique(POTS)))
#     #     print("Warning, not all POTS have unique elements!")
#     # else:
#     POTS_labels=np.full(len(gff_file['type']),np.nan)
#     for i in range(len(POTS)):
#         POTS_labels[POTS[i]]=i
#     return POTS,POTS_labels

def find_pots_and_integration_site(gff_file):
    """Find pots and  integration sites and add to gff dataframe as a labelling """
    POTS = []

    # First the integration sites will be added as the first label (label=0)
    intsites=["Sc_INT_028v1_FLANK5.acc","Sc_INT69B_FLANK5.acc","Sc_INT72v1_FLANK5.acc"]
    # integration_site=np.where(gff_file['Name']=="Sc_INT_028v1_FLANK5.acc")[0]
    integration_sites=[]
    multiple_copies=[]
    for i in intsites:
        integration_site=np.where(gff_file['Name']==i)[0]
        if len(integration_site)>1:
            multiple_copies.append(integration_site)
        integration_sites.append(integration_site)
    integration_sites=[item for sublist in integration_sites for item in sublist]
    
    POTS.append(integration_sites)
    
    #TO DO: add additional check that pro, ter en cds are not too far apart (max 25 bp)
    for i in range(len(gff_file['type']) - 2):
        # Check for terminator, CDS, and promoter in consecutive entries
        if (gff_file['type'][i] == "terminator" and
                gff_file['type'][i + 1] == "CDS" and
                gff_file['type'][i + 2] == "promoter"):
            # extra check dat Pro en Orf en Orf en Ter dicht genoeg bij elkaar zijn
            pro_orf_dist=gff_file['start'][i+2]-gff_file['end'][i+1] 
            orf_ter_dist=gff_file['start'][i+1]-gff_file['end'][i]
            if pro_orf_dist<25 and orf_ter_dist<25:
                POTS.append([i+2, i + 1, i])

    for i in range(len(gff_file['type']) - 2):
        # Check for terminator, CDS, and promoter in consecutive entries but in reverse order
        if (gff_file['type'][i + 2] == "terminator" and
                gff_file['type'][i + 1] == "CDS" and
                gff_file['type'][i] == "promoter"):
            pro_orf_dist=gff_file['start'][i]-gff_file['end'][i+1] 
            orf_ter_dist=gff_file['start'][i+1]-gff_file['end'][i+2]
            if pro_orf_dist<25 and orf_ter_dist<25: #check whether distance between elements is smaller than 25
                POTS.append([i, i + 1, i+2])

    POTS_labels=np.full(len(gff_file['type']),np.nan)
    for i in range(len(POTS)):
        POTS_labels[POTS[i]]=i

    return POTS,POTS_labels, multiple_copies


def get_gap_length(gff_file):
    """For each contig, get the length between two elements"""
    grouped_gff=gff_file.groupby(by=['seqid'])
    gaps=[]
    for name, group in grouped_gff:
        gaps.append(0)
        for i in range(1,np.shape(group)[0]):
            gap=group.iloc[i,:]['start']-group.iloc[i-1,:]['end']
            gaps.append(gap)
    return gaps

# def get_int_site_labeling2(gff_file,max_gap_size):
#     int_sites=np.where(gff_file['POT_label']==0)[0]
#     in_int_site=np.empty(np.shape(gff_file)[0])
#     in_int_site[:]=np.nan

#     for i in int_sites:
#         # print(i,gff_file.iloc[i,:]['strand'])
#         if gff_file.iloc[i,:]['strand']=="+":
#             for j in range(i, np.shape(gff_file)[0]):
#                 if gff_file.iloc[j, :]['gaps'] < max_gap_size:
#                     in_int_site[j] = j
#                 else:
#                     break  # Exit the loop when the condition is not met
#         if gff_file.iloc[i,:]['strand']=="-":
#             for j in range(i,0,-1):
#                 # print(j)
#                 if gff_file.iloc[j,:]['gaps'] < max_gap_size:
#                     in_int_site[j] = j
#                 else:
#                     break
#     in_int_site_labeling=in_int_site
#     return in_int_site_labeling
        
# for i in range(len(int_site)): #loop over all integration sites that we find.
#     # gff_file.iloc[int_site[i]]
#     if group.iloc[int_site[i]]['strand']=="+": #checks if we should check forward or not   
#         for i in range(int_site[i], np.shape(group)[0]):
#             #If condition is not met on contig, break 
#             if group.iloc[i, :]['gaps'] < max_gap_size:
#                 in_int_site[i] = i
#             else:
#                 break  # Exit the loop when the condition is not met
#     elif group.iloc[int_site[i]]['strand']=="-":
#         #reverse the search for maximum gap length. We have to CHECK FOR THE ONE ON THE RIGHT, THEREFORE i+1
#         for i in range(int_site[i]-1,0,-1):

#             if group.iloc[i+1,:]['gaps'] < max_gap_size:
#                 in_int_site[i] = i

            
def get_int_site_labeling(POT_matrix,max_gap_size):
    """Synteny analysis of the POTs
    INput: gff file, and a threshold for what the maximum gap size can be"""
    grouped_gff=POT_matrix.groupby(['seqid','strand'])
    in_int_site_labeling=[]
    list_of_groups=[]
    hits=0
    for name, group in grouped_gff:
        in_int_site=np.empty(np.shape(group)[0])
        in_int_site[:]=np.nan
        
        int_site=np.where(group['POT_label']==0)[0]
        if len(int_site) != 0:
            hits+=1
            for i in int_site:
                if group.iloc[i]['strand']=="+": 
                    #If int_site is positive, skip to first POT
                    for j in range(i+1,np.shape(group)[0]):#!
                        if group.iloc[j, :]['gaps'] < max_gap_size:
                            in_int_site[j] = j
                        else:
                            break  # Exit the loop when the condition is not met
                elif group.iloc[i]['strand']=="-":
                    #reverse the search for maximum gap length. We have to CHECK FOR THE ONE ON THE RIGHT, THEREFORE i+1

                    for j in range(i-1,-1,-1):
                        # print(j)
                        if group.iloc[j+1,:]['gaps'] < max_gap_size:
                            in_int_site[j] = j
                        else:
                            break
        group['synteny_int_site']=in_int_site
        list_of_groups.append(group)
    merged_df=pd.concat(list_of_groups,ignore_index=True)
    return merged_df
    #     in_int_site_labeling.append(in_int_site)
    # print(in_int_site_labeling) 
    # in_int_site_labeling=[item for sublist in in_int_site_labeling for item in sublist]
    # # print(len(in_int_site_labeling))
    

def get_POT_count_matrix(gff_file):
    """Generates POT count matrix, along with POT length"""
    # Define the sorting order for sorting by type
    sorting_order = {'promoter': 0, 'CDS': 1, 'terminator': 2}
    
    # Sort each group into Promoter-CDS-Terminator
    gff_file = gff_file.groupby('POT_label').apply(
        lambda x: x.sort_values(by='type', key=lambda y: y.map(sorting_order))
    ).reset_index(drop=True)
    
    # Group by POT_label and concatenate the 'Name' column values
    POTS_df = gff_file.groupby("POT_label").agg({'Name': lambda x: '-'.join(x)}).reset_index(drop=True)
    
    # Create a counts column to count the number of unique elements
    POT_counts = np.ones(np.shape(POTS_df)[0])
    POTS_df['counts'] = POT_counts

    # Calculate the start and end for each POT_label and add to df
    start = gff_file.groupby('POT_label')[['start']].min().reset_index(drop=True)
    end = gff_file.groupby('POT_label')[['end']].max().reset_index(drop=True)
    total_POT_length = end['end'] - start['start']
    POTS_df['POT_length'] = total_POT_length
    
    # Group by 'Name' and aggregate counts by sum and POT_length by mean
    POT_counts_df = POTS_df.groupby('Name').agg({"counts": 'sum', 'POT_length': 'mean'}).reset_index(level="Name")
    
    return POT_counts_df


def merge_pots(gff_file):
    # Sort each group into Promoter-CDS-Terminator
    sorting_order = {'promoter': 0, 'CDS': 1, 'terminator': 2}
    gff_file = gff_file.groupby('POT_label').apply(
        lambda x: x.sort_values(by='type', key=lambda y: y.map(sorting_order))).reset_index(drop=True)
    grouped_gff=gff_file.groupby(by=['seqid','POT_label'])

    POT_matrix=[]
    for name, group in grouped_gff:
        if 0 not in group['POT_label'].values: #make sure it is not integration site
            min_strand_pos=np.min(group['start'])
            max_strand_pos=np.max(group['end'])
            POT_length=max_strand_pos-min_strand_pos
            name_of_POT="-".join(group['Name'].values)
            temp_dict={'seqid':group['seqid'].values[0],'strand':group['strand'].values[0],'Name':name_of_POT,'start':min_strand_pos,
                    'end':max_strand_pos,"POT_length":POT_length,'contig_lengths':group['contig_lengths'].values[0],
                    "mean_depth":group['mean_depth'].values[0],'POT_label':group['POT_label'].values[0]}
            POT_matrix.append(temp_dict)
    POT_matrix=pd.DataFrame(POT_matrix)
    return POT_matrix



def sum_similar_POTs(POT_count_df):
    """A postprocessing step. We do not expect much effect from using a certain terminator.
      To decrease dimensionality, we can therefore sum elements that are already described by a unique promoter orf"""
    newnames=[]
    for POT_name in POT_count_df['Name']:
        temp=POT_name.split("-")
        new_name=[]
        for element in temp:
            if ".pro" in element:
                new_name.append(element)
            elif ".orf" in element: 
                new_name.append(element)
        new_name="-".join(new_name)
        newnames.append(new_name)
    POT_count_df['Newname']=newnames
    POT_counts_df = POT_count_df.groupby('Newname').agg({"counts": 'sum', 'POT_length': 'mean'}).reset_index(level="Newname")
    return POT_counts_df