###############################################
##Dmitry Sutormin, 2021##
##CRISPR spacers data analysis##

#1) Script parses MMseqs2 clustering results, identifies clusters of CRISPR-Cas spacers.
#2) Outputs spacer clusters and intersection of spacer sets derived from different samples.
###############################################


#######
#Packages to be imported.
#######

import os
import numpy as np
import Bio
from Bio import SeqIO
import re
import matplotlib.pyplot as plt
from matplotlib_venn import venn2, venn2_circles, venn3, venn3_circles
import csv
import pandas as pd
import seaborn as sns


#######
## Part 1.
## Cd-hit-based clustering.
## Identify and count spacers shared between different datasets, prepare heatmap and link data for Circos.
#######

#Cd-hit clustering results for WS datasets.
Cdhit_data_path="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\Spongy\CRISPR_spacers\\2021\Cctyper_Spacers\Spacers_clustering_cdhit_2018\\all_spacers_cdhit_id95_c95_2018.clstr"

#Output for WS datasets:
Cdhit_output_path="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\Spongy\CRISPR_spacers\\2021\Cctyper_Spacers\Spacers_clustering_cdhit_2018\\"

#Cd-hit clustering results: WS vs CRISPRCasDB.
Cdhit_WS_vs_DB_data_path="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\Spongy\CRISPR_spacers\\2021\CRISPRCasDB_2021\Cctyper\Cdhit_2018\\CRISPRCasDB_and_all_spacers_cdhit_id95_c95_2018.clstr"

#Output for WS vs CRISPRCasDB:
Cdhit_WS_vs_DB_output_path="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\Spongy\CRISPR_spacers\\2021\CRISPRCasDB_2021\Cdhit_2018\\"


#######
#Read spacers fasta.
#######

def read_spacers_fasta_dict(spacers_fasta_path):
    spacers_dict={}
    spacers=open(spacers_fasta_path, 'r')
    for record in SeqIO.parse(spacers, "fasta"):
        spacer_id=str(record.name) 
        spacer_seq=str(record.seq) 
        spacers_dict[spacer_id]=spacer_seq
    spacers.close()
    
    
    print(len(spacers_dict))
    
    return spacers_dict

#######
#Read clustering results prepared by cd-hit, prepare dictionary of clusters.
#######

def read_cdhit_clusters_data(datapath):
    filein=open(datapath, 'r')
    Clusters_dict={}
    initiation=1
    for line in filein:
        line=line.rstrip('\r\n')
        if line[0]==">":
            if initiation!=1:
                cluster_list_sorted=sorted(cluster_list, key=lambda x: x[1])
                spacer_names=[]
                for spacer_name in cluster_list_sorted:
                    spacer_names.append(spacer_name[0])
                    
                Clusters_dict[cluster_name]=spacer_names        
            cluster_name=line.lstrip('>')  
            cluster_list=[]
            initiation=0
        else:
            line=line.split(' ')
            line[1]=line[1].lstrip('>').rstrip('...')
            cluster_list.append(line[1:3])
    
    print('Number of clusters: ' + str(len(Clusters_dict)))   
    #print(Clusters_dict)
        
    filein.close()
    
    return Clusters_dict

#######
#Analyse dictionary of clusters: clusters size distribution.
#######

def analyse_clusters_dict(Clusters_dict, name, output_path):
    Clusters_len=[]
    for cluster_rep, cluster_dep_ar in Clusters_dict.items():
        Clusters_len.append(len(cluster_dep_ar))
        if len(cluster_dep_ar)>2:
            print(f'{cluster_rep} : {cluster_dep_ar}')
    
    
    #Plot distribution of clusters size.
    plt.hist(Clusters_len)
    plt.yscale('log')
    plt.xlabel('Size of clusters, 95% id')
    plt.ylabel('Number of clusters')
    plt.show()
    plt.savefig(f'{output_path}Cdhit_clusters_size_distribution_{name}.png')    

    return

#######
#Draw heatmap.
#######  

def draw_heatmap(Matrix_of_shared, keys_list, name, outpath, size):
    #Visualize data with heatmap.
    fig=plt.figure(figsize=(size,size), dpi=100)
    ax=fig.add_subplot(111)
    #Based on: https://matplotlib.org/gallery/images_contours_and_fields/image_annotated_heatmap.html#sphx-glr-gallery-images-contours-and-fields-image-annotated-heatmap-py
    im=ax.imshow(Matrix_of_shared, cmap='gist_yarg')
    ax.set_xticks(np.arange(len(keys_list)))
    ax.set_yticks(np.arange(len(keys_list)))
    ax.set_xticklabels(keys_list, size=10, rotation=90)
    ax.set_yticklabels(keys_list, size=10)  
    ax.set_ylim(sorted(ax.get_xlim(), reverse=True)) #Solves a bug in matplotlib 3.1.1 discussed here: https://stackoverflow.com/questions/56942670/matplotlib-seaborn-first-and-last-row-cut-in-half-of-heatmap-plot
    
    # Minor ticks taken from https://stackoverflow.com/questions/38973868/adjusting-gridlines-and-ticks-in-matplotlib-imshow
    ax.set_xticks(np.arange(-.5, len(keys_list), 1), minor=True)
    ax.set_yticks(np.arange(-.5, len(keys_list), 1), minor=True)
    
    ax.tick_params(axis='both', which='minor', length=0)
    
    # Gridlines based on minor ticks
    ax.grid(which='minor', color='grey', linestyle='-', linewidth=1)    
    
    for i in range(len(keys_list)):
        for j in range(len(keys_list)):
            if Matrix_of_shared.at[keys_list[i], keys_list[j]]>0:
                if isinstance(Matrix_of_shared.at[keys_list[i], keys_list[j]], int):
                    text = ax.text(j, i, Matrix_of_shared.at[keys_list[i], keys_list[j]], ha="center", va="center", color="k", size=7)
                else:
                    text = ax.text(j, i, int(Matrix_of_shared.at[keys_list[i], keys_list[j]]), ha="center", va="center", color="k", size=7)
    fig.tight_layout()
    plt.show()
    plt.savefig(f'{outpath}{name}_heatmap.png', dpi=400, figsize=(size, size))
    plt.savefig(f'{outpath}{name}_heatmap.svg', dpi=400, figsize=(size, size))
    return

#######
#Draw venn diagram.
#######  

def draw_venn(Matrix_of_shared, keys_list, name, outpath, size):
    #Visualize data with heatmap.
    fig=plt.figure(figsize=(size,size), dpi=100)
    ax=fig.add_subplot(111)
    
    venn2(subsets=(Matrix_of_shared.at[keys_list[0], keys_list[0]]-Matrix_of_shared.at[keys_list[1], keys_list[0]], Matrix_of_shared.at[keys_list[1], keys_list[1]]-Matrix_of_shared.at[keys_list[1], keys_list[0]], Matrix_of_shared.at[keys_list[1], keys_list[0]]), set_labels=set(keys_list))
    venn2_circles(subsets=(Matrix_of_shared.at[keys_list[0], keys_list[0]]-Matrix_of_shared.at[keys_list[1], keys_list[0]], Matrix_of_shared.at[keys_list[1], keys_list[1]]-Matrix_of_shared.at[keys_list[1], keys_list[0]], Matrix_of_shared.at[keys_list[1], keys_list[0]]))    
    
    fig.tight_layout()
    plt.show()
    plt.savefig(f'{outpath}{name}_venn.png', dpi=400, figsize=(size, size))
    plt.savefig(f'{outpath}{name}_venn.svg', dpi=400, figsize=(size, size))    
    return

#######
#For any number of datasets.
#Count overlap of datasets.
#######

def count_spacer_set_overlap(Clusters_dict, outpath):
    
    #Grouping by sample.
    #List_of_all_sets=['H_palmata_2016', 'H_palmata_2018', 'VRS_H_palmata_2018', 'H_panicea_2016', 'H_panicea_2018', 'VRS_H_panicea_2018', 'H_sitiens_2016', 'H_sitiens_2018', 'VRS_H_sitiens_2018', 'MW_2016', 'MW_2018', 'VRS_MW_2018']
    List_of_all_sets=[ 'H_panicea_2018', 'VRS_H_panicea_2018', 'H_sitiens_2018', 'VRS_H_sitiens_2018', 'H_palmata_2018', 'VRS_H_palmata_2018','MW_2018', 'VRS_MW_2018']
    
    
    spacer_all_sets_overlap_matrix_non_clustered=pd.DataFrame(0, index=List_of_all_sets, columns=List_of_all_sets)
    spacer_all_sets_overlap_matrix_clustered=pd.DataFrame(0,     index=List_of_all_sets, columns=List_of_all_sets)
    
    for cluster_id, spacer_set in Clusters_dict.items():
        sample_id_list=[]
        for spacer_id in spacer_set:
            sample_id=spacer_id.split('_NODE_')[0]
            sample_id_list.append(sample_id)
        sample_id_set=list(set(sample_id_list))  
        
        for i in range(len(sample_id_set)):
            #Different spacers counting (with clustering).
            spacer_all_sets_overlap_matrix_clustered.at[sample_id_set[i], sample_id_set[i]]+=1
            #All spacers counting (without clustering).
            spacer_all_sets_overlap_matrix_non_clustered.at[sample_id_set[i], sample_id_set[i]]+=sample_id_list.count(sample_id_set[i])
            
            for j in range(len(sample_id_set)):
                if j>i:
                    #Sort elements in the pair to match the order of elements in the matrix to prepare a triangle matrix.
                    sample_id_pair=[[sample_id_set[i], List_of_all_sets.index(sample_id_set[i])], [sample_id_set[j], List_of_all_sets.index(sample_id_set[j])]]
                    sample_id_pair_sorted=sorted(sample_id_pair, key=lambda x: x[1], reverse=True)
                    #Different spacers counting (with clustering).
                    spacer_all_sets_overlap_matrix_clustered.at[sample_id_pair_sorted[0][0], sample_id_pair_sorted[1][0]]+=1
                    #All spacers counting (without clustering).
                    spacer_all_sets_overlap_matrix_non_clustered.at[sample_id_pair_sorted[0][0], sample_id_pair_sorted[1][0]]+=(sample_id_list.count(sample_id_pair_sorted[0][0])+sample_id_list.count(sample_id_pair_sorted[1][0]))
                            
    draw_heatmap(spacer_all_sets_overlap_matrix_clustered,     List_of_all_sets, 'All_datasets_clustered',     outpath, 5)
    draw_heatmap(spacer_all_sets_overlap_matrix_non_clustered, List_of_all_sets, 'All_datasets_not_clustered', outpath, 5)
        
        
    #Grouping by location.
    List_of_loc_sets=['H_panicea', 'H_sitiens', 'H_palmata', 'MW']
    
    spacer_loc_sets_overlap_matrix_non_clustered=pd.DataFrame(0, index=List_of_loc_sets, columns=List_of_loc_sets)
    spacer_loc_sets_overlap_matrix_clustered=pd.DataFrame(0,     index=List_of_loc_sets, columns=List_of_loc_sets)
    
    for cluster_id, spacer_set in Clusters_dict.items():
        sample_id_list=[]
        for spacer_id in spacer_set:
            sample_id=spacer_id.split('_20')[0]
            sample_id=sample_id.lstrip('VRS_')
            sample_id_list.append(sample_id)
        sample_id_set=list(set(sample_id_list))  
        
        for i in range(len(sample_id_set)):
            #Different spacers counting (with clustering).
            spacer_loc_sets_overlap_matrix_clustered.at[sample_id_set[i], sample_id_set[i]]+=1
            #All spacers counting (without clustering).
            spacer_loc_sets_overlap_matrix_non_clustered.at[sample_id_set[i], sample_id_set[i]]+=sample_id_list.count(sample_id_set[i])
            
            for j in range(len(sample_id_set)):
                if j>i:
                    #Sort elements in the pair to match the order of elements in the matrix to prepare a triangle matrix.
                    sample_id_pair=[[sample_id_set[i], List_of_loc_sets.index(sample_id_set[i])], [sample_id_set[j], List_of_loc_sets.index(sample_id_set[j])]]
                    sample_id_pair_sorted=sorted(sample_id_pair, key=lambda x: x[1], reverse=True)
                    #Different spacers counting (with clustering).
                    spacer_loc_sets_overlap_matrix_clustered.at[sample_id_pair_sorted[0][0], sample_id_pair_sorted[1][0]]+=1
                    #All spacers counting (without clustering).
                    spacer_loc_sets_overlap_matrix_non_clustered.at[sample_id_pair_sorted[0][0], sample_id_pair_sorted[1][0]]+=(sample_id_list.count(sample_id_pair_sorted[0][0])+sample_id_list.count(sample_id_pair_sorted[1][0]))
                            
    draw_heatmap(spacer_loc_sets_overlap_matrix_clustered,     List_of_loc_sets, 'Locations_clustered',     outpath, 3)
    draw_heatmap(spacer_loc_sets_overlap_matrix_non_clustered, List_of_loc_sets, 'Locations_not_clustered', outpath, 3)
          
    
    #Grouping by year.
    #List_of_year_sets=['2016', '2018']
    #
    #spacer_year_sets_overlap_matrix_non_clustered=pd.DataFrame(0, index=List_of_year_sets, columns=List_of_year_sets)
    #spacer_year_sets_overlap_matrix_clustered=pd.DataFrame(0,     index=List_of_year_sets, columns=List_of_year_sets)
    #
    #for cluster_id, spacer_set in Clusters_dict.items():
    #    sample_id_list=[]
    #    for spacer_id in spacer_set:
    #        sample_id=spacer_id.split('_NODE_')[0].split('_')[-1]
    #        sample_id_list.append(sample_id)
    #    sample_id_set=list(set(sample_id_list))  
    #    
    #    for i in range(len(sample_id_set)):
    #        #Different spacers counting (with clustering).
    #        spacer_year_sets_overlap_matrix_clustered.at[sample_id_set[i], sample_id_set[i]]+=1
    #        #All spacers counting (without clustering).
    #        spacer_year_sets_overlap_matrix_non_clustered.at[sample_id_set[i], sample_id_set[i]]+=sample_id_list.count(sample_id_set[i])
    #        
    #        for j in range(len(sample_id_set)):
    #            if j>i:
    #                #Sort elements in the pair to match the order of elements in the matrix to prepare a triangle matrix.
    #                sample_id_pair=[[sample_id_set[i], List_of_year_sets.index(sample_id_set[i])], [sample_id_set[j], List_of_year_sets.index(sample_id_set[j])]]
    #                sample_id_pair_sorted=sorted(sample_id_pair, key=lambda x: x[1], reverse=True)
    #                #Different spacers counting (with clustering).
    #                spacer_year_sets_overlap_matrix_clustered.at[sample_id_pair_sorted[0][0], sample_id_pair_sorted[1][0]]+=1
    #                #All spacers counting (without clustering).
    #                spacer_year_sets_overlap_matrix_non_clustered.at[sample_id_pair_sorted[0][0], sample_id_pair_sorted[1][0]]+=(sample_id_list.count(sample_id_pair_sorted[0][0])+sample_id_list.count(sample_id_pair_sorted[1][0]))
    #                        
    #draw_heatmap(spacer_year_sets_overlap_matrix_clustered,     List_of_year_sets, 'Years_clustered',     outpath, 2)
    #draw_heatmap(spacer_year_sets_overlap_matrix_non_clustered, List_of_year_sets, 'Years_not_clustered', outpath, 2) 
    #draw_venn(spacer_year_sets_overlap_matrix_clustered,        List_of_year_sets, 'Years_clustered',     outpath, 3)
    #draw_venn(spacer_year_sets_overlap_matrix_non_clustered,    List_of_year_sets, 'Years_not_clustered', outpath, 3)
    
    #Grouping by bact/vir.
    List_of_BV_sets=['Bacterial', 'Viral']
    
    spacer_BV_sets_overlap_matrix_non_clustered=pd.DataFrame(0, index=List_of_BV_sets, columns=List_of_BV_sets)
    spacer_BV_sets_overlap_matrix_clustered=pd.DataFrame(0,     index=List_of_BV_sets, columns=List_of_BV_sets)
    
    for cluster_id, spacer_set in Clusters_dict.items():
        sample_id_list=[]
        for spacer_id in spacer_set:
            if 'VRS_' in spacer_id:
                sample_id_list.append('Viral')
            else:
                sample_id_list.append('Bacterial')
        sample_id_set=list(set(sample_id_list))  
        
        for i in range(len(sample_id_set)):
            #Different spacers counting (with clustering).
            spacer_BV_sets_overlap_matrix_clustered.at[sample_id_set[i], sample_id_set[i]]+=1
            #All spacers counting (without clustering).
            spacer_BV_sets_overlap_matrix_non_clustered.at[sample_id_set[i], sample_id_set[i]]+=sample_id_list.count(sample_id_set[i])
            
            for j in range(len(sample_id_set)):
                if j>i:
                    #Sort elements in the pair to match the order of elements in the matrix to prepare a triangle matrix.
                    sample_id_pair=[[sample_id_set[i], List_of_BV_sets.index(sample_id_set[i])], [sample_id_set[j], List_of_BV_sets.index(sample_id_set[j])]]
                    sample_id_pair_sorted=sorted(sample_id_pair, key=lambda x: x[1], reverse=True)
                    #Different spacers counting (with clustering).
                    spacer_BV_sets_overlap_matrix_clustered.at[sample_id_pair_sorted[0][0], sample_id_pair_sorted[1][0]]+=1
                    #All spacers counting (without clustering).
                    spacer_BV_sets_overlap_matrix_non_clustered.at[sample_id_pair_sorted[0][0], sample_id_pair_sorted[1][0]]+=(sample_id_list.count(sample_id_pair_sorted[0][0])+sample_id_list.count(sample_id_pair_sorted[1][0]))
                            
    draw_heatmap(spacer_BV_sets_overlap_matrix_clustered,     List_of_BV_sets, 'Bact_Vir_clustered',     outpath, 2)
    draw_heatmap(spacer_BV_sets_overlap_matrix_non_clustered, List_of_BV_sets, 'Bact_Vir_not_clustered', outpath, 2)     
    draw_venn(spacer_BV_sets_overlap_matrix_clustered,        List_of_BV_sets, 'Bact_Vir_clustered',     outpath, 3)
    draw_venn(spacer_BV_sets_overlap_matrix_non_clustered,    List_of_BV_sets, 'Bact_Vir_not_clustered', outpath, 3)    
                
    
    return 

#######
#Count overlap of sponge and CRISPRCasDB datasets.
#######

def count_spacer_set_overlap_sponge_CCDB(Clusters_dict, outpath):
        
    #Grouping by dataset type.
    List_of_DB_sets=['WS', 'CRISPRCasDB']
    
    spacer_DB_sets_overlap_matrix_non_clustered=pd.DataFrame(0, index=List_of_DB_sets, columns=List_of_DB_sets)
    spacer_DB_sets_overlap_matrix_clustered=pd.DataFrame(0,     index=List_of_DB_sets, columns=List_of_DB_sets)
    
    for cluster_id, spacer_set in Clusters_dict.items():
        sample_id_list=[]
        for spacer_id in spacer_set:
            if '_NODE_' in spacer_id:
                sample_id_list.append('WS')
            else:
                sample_id_list.append('CRISPRCasDB')
        sample_id_set=list(set(sample_id_list))  
        
        for i in range(len(sample_id_set)):
            #Different spacers counting (with clustering).
            spacer_DB_sets_overlap_matrix_clustered.at[sample_id_set[i], sample_id_set[i]]+=1
            #All spacers counting (without clustering).
            spacer_DB_sets_overlap_matrix_non_clustered.at[sample_id_set[i], sample_id_set[i]]+=sample_id_list.count(sample_id_set[i])
            
            for j in range(len(sample_id_set)):
                if j>i:
                    #Sort elements in the pair to match the order of elements in the matrix to prepare a triangle matrix.
                    sample_id_pair=[[sample_id_set[i], List_of_DB_sets.index(sample_id_set[i])], [sample_id_set[j], List_of_DB_sets.index(sample_id_set[j])]]
                    sample_id_pair_sorted=sorted(sample_id_pair, key=lambda x: x[1], reverse=True)
                    #Different spacers counting (with clustering).
                    spacer_DB_sets_overlap_matrix_clustered.at[sample_id_pair_sorted[0][0], sample_id_pair_sorted[1][0]]+=1
                    #All spacers counting (without clustering).
                    spacer_DB_sets_overlap_matrix_non_clustered.at[sample_id_pair_sorted[0][0], sample_id_pair_sorted[1][0]]+=(sample_id_list.count(sample_id_pair_sorted[0][0])+sample_id_list.count(sample_id_pair_sorted[1][0]))
                            
    draw_heatmap(spacer_DB_sets_overlap_matrix_clustered,     List_of_DB_sets, 'WS_vs_CRISPRCasDB_clustered',     outpath, 2)
    draw_heatmap(spacer_DB_sets_overlap_matrix_non_clustered, List_of_DB_sets, 'WS_vs_CRISPRCasDB_not_clustered', outpath, 2)     
    draw_venn(spacer_DB_sets_overlap_matrix_clustered,        List_of_DB_sets, 'WS_vs_CRISPRCasDB_clustered',     outpath, 3)
    draw_venn(spacer_DB_sets_overlap_matrix_non_clustered,    List_of_DB_sets, 'WS_vs_CRISPRCasDB_not_clustered', outpath, 3)    
                
    
    return 

#Cdhit_clusters=read_cdhit_clusters_data(Cdhit_data_path)
#analyse_clusters_dict(Cdhit_clusters, 'id95_c95_all_sponge_spacers_clust', Cdhit_output_path)
#count_spacer_set_overlap(Cdhit_clusters, Cdhit_output_path)

Cdhit_WS_vs_DB_clusters=read_cdhit_clusters_data(Cdhit_WS_vs_DB_data_path)
analyse_clusters_dict(Cdhit_WS_vs_DB_clusters, 'id95_c95_all_sponge_spacers_and_CRISPRCasDB_clust', Cdhit_WS_vs_DB_output_path)
count_spacer_set_overlap_sponge_CCDB(Cdhit_WS_vs_DB_clusters, Cdhit_WS_vs_DB_output_path)




#######
## Part 2.
## For mmseqs2-based clustering.
## Identify and count spacers shared between different datasets, prepare heatmap and link data for Circos.
#######

#MMseqs2 clustering results.
MMseq_data_path="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\Spongy\CRISPR_spacers\\2021\CRISPRCasDB_2021\Ccfinder\CrisprCasDB_and_sponge_all_spacers_clust_by_sample_res.tsv"

#Fasta files dict.
Spacers_fasta_dict={'I_palmata_2016' :      'C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\Spongy\CRISPR_spacers\\2021\CrisprCasFinder_Spacers\I_palmata_2016_CrisprResult_del1_less_than_72.fasta',
                    'I_palmata_2018' :      'C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\Spongy\CRISPR_spacers\\2021\CrisprCasFinder_Spacers\I_palmata_2018_CrisprResult_del1_less_than_72.fasta',
                    'VRS_I_palmata_2018' :  'C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\Spongy\CRISPR_spacers\\2021\CrisprCasFinder_Spacers\VRS_I_palmata_2018_CrisprResult_del1_less_than_72.fasta',
                    'H_panicea_2016' :      'C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\Spongy\CRISPR_spacers\\2021\CrisprCasFinder_Spacers\H_panicea_2016_CrisprResult_del1_less_than_72.fasta',
                    'H_panicea_2018' :      'C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\Spongy\CRISPR_spacers\\2021\CrisprCasFinder_Spacers\H_panicea_2018_CrisprResult_del1_less_than_72.fasta',
                    'VRS_H_panicea_2018' :  'C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\Spongy\CRISPR_spacers\\2021\CrisprCasFinder_Spacers\VRS_H_panicea_2018_CrisprResult_del1_less_than_72.fasta',
                    'H_sitiens_2016' :      'C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\Spongy\CRISPR_spacers\\2021\CrisprCasFinder_Spacers\H_sitiens_2016_CrisprResult_del1_less_than_72.fasta',
                    'H_sitiens_2018' :      'C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\Spongy\CRISPR_spacers\\2021\CrisprCasFinder_Spacers\H_sitiens_2018_CrisprResult_del1_less_than_72.fasta',
                    'VRS_H_sitiens_2018' :  'C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\Spongy\CRISPR_spacers\\2021\CrisprCasFinder_Spacers\VRS_H_sitiens_2018_CrisprResult_del1_less_than_72.fasta',
                    'MW_2016' :             'C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\Spongy\CRISPR_spacers\\2021\CrisprCasFinder_Spacers\MW_2016_CrisprResult_del1_less_than_72.fasta',
                    'MW_2018' :             'C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\Spongy\CRISPR_spacers\\2021\CrisprCasFinder_Spacers\MW_2018_CrisprResult_del1_less_than_72.fasta',
                    'VRS_MW_2018' :         'C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\Spongy\CRISPR_spacers\\2021\CrisprCasFinder_Spacers\VRS_MW_2018_CrisprResult_del1_less_than_72.fasta',}

#Fasta files dict.
Spacers_fasta_dict_2={'I_palmata' :  'C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\Spongy\CRISPR_spacers\\2021\CrisprCasFinder_Spacers\Spacers_clustering\I_palmata_all_spacers_clust.fasta',
                      'H_panicea' :  'C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\Spongy\CRISPR_spacers\\2021\CrisprCasFinder_Spacers\Spacers_clustering\H_panicea_all_spacers_clust.fasta',
                      'H_sitiens' :  'C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\Spongy\CRISPR_spacers\\2021\CrisprCasFinder_Spacers\Spacers_clustering\H_sitiens_all_spacers_clust.fasta',
                      'MW' :         'C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\Spongy\CRISPR_spacers\\2021\CrisprCasFinder_Spacers\Spacers_clustering\MW_all_spacers_clust.fasta',}

#Fasta files dict.
Spacers_fasta_dict_3={'CRISPRCasDB'  :  'C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\Spongy\CRISPR_spacers\\2021\CRISPRCasDB_2021\CRISPRCasDB_20210121_spacer_34.fasta',
                      'Sponge'       :  'C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\Spongy\CRISPR_spacers\\2021\CrisprCasFinder_Spacers\Spacers_clustering\All_spacers_clust_by_sample_type_clust.fasta',}

#Output:
Output_path="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\Spongy\CRISPR_spacers\\2021\CRISPRCasDB_2021\Ccfinder\\"



#######
#Read spacers fasta.
#######

def read_spacers_fasta_dict(spacers_fasta_dict):
    spacers_id_dict={}
    for set_name, file_path in spacers_fasta_dict.items():
        spacers_id_list=[]
        spacers=open(file_path, 'r')
        for record in SeqIO.parse(spacers, "fasta"):
            spacer_id=str(record.name) 
            spacers_id_list.append(spacer_id)
        spacers.close()
        
        spacers_id_dict[set_name]=spacers_id_list
    
    print(spacers_id_dict)
    
    return spacers_id_dict

#######
#Read clustering results, prepare dictionary of clusters.
#######

def read_clusters_data(datapath):
    filein=open(datapath, 'r')
    Clusters_dict={}
    for line in filein:
        line=line.rstrip('\r\n').split('\t')
        cluster_rep=line[0]
        cluster_dep=line[1]
        if cluster_rep in Clusters_dict:
            current_cluster_ar=Clusters_dict[cluster_rep]
            current_cluster_ar.append(cluster_dep)
            Clusters_dict[cluster_rep]=current_cluster_ar
        else:
            Clusters_dict[cluster_rep]=[cluster_dep]
    filein.close()
    
    print('Number of clusters: ' + str(len(Clusters_dict)))
    return Clusters_dict


#######
#Analyse dictionary of clusters: clusters size distribution.
#######

def analyse_clusters_dict(Clusters_dict, output_path):
    Clusters_len=[]
    for cluster_rep, cluster_dep_ar in Clusters_dict.items():
        Clusters_len.append(len(cluster_dep_ar))
        if len(cluster_dep_ar)>2:
            print(f'{cluster_rep} : {cluster_dep_ar}')
    
    
    #Plot distribution of clusters size.
    plt.hist(Clusters_len)
    plt.yscale('log')
    plt.xlabel('Size of clusters, 95% id')
    plt.ylabel('Number of clusters')
    plt.show()
    plt.savefig(output_path + 'MMseqs2_clusters_size_distribution_id95_c95_all_sponge_spacers_clust_and_CrisprCasDB.png')    

    return


#######
#For any number of datasets.
#Count overlap of datasets.
#######

def count_spacer_set_overlap(Clusters_dict):
    
    #List_of_sets=['VRS_MW_2018', 'VRS_H_sitiens_2018', 'VRS_H_panicea_2018', 'VRS_H_palmata_2018', 'MW_2018', 'MW_2016', 'H_sitiens_2018', 'H_sitiens_2016', 'H_panicea_2018', 'H_panicea_2016', 'H_palmata_2018', 'H_palmata_2016'] #For inital samples.
    #List_of_sets=['MW', 'H_sitiens', 'H_panicea', 'H_palmata']
    #List_of_sets=['VRS_MW_2018', 'VRS_H_sitiens_2018', 'VRS_H_panicea_2018', 'VRS_I_palmata_2018', 'MW_2018', 'MW_2016', 'H_sitiens_2018', 'H_sitiens_2016', 'H_panicea_2018', 'H_panicea_2016', 'I_palmata_2018', 'I_palmata_2016'] #For inital samples.    
    List_of_sets=['MW', 'H_sitiens', 'H_panicea', 'I_palmata']
    
    Dict_of_set_pairs={}
    for i in range(len(List_of_sets)):
        for j in range(len(List_of_sets)):
            if j>i:
                Pair=List_of_sets[i]+'@'+List_of_sets[j]
                Dict_of_set_pairs[Pair]=0
    
    print(Dict_of_set_pairs) 
    
    for rep_spacer, dep_spacers in Clusters_dict.items():
        for dep_spacer in dep_spacers:
            rep_spacer_set_id=rep_spacer.split('__')[0] #split by '_NODE' for initial samples. '__' for pre-clustered samples. 
            dep_spacer_set_id=dep_spacer.split('__')[0] #split by '_NODE' for initial samples. '__' for pre-clustered samples.
            if rep_spacer_set_id!=dep_spacer_set_id:
                for pair in Dict_of_set_pairs.keys():
                    pair_ar=pair.split('@')
                    if rep_spacer_set_id in pair_ar:
                        if dep_spacer_set_id in pair_ar:
                            Dict_of_set_pairs[pair]+=1
    
    print(Dict_of_set_pairs)            
    
    return Dict_of_set_pairs


#######
#For 2 datasets.
#Count overlap of datasets.
#######

def count_spacer_set_overlap_2(Clusters_dict):

    i=0
    for rep_spacer, dep_spacers in Clusters_dict.items():
        for dep_spacer in dep_spacers:
            rep_spacer_set_id=rep_spacer.split('_NODE')
            dep_spacer_set_id=dep_spacer.split('_NODE')
            if len(rep_spacer_set_id)!=len(dep_spacer_set_id):
                print(rep_spacer, dep_spacers)
                i+=1
                
    print(i)            
    
    return {'CRISPRCasDB@Sponge' : i}


#######
#Identify matching spacers belong to different spacer sets. Prepare links for circos.
#######

def get_matching_spacers_make_circos_links(Clusters_dict, spacers_id_dict, output_path):
    
    links_file=open(output_path+'Circos\Sponge_ccfinder_all_spacers_clust_by_sample_links_circos.txt', 'w')
    
    for rep_spacer, dep_spacers in Clusters_dict.items(): 
        for dep_spacer in dep_spacers:
            rep_spacer_set_id=rep_spacer.split('__')[0] #split by '_NODE' for initial samples. '__' for pre-clustered samples.
            dep_spacer_set_id=dep_spacer.split('__')[0] #split by '_NODE' for initial samples. '__' for pre-clustered samples.  
            if rep_spacer_set_id!=dep_spacer_set_id:
                spacers_set_rep=spacers_id_dict[rep_spacer_set_id]
                spacers_set_dep=spacers_id_dict[dep_spacer_set_id]
                
                spacers_set_rep_coord=spacers_set_rep.index(rep_spacer)
                spacers_set_dep_coord=spacers_set_dep.index(dep_spacer)
                
                links_file.write(f'{rep_spacer_set_id} {spacers_set_rep_coord} {spacers_set_rep_coord+1} {dep_spacer_set_id} {spacers_set_dep_coord} {spacers_set_dep_coord+1}\n')
    
    links_file.close()
    
    return


#######
#Draw heatmap.
#######  

def draw_heatmap(Matrix_of_shared, keys_list, outpath):
    #Visualize data with heatmap.
    fig=plt.figure(figsize=(2,2), dpi=100)
    ax=fig.add_subplot(111)
    #Based on: https://matplotlib.org/gallery/images_contours_and_fields/image_annotated_heatmap.html#sphx-glr-gallery-images-contours-and-fields-image-annotated-heatmap-py
    im=ax.imshow(Matrix_of_shared, cmap='gist_yarg')
    ax.set_xticks(np.arange(len(keys_list)))
    ax.set_yticks(np.arange(len(keys_list)))
    ax.set_xticklabels(keys_list, size=10, rotation=90)
    ax.set_yticklabels(keys_list, size=10)  
    ax.set_ylim(sorted(ax.get_xlim(), reverse=True)) #Solves a bug in matplotlib 3.1.1 discussed here: https://stackoverflow.com/questions/56942670/matplotlib-seaborn-first-and-last-row-cut-in-half-of-heatmap-plot
    
    # Minor ticks taken from https://stackoverflow.com/questions/38973868/adjusting-gridlines-and-ticks-in-matplotlib-imshow
    ax.set_xticks(np.arange(-.5, len(keys_list), 1), minor=True)
    ax.set_yticks(np.arange(-.5, len(keys_list), 1), minor=True)
    
    ax.tick_params(axis='both', which='minor', length=0)
    
    # Gridlines based on minor ticks
    ax.grid(which='minor', color='grey', linestyle='-', linewidth=1)    
    
    for i in range(len(keys_list)):
        for j in range(len(keys_list)):
            if Matrix_of_shared.at[keys_list[i], keys_list[j]]>0:
                if isinstance(Matrix_of_shared.at[keys_list[i], keys_list[j]], int):
                    text = ax.text(j, i, Matrix_of_shared.at[keys_list[i], keys_list[j]], ha="center", va="center", color="k", size=7)
                else:
                    text = ax.text(j, i, int(Matrix_of_shared.at[keys_list[i], keys_list[j]]), ha="center", va="center", color="k", size=7)
    fig.tight_layout()
    plt.show()
    plt.savefig(f'{outpath}Sponge_ccfinder_all_spacers_clust_and_CrisprCasFinder_heatmap.png', dpi=400, figsize=(2, 2))
    plt.savefig(f'{outpath}Sponge_ccfinder_all_spacers_clust_and_CrisprCasFinder_heatmap.svg', dpi=400, figsize=(2, 2))
    return


#######
#Identify matching spacers belong to different spacer sets. Prepare links for circos.
#######

def spacer_sets_overlap_heatmap(Dict_of_set_pairs, spacers_id_dict, outpath):
    
    spacer_sets_overlap_matrix=pd.DataFrame(0, index=spacers_id_dict.keys(), columns=spacers_id_dict.keys())
    
    for set_id, spacer_set in spacers_id_dict.items():
        spacer_sets_overlap_matrix.at[set_id, set_id]=len(spacer_set)
    
    for pair, num_spacers in Dict_of_set_pairs.items():
        pair_ar=pair.split('@')
        spacer_sets_overlap_matrix.at[pair_ar[0], pair_ar[1]]=num_spacers
    
    print(spacer_sets_overlap_matrix)  
    draw_heatmap(spacer_sets_overlap_matrix, list(spacers_id_dict.keys()), outpath)
    
    
    return

#Spacers_id_dict=read_spacers_fasta_dict(Spacers_fasta_dict_3)
#MMseq_Clusters_dict=read_clusters_data(MMseq_data_path)
#analyse_clusters_dict(MMseq_Clusters_dict, Output_path)
#Dict_of_set_pairs=count_spacer_set_overlap(MMseq_Clusters_dict)
#Dict_of_set_pairs=count_spacer_set_overlap_2(MMseq_Clusters_dict) #To find spacer shared between 2 datasets.
#get_matching_spacers_make_circos_links(MMseq_Clusters_dict, Spacers_id_dict, Output_path)
#spacer_sets_overlap_heatmap(Dict_of_set_pairs, Spacers_id_dict, Output_path)




#######
## Part 3.
## For mmseqs2-based clustering.
## Read output of mmseqs2, get representative sequences, return fasta for representative sequences.
#######

Sample_set_ID=""
#MMseqs2 clustering results.
MMseq_data_inpath=f'C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\Spongy\CRISPR_spacers\\2021\CrisprCasFinder_Spacers\Spacers_clustering\\{Sample_set_ID}All_spacers_clust_by_sample_res.tsv'


#######
#Read all input fasta files, prepare dict of fasta sequences.
#######

def read_all_fasta(spacers_fasta_dict):
    
    spacers_seq_sets_dict={}
    for set_name, file_path in spacers_fasta_dict.items():
        spacers_seq_dict={}
        spacers=open(file_path, 'r')
        for record in SeqIO.parse(spacers, "fasta"):
            spacer_id=str(record.name) 
            spacer_seq=str(record.seq) 
            spacers_seq_dict[spacer_id]=spacer_seq
        spacers.close()
        
        spacers_seq_sets_dict[set_name]=spacers_seq_dict
    
    return spacers_seq_sets_dict


#######
#Return sequence by spacer ID, write final fasta file of representative spacer sequences.
#######

def write_final_spacer_seqs(Clusters_dict, spacers_seq_sets_dict, sample_set_ID, output_path):
    
    fileout=open(f'{output_path}{sample_set_ID}_all_spacers_clust_by_sample_type_clust.fasta', 'w')

    for Rep_seq_id, Dep_seq_ar in Clusters_dict.items():
        rep_spacer_set_id=Rep_seq_id.split('__')[0] #Split by '_NODE' for initial files. '__' for pre-clustered files.
        rep_spacer_seq=spacers_seq_sets_dict[rep_spacer_set_id][Rep_seq_id]
        
        #fileout.write(f'>{sample_set_ID}__{Rep_seq_id}\n{rep_spacer_seq}\n') #For initial files.
        fileout.write(f'>{Rep_seq_id}\n{rep_spacer_seq}\n') #For pre-clustered files.
            
    fileout.close()
    
    return 

#Spacers_seq_sets_dict=read_all_fasta(Spacers_fasta_dict_2)
#MMseq_Clusters_dict=read_clusters_data(MMseq_data_inpath)
#write_final_spacer_seqs(MMseq_Clusters_dict, Spacers_seq_sets_dict, Sample_set_ID, Output_path)



#######
## Part 4. 
## For mmseqs2-based clustering.
## Required for proper clustering of sequences (mmseqs2 does not match the reverse complement sequences).
## Prepare files containing reverse-complement sequences (R).
## Merge F and R files.
#######

#Input folder to prepare files with F and R sequences.
Inoutpath="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\Spongy\CRISPR_spacers\\2021\CRISPRCasDB_2021\Cctyper\\"

def read_F_make_R_merge(pathin):
    for file in os.listdir(pathin):
        print(file)
        if '.fasta' in file:
            spacers_F=open(f'{pathin}\{file}', 'r')
            spacers_R=open(f'{pathin}\{file.strip(".fasta")}_R.fasta', 'w')
            spacers_FR=open(f'{pathin}\{file.strip(".fasta")}_FR.fasta', 'w')
            
            for record in SeqIO.parse(spacers_F, "fasta"):
                spacer_id_F=str(record.name)+'_F' 
                spacer_seq_F=str(record.seq)
                spacer_id_R=str(record.name)+'_R' 
                spacer_seq_R=str(record.seq.reverse_complement())
                
                spacers_R.write(f'>{spacer_id_R}\n{spacer_seq_R}\n')
                spacers_FR.write(f'>{spacer_id_F}\n{spacer_seq_F}\n>{spacer_id_R}\n{spacer_seq_R}\n')
            
            spacers_F.close()
            spacers_R.close()
            spacers_FR.close()
            
    return

#read_F_make_R_merge(Inoutpath)
