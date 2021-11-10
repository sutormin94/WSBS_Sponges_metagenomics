###############################################
##Dmitry Sutormin, 2021##
##Prediction of viral contigs data##

#1) Script compares the numbers of filtered predicted viral contigs across datasets.
#2) Reads results of contigs clustering and draws heatmap representing the datasets overlappings.
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
#Data to be used.
#######

#Folder with filtered datasets.
Datasets_path="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\Spongy\Viral_contigs\Viral_contigs_prediction_and_filtering\\2018_only\\"

#Output folder.
Output_path="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\Spongy\Viral_contigs\Viral_contigs_analysis\\"


def compare_number_of_viral_contigs(inpath, outpath):
    
    Set_name_ar=[]
    Cont_num_ar=[]
    for set_folder in os.listdir(inpath):
        print(set_folder)
        for file in os.listdir(f'{inpath}\{set_folder}'):
            if '.fasta' in file:
                Set_name_ar.append(set_folder)
                viral_contigs=open(f'{inpath}\{set_folder}\{file}', 'r')   
                i=0
                for record in SeqIO.parse(viral_contigs, "fasta"):
                    i+=1
                Cont_num_ar.append(i)
    
    print(Set_name_ar)
    print(Cont_num_ar)
    
    #Visualize data with barplot.
    size1=5
    size2=3
    X_coords=range(len(Set_name_ar))
    
    fig=plt.figure(figsize=(size1,size2), dpi=100)
    ax=fig.add_subplot(111)
    
    ax.bar(X_coords, Cont_num_ar, width=0.6, edgecolor='k', linewidth=0.6)
    ax.set_ylabel('Number of\nviral contigs', size=15)
    ax.set_xticks(X_coords)
    ax.set_xticklabels(Set_name_ar, rotation=90, size=8)   
    ax.tick_params(axis='x', which='major', pad=0.5)    
    
    fig.tight_layout()
    plt.show()
    plt.savefig(f'{outpath}Number_of_viral_contigs_barplot_2018.png', dpi=400, figsize=(size1, size2))
    plt.savefig(f'{outpath}Number_of_viral_contigs_barplot_2018.svg', dpi=400, figsize=(size1, size2))   
            
            
    return

#compare_number_of_viral_contigs(Datasets_path, Output_path)




#Folder with clustered contigs.
Clust_path="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\Spongy\Viral_contigs\Viral_contigs_analysis\\Cdhit_3_all_full_name_2018\All_viral_contigs_cdhit_id95_c80_full_name_2018.clstr"


#######
#Read clustering results prepared by cd-hit, prepare dictionary of clusters.
#######

def read_cdhit_clusters_data(datapath):
    filein=open(datapath, 'r')
    Clusters_dict={}
    num_seqs=0
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
                num_seqs+=len(spacer_names)
            cluster_name=line.lstrip('>')  
            cluster_list=[]
            initiation=0
        else:
            line=line.split(' ')
            line[1]=line[1].lstrip('>').rstrip('...')
            cluster_list.append(line[1:3])
    
    #For the last cluster.
    cluster_list_sorted=sorted(cluster_list, key=lambda x: x[1])
    spacer_names=[]
    for spacer_name in cluster_list_sorted:
        spacer_names.append(spacer_name[0])
    Clusters_dict[cluster_name]=spacer_names 
    num_seqs+=len(spacer_names)
    
    print(f'Number of clusters: {len(Clusters_dict)}')   
    print(f'Number of viral contigs: {num_seqs}')
        
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
    plt.savefig(f'{output_path}2018_Cdhit_clusters_size_distribution{name}.png')    

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
            if Matrix_of_shared.at[keys_list[i], keys_list[j]]!=0:
                text = ax.text(j, i, Matrix_of_shared.at[keys_list[i], keys_list[j]], ha="center", va="center", color="k", size=7)
    
    fig.tight_layout()
    plt.show()
    plt.savefig(f'{outpath}2018_{name}_heatmap.png', dpi=400, figsize=(size, size))
    plt.savefig(f'{outpath}2018_{name}_heatmap.svg', dpi=400, figsize=(size, size))
    return

#######
#Draw venn diagram.
#######  

def draw_venn(Matrix_of_shared, keys_list, name, outpath, size):
    #Visualize data with heatmap.
    fig=plt.figure(figsize=(size,size), dpi=100)
    ax=fig.add_subplot(111)
    
    if 'not_clustered' in name:
        venn2(subsets=(Matrix_of_shared[0]-Matrix_of_shared[2][0], Matrix_of_shared[1]-Matrix_of_shared[2][1], Matrix_of_shared[2][0]+Matrix_of_shared[2][1]), set_labels=set(keys_list))
        venn2_circles(subsets=(Matrix_of_shared[0]-Matrix_of_shared[2][0], Matrix_of_shared[1]-Matrix_of_shared[2][1], Matrix_of_shared[2][0]+Matrix_of_shared[2][1]))    
    else:
        venn2(subsets=(Matrix_of_shared.at[keys_list[0], keys_list[0]]-Matrix_of_shared.at[keys_list[1], keys_list[0]], Matrix_of_shared.at[keys_list[1], keys_list[1]]-Matrix_of_shared.at[keys_list[1], keys_list[0]], Matrix_of_shared.at[keys_list[1], keys_list[0]]), set_labels=set(keys_list))
        venn2_circles(subsets=(Matrix_of_shared.at[keys_list[0], keys_list[0]]-Matrix_of_shared.at[keys_list[1], keys_list[0]], Matrix_of_shared.at[keys_list[1], keys_list[1]]-Matrix_of_shared.at[keys_list[1], keys_list[0]], Matrix_of_shared.at[keys_list[1], keys_list[0]]))                 
        
    fig.tight_layout()
    plt.show()
    plt.savefig(f'{outpath}2018_{name}_venn.png', dpi=400, figsize=(size, size))
    plt.savefig(f'{outpath}2018_{name}_venn.svg', dpi=400, figsize=(size, size))    
    return


#######
#For any number of datasets.
#Count overlap of datasets.
#######

def count_contig_set_overlap(Clusters_dict, outpath):
    
    #Grouping by sample.
    List_of_all_sets=['BCT_H_panicea_2018', 'BCT_H_sitiens_2018', 'BCT_I_palmata_2018', 'BCT_MW_2018', 'VRS_H_panicea_2018',  'VRS_H_sitiens_2018', 'VRS_I_palmata_2018', 'VRS_MW_2018']
    
    spacer_all_sets_overlap_matrix_non_clustered=pd.DataFrame(0, index=List_of_all_sets, columns=List_of_all_sets)
    spacer_all_sets_overlap_matrix_clustered=pd.DataFrame(0,     index=List_of_all_sets, columns=List_of_all_sets)
    
    for cluster_id, spacer_set in Clusters_dict.items():
        sample_id_list=[]
        for spacer_id in spacer_set:
            if '_NODE_' in spacer_id:
                sample_id=spacer_id.split('_NODE_')[0]
            elif '_CUTOFF_' in spacer_id:
                sample_id=spacer_id.split('_CUTOFF_')[0]
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
    List_of_loc_sets=['H_panicea', 'H_sitiens', 'I_palmata', 'MW']
    
    spacer_loc_sets_overlap_matrix_non_clustered=pd.DataFrame(0, index=List_of_loc_sets, columns=List_of_loc_sets)
    spacer_loc_sets_overlap_matrix_clustered=pd.DataFrame(0,     index=List_of_loc_sets, columns=List_of_loc_sets)
    
    for cluster_id, spacer_set in Clusters_dict.items():
        sample_id_list=[]
        for spacer_id in spacer_set:
            sample_id=spacer_id.split('_20')[0]
            sample_id=sample_id.lstrip('VRS_').lstrip('BTC_')
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
    print('Here')
    
    #Grouping by year.
    #List_of_year_sets=['2016', '2018']
    #
    #spacer_year_sets_overlap_matrix_non_clustered=pd.DataFrame(0, index=List_of_year_sets, columns=List_of_year_sets)
    #spacer_year_sets_overlap_matrix_non_clustered_venn=[0,0,[0,0]]
    #spacer_year_sets_overlap_matrix_clustered=pd.DataFrame(0,     index=List_of_year_sets, columns=List_of_year_sets)
    #
    #for cluster_id, spacer_set in Clusters_dict.items():
    #    sample_id_list=[]
    #    for spacer_id in spacer_set:
    #        if '_NODE_' in spacer_id:
    #            sample_id=spacer_id.split('_NODE_')[0].split('_')[-1]
    #        elif '_CUTOFF_' in spacer_id:
    #            sample_id=spacer_id.split('_CUTOFF_')[0].split('_')[-1]          
    #        
    #        sample_id_list.append(sample_id)
    #    sample_id_set=list(set(sample_id_list))  
    #    
    #    for i in range(len(sample_id_set)):
    #        #Different spacers counting (with clustering).
    #        spacer_year_sets_overlap_matrix_clustered.at[sample_id_set[i], sample_id_set[i]]+=1
    #        #All spacers counting (without clustering).
    #        spacer_year_sets_overlap_matrix_non_clustered.at[sample_id_set[i], sample_id_set[i]]+=sample_id_list.count(sample_id_set[i])
    #        spacer_year_sets_overlap_matrix_non_clustered_venn[List_of_year_sets.index(sample_id_set[i])]+=sample_id_list.count(sample_id_set[i])
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
    #                spacer_year_sets_overlap_matrix_non_clustered_venn[2][0]+=sample_id_list.count(sample_id_pair_sorted[1][0])
    #                spacer_year_sets_overlap_matrix_non_clustered_venn[2][1]+=sample_id_list.count(sample_id_pair_sorted[0][0])
    #                        
    #draw_heatmap(spacer_year_sets_overlap_matrix_clustered,          List_of_year_sets, 'Years_clustered',     outpath, 2)
    #draw_heatmap(spacer_year_sets_overlap_matrix_non_clustered,      List_of_year_sets, 'Years_not_clustered', outpath, 2) 
    #draw_venn(spacer_year_sets_overlap_matrix_clustered,             List_of_year_sets, 'Years_clustered',     outpath, 3)
    #draw_venn(spacer_year_sets_overlap_matrix_non_clustered_venn,    List_of_year_sets, 'Years_not_clustered', outpath, 3)
    
    #Grouping by bact/vir.
    List_of_BV_sets=['Bacterial', 'Viral']
    
    spacer_BV_sets_overlap_matrix_non_clustered=pd.DataFrame(0, index=List_of_BV_sets, columns=List_of_BV_sets)
    spacer_BV_sets_overlap_matrix_non_clustered_venn=[0,0,[0,0]]
    spacer_BV_sets_overlap_matrix_clustered=pd.DataFrame(0,     index=List_of_BV_sets, columns=List_of_BV_sets)
    
    for cluster_id, spacer_set in Clusters_dict.items():
        sample_id_list=[]
        for spacer_id in spacer_set:
            if 'VRS_' in spacer_id:
                sample_id_list.append('Viral')
            elif 'BCT_' in spacer_id:
                sample_id_list.append('Bacterial')
        sample_id_set=list(set(sample_id_list))  
        
        for i in range(len(sample_id_set)):
            #Different spacers counting (with clustering).
            spacer_BV_sets_overlap_matrix_clustered.at[sample_id_set[i], sample_id_set[i]]+=1
            #All spacers counting (without clustering).
            spacer_BV_sets_overlap_matrix_non_clustered.at[sample_id_set[i], sample_id_set[i]]+=sample_id_list.count(sample_id_set[i])
            spacer_BV_sets_overlap_matrix_non_clustered_venn[List_of_BV_sets.index(sample_id_set[i])]+=sample_id_list.count(sample_id_set[i])  
            
            for j in range(len(sample_id_set)):
                if j>i:
                    #Sort elements in the pair to match the order of elements in the matrix to prepare a triangle matrix.
                    sample_id_pair=[[sample_id_set[i], List_of_BV_sets.index(sample_id_set[i])], [sample_id_set[j], List_of_BV_sets.index(sample_id_set[j])]]
                    sample_id_pair_sorted=sorted(sample_id_pair, key=lambda x: x[1], reverse=True)
                    #Different spacers counting (with clustering).
                    spacer_BV_sets_overlap_matrix_clustered.at[sample_id_pair_sorted[0][0], sample_id_pair_sorted[1][0]]+=1
                    #All spacers counting (without clustering).
                    spacer_BV_sets_overlap_matrix_non_clustered.at[sample_id_pair_sorted[0][0], sample_id_pair_sorted[1][0]]+=(sample_id_list.count(sample_id_pair_sorted[0][0])+sample_id_list.count(sample_id_pair_sorted[1][0]))
                    spacer_BV_sets_overlap_matrix_non_clustered_venn[2][0]+=sample_id_list.count(sample_id_pair_sorted[1][0])
                    spacer_BV_sets_overlap_matrix_non_clustered_venn[2][1]+=sample_id_list.count(sample_id_pair_sorted[0][0])
                            
    draw_heatmap(spacer_BV_sets_overlap_matrix_clustered,          List_of_BV_sets, 'Bact_Vir_clustered',     outpath, 2)
    draw_heatmap(spacer_BV_sets_overlap_matrix_non_clustered,      List_of_BV_sets, 'Bact_Vir_not_clustered', outpath, 2)     
    draw_venn(spacer_BV_sets_overlap_matrix_clustered,             List_of_BV_sets, 'Bact_Vir_clustered',     outpath, 3)
    draw_venn(spacer_BV_sets_overlap_matrix_non_clustered_venn,    List_of_BV_sets, 'Bact_Vir_not_clustered', outpath, 3)    
                
    
    return 


Viral_contigs_clust=read_cdhit_clusters_data(Clust_path)
analyse_clusters_dict(Viral_contigs_clust, 'viral_contigs_cdhit_id95_c95', Output_path)
count_contig_set_overlap(Viral_contigs_clust, Output_path)



