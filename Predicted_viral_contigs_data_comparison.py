###############################################
##Dmitry Sutormin, 2023##
##Prediction of viral contigs data##

#1) Script compares the numbers of filtered predicted viral contigs across datasets.
#2) The script merges fasta files with viral contigs into one fasta file for clustering.
#3) Reads results of contigs clustering and draws heatmap representing the datasets overlappings.
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
import matplotlib.patches as mpatches
import csv
import pandas as pd
import seaborn as sns


#######
## Part 1. Count viral contigs and prepare merged file for clustering.
#######

#######
#Data to be used.
#######

#Folder with filtered datasets.
Datasets_path="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\Spongy\Viral_contigs\Viral_contigs_prediction_and_filtering\\"

#DeepVirFinder score threshold.
Score_thr=0.6

#DeepVirFinder p-value threshold.
Pvalue_thr=0.05

#Output folder.
Output_path="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\Spongy\Viral_contigs\Viral_contigs_prediction_and_filtering\\"


#######
#Fasta file parser. Makes contigs names shorter.
#######

def read_fasta_and_rename(inpath, set_name, all_contigs_dict):
    
    fasta_dict={}
    
    filein=open(inpath, 'r')
    i=1
    for record in SeqIO.parse(filein, "fasta"):
        record_name=f'{set_name}_{i}'
        record_name=record_name.replace('illumina', 'Il')
        record_name=record_name.replace('nanopore', 'Na')
        record_name=record_name.replace('sitiens', 'sit')
        record_name=record_name.replace('panicea', 'pan')
        record_name=record_name.replace('palmata', 'pal')
        record_seq=str(record.seq)
        
        fasta_dict[record_name]=record_seq

        i+=1
    
    all_contigs_dict[set_name]=fasta_dict
        
    filein.close()
    
    print(f'Number of viral contigs found for {set_name}: {len(fasta_dict)}')
    
    return all_contigs_dict

#######
#Make barplot of the number of viral contigs.
#######

def barplot_viral_contigs(All_contigs_dict, Color_dict, Barplot_name):
    
    fig, ax = plt.subplots()
    
    j=0
    x_coords=[]
    Num_cumul=0
    Set_name_base_dict={}
    Set_name_base_list=[]
    for Set_name, Contig_dict in All_contigs_dict.items():
        
        Set_name_base=Set_name.rstrip('_illumina').rstrip('_nanopore')
        Seq_type=Set_name.split('_')[-1]
        Num_cont=len(Contig_dict)
        
        if Set_name_base not in Set_name_base_dict:
            Num_cumul=0
            Bar_coord=j
            if Bar_coord==0:
                ax.bar(Bar_coord, Num_cont, 0.65, color=Color_dict[Seq_type], bottom=Num_cumul, tick_label=Set_name_base, label=f'{Seq_type}')
            else:
                ax.bar(Bar_coord, Num_cont, 0.65, color=Color_dict[Seq_type], bottom=Num_cumul, tick_label=Set_name_base)  
            
            Set_name_base_dict[Set_name_base]=[Bar_coord, Num_cont] 
            Set_name_base_list.append(Set_name_base)
            x_coords.append(j)
            j+=1
        
        else:
            Num_cumul=Set_name_base_dict[Set_name_base][1]
            Bar_coord=Set_name_base_dict[Set_name_base][0]
            if Bar_coord==1:
                ax.bar(Bar_coord, Num_cont, 0.65, color=Color_dict[Seq_type], bottom=Num_cumul, tick_label=Set_name_base, label=f'{Seq_type}')
            else:
                ax.bar(Bar_coord, Num_cont, 0.65, color=Color_dict[Seq_type], bottom=Num_cumul, tick_label=Set_name_base)  
            

    #ax.set_xlim([-0.5,0.5])
    ax.tick_params(axis='x', which='major', labelsize=12)
    ax.set_xticks(x_coords, labels=Set_name_base_list, minor=False, size=12, rotation=90)
    ax.set_ylabel('Number of contigs', size=12)
    ax.set_xlabel('Metagenomes', size=12)
    ax.set_title('Before clustering', size=15)
    ax.spines['top'].set_color(None)
    ax.spines['right'].set_color(None)
    
    #Place legend outside of a graph. Taken from: https://stackoverflow.com/questions/4700614/how-to-put-the-legend-out-of-the-plot
    box=ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    ax.legend(fontsize=12, ncol=1, loc='upper left', frameon=False, markerscale=1, handlelength=0.7, handletextpad=0.3, columnspacing=0.7, bbox_to_anchor=(1, 0.95))
    plt.tight_layout(rect=[0,0,0.7,1])
    plt.show()
    plt.savefig(f'{Barplot_name}.png', dpi=300, figsize=(2, 5))
    plt.savefig(f'{Barplot_name}.svg', dpi=300, figsize=(2, 5))    
    
    return


def compare_number_of_viral_contigs(inpath, score_thr, pvalue_thr, outpath):
    
    #Read input fasta files.
    All_contigs_dict={}
    Contigs_num_dict={}
    
    for year_dir in os.listdir(inpath):
        if os.path.isdir(f'{inpath}\{year_dir}'):
            #print(year_dir)
            for sample_dir in os.listdir(f'{inpath}\{year_dir}'):
                if os.path.isdir(f'{inpath}\\{year_dir}\\{sample_dir}'):
                    #print(sample_dir)  
                    Set_name=sample_dir
                    Input_fasta_path=f'{inpath}\\{year_dir}\\{sample_dir}\\{sample_dir}_sc_{score_thr}_pval_{pvalue_thr}_vir_contig_set_final.fasta'
                    All_contigs_dict=read_fasta_and_rename(Input_fasta_path, Set_name, All_contigs_dict)
                    Contigs_num_dict[Set_name]=len(All_contigs_dict[Set_name])
                    
    
    #Write merged fasta file.
    fasta_fileout=open(f'{outpath}\\All_viral_contigs_sc_{score_thr}_pval_{pvalue_thr}_unclust.fasta', 'w')
    for set_name, set_cont_dict in All_contigs_dict.items():
        for contig_name, contig_seq in set_cont_dict.items():
            fasta_fileout.write(f'>{contig_name}\n{contig_seq}\n')
    fasta_fileout.close()
    
    #Count number of contigs and visualize with a barplot.
    Barplot_name=f'{outpath}\\All_viral_contigs_sc_{score_thr}_pval_{pvalue_thr}_unclust'
    Color_dict={'illumina' : '#3cd7eb', 'nanopore' : '#ed9892'}
    barplot_viral_contigs(All_contigs_dict, Color_dict, Barplot_name)
            
    return

#compare_number_of_viral_contigs(Datasets_path, Score_thr, Pvalue_thr, Output_path)




#######
## Part 2. Reads and analyses results of clustering.
#######

#######
#Data to be used.
#######

#Dataset name.
Dataset_name="All_viral_contigs_2016_2018_2022"

#Clustering identity level.
Clust_identity_level='0.90'

#Fasta with all contigs before clustering.
All_contigs_fasta="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\Spongy\Viral_contigs\Viral_contigs_analysis\All_viral_contigs_clustering\\All_viral_contigs_sc_0.6_pval_0.05_unclust.fasta"

#Folder with clustered contigs.
Clust_res_path=f'C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\Spongy\Viral_contigs\Viral_contigs_analysis\All_viral_contigs_clustering\\All_viral_contigs_sc_0.6_pval_0.05_gclust_{Clust_identity_level}.txt'

#Output folder.
Output_path="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\Spongy\Viral_contigs\Viral_contigs_analysis\All_viral_contigs_clustering\\"


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
#Read clustering results prepared by gclust, prepare dictionary of clusters.
#######

def read_gclust_clusters_data(datapath):
    filein=open(datapath, 'r')
    Clusters_dict={}
    num_seqs=0
    initiation=1
    for line in filein:
        line=line.rstrip('\r\n')
        if line[0]==">":
            if initiation!=1:
                Clusters_dict[cluster_name]=cluster_list 
                num_seqs+=len(cluster_list)
            cluster_name=line.lstrip('>')  
            cluster_list=[]
            initiation=0
        elif line[0:3]!="...":
            line=line.split(' ')
            contig_name=line[1].lstrip('>').rstrip('')
            contig_len=int(line[0].split('\t')[1].rstrip('nt,'))
            cluster_list.append([contig_name, contig_len])
    
    #For the last cluster.
    Clusters_dict[cluster_name]=cluster_list 
    num_seqs+=len(cluster_list)
     
    filein.close()
    
    print(f'Number of clusters: {len(Clusters_dict)}')   
    print(f'Number of viral contigs: {num_seqs}')    
    
    return Clusters_dict

#######
#Simple fasta file parser.
#######

def read_fasta(inpath):
    
    fasta_dict={}
    
    filein=open(inpath, 'r')
    for record in SeqIO.parse(filein, "fasta"):
        record_name=str(record.id)
        record_seq=str(record.seq)
        fasta_dict[record_name]=record_seq
        
    filein.close()
    
    print(f'Initial number of viral contigs (unclustered): {len(fasta_dict)}')
    
    return fasta_dict

#######
#Analyse dictionary of clusters: clusters size distribution.
#######

def analyse_clusters_size(Clusters_dict, dataset_name, clust_identity_level, output_path):
    
    Clusters_len=[]
    for cluster_name, cluster_dep_ar in Clusters_dict.items():
        Clusters_len.append(len(cluster_dep_ar))
        #if len(cluster_dep_ar)>2:
        #    print(f'{cluster_name} : {cluster_dep_ar}')
    
    #Plot distribution of clusters size.
    plt.hist(Clusters_len, rwidth=0.95, lw=1, ec="black")
    plt.yscale('log')
    plt.xticks(fontsize=13)
    plt.yticks(fontsize=13)
    plt.title(f'Size of clusters distribution after gclust with id {int(float(clust_identity_level)*100)}%', size=13)
    plt.xlabel(f'Size of clusters', size=13)
    plt.ylabel('Number of clusters, log', size=13)
    plt.show()
    plt.savefig(f'{output_path}{dataset_name}_clusters_size_distribution_gclust_id_{clust_identity_level}.png')    
    plt.savefig(f'{output_path}{dataset_name}_clusters_size_distribution_gclust_id_{clust_identity_level}.svg')

    return

#######
#Select representative contig - the longest one.
####### 

def select_repres_cont(Clusters_dict, sample_set_name):
    
    Clusters_dict_upd={}
    
    for cluster_name, cluster_list in Clusters_dict.items():
        longest_cont=0
        if len(cluster_list)==1:
            Clusters_dict_upd[cluster_name]=cluster_list
        else:
            cluster_list_sorted=sorted(cluster_list, key=lambda x: x[1], reverse=True)
            #print(cluster_list_sorted)
            cluster_list_upd=[]
            for i in range(len(cluster_list_sorted)):
                if sample_set_name not in cluster_list_sorted[i][0]:
                    cluster_list_upd.append(cluster_list_sorted[i])
                elif (sample_set_name in cluster_list_sorted[i][0]) and (longest_cont==0):
                    cluster_list_upd.append(cluster_list_sorted[i])
                    longest_cont=1
                elif (sample_set_name in cluster_list_sorted[i][0]) and (longest_cont==1):
                    continue
            #print(cluster_list_upd)
            Clusters_dict_upd[cluster_name]=cluster_list_upd
                
    
    return Clusters_dict_upd

#######
#Make a non-redundant sets of contigs within sample types (a particular sponge or a marine water).
#E.g., merge illumina with nanopore and VRS with BCT.
#######

def merg_meth_vars(Clusters_dict, sample_set_list):
    
    Clusters_dict_upd={}
    
    if len(sample_set_list)==1:
        Clusters_dict_upd=Clusters_dict
    
    else: 
        for cluster_name, cluster_list in Clusters_dict.items():
            #print(cluster_name, cluster_list)
            if len(cluster_list)==1:
                Clusters_dict_upd[cluster_name]=cluster_list   
            else:
                cluster_list_sorted=sorted(cluster_list, key=lambda x: x[1], reverse=True)
                cluster_list_upd=[]
                sample_set_name_base=sample_set_list[0].rstrip('_Il').rstrip('_Na').lstrip('BCT_').lstrip('VRS_')
                #Merge BCT and VRS or Illumina and Nanopore contigs.
                tmp_list=[]
                for i in range(len(cluster_list_sorted)):
                    if sample_set_name_base not in cluster_list_sorted[i][0]:
                        cluster_list_upd.append(cluster_list_sorted[i])
                    elif sample_set_name_base in cluster_list_sorted[i][0]:
                        tmp_list.append(cluster_list_sorted[i])
                if len(tmp_list)==1:
                    cluster_list_upd.append(tmp_list[0])
                elif len(tmp_list)>1:
                    Selected_contig=tmp_list[0]
                    Selected_contig_upd=Selected_contig+['both']
                    cluster_list_upd.append(Selected_contig_upd)

                        
                Clusters_dict_upd[cluster_name]=cluster_list_upd                

    return Clusters_dict_upd


#######
#Make barplot of the number of viral contigs for all samples.
#######

def barplot_viral_contigs_all(Clusters_dict, List_of_all_dict, Color_dict, Barplot_name, clust_identity_level):
    
    #Calculate numbers of contigs.
    Contigs_num_dict={}
    for sample_year, sample_dict in List_of_all_dict.items():
        for sample_type, sample_name_list in sample_dict.items():
            
            sample_set_name_base=sample_name_list[0].rstrip('_Il').rstrip('_Na').lstrip('BCT_').lstrip('VRS_')
            if ('VRS' in sample_name_list[0]) or ('BCT' in sample_name_list[0]):
                Contigs_num_dict[f'{sample_type} {sample_year}']={'Bac fraction contigs' : 0, 'In both fractions (Bac+Vir)' : 0, 'Vir fraction contigs' : 0}
                for sample_name in sample_name_list:  
                    for cluster_name, cluster_ar in Clusters_dict.items():
                        for cluster_data in cluster_ar:
                            if (sample_name in cluster_data[0]) and (len(cluster_data)>2):
                                Contigs_num_dict[f'{sample_type} {sample_year}']['In both fractions (Bac+Vir)']+=1
                            elif (sample_name in cluster_data[0]) and (len(cluster_data)<=2) and ('VRS_' in sample_name):
                                Contigs_num_dict[f'{sample_type} {sample_year}']['Vir fraction contigs']+=1
                            elif (sample_name in cluster_data[0]) and (len(cluster_data)<=2) and ('BCT_' in sample_name):
                                Contigs_num_dict[f'{sample_type} {sample_year}']['Bac fraction contigs']+=1                                                
            else:
                Contigs_num_dict[f'{sample_type} {sample_year}']={'Illum contigs' : 0, 'In both fractions (Illum+Nano)' : 0, 'Nano contigs' : 0}
                for sample_name in sample_name_list:  
                    for cluster_name, cluster_ar in Clusters_dict.items():
                        for cluster_data in cluster_ar:
                            if (sample_name in cluster_data[0]) and (len(cluster_data)>2):
                                Contigs_num_dict[f'{sample_type} {sample_year}']['In both fractions (Illum+Nano)']+=1
                            elif (sample_name in cluster_data[0]) and (len(cluster_data)<=2) and ('_Il' in sample_name):
                                Contigs_num_dict[f'{sample_type} {sample_year}']['Illum contigs']+=1
                            elif (sample_name in cluster_data[0]) and (len(cluster_data)<=2) and ('_Na' in sample_name):
                                Contigs_num_dict[f'{sample_type} {sample_year}']['Nano contigs']+=1                
            
    print(Contigs_num_dict)                
            
    
    fig, ax = plt.subplots()
    
    j=0
    x_coords=[]
    Set_name_base_list=[]
    for Set_name, Subset_cont_num_dict in Contigs_num_dict.items():
        
        Num_cumul=0
        for Subset_name, Cont_num in Subset_cont_num_dict.items():
            
            Bar_coord=j
            if Bar_coord==0:
                ax.bar(Bar_coord, Cont_num, 0.65, color=Color_dict[Subset_name], bottom=Num_cumul, tick_label=Set_name, label=f'{Subset_name}')
            else:
                ax.bar(Bar_coord, Cont_num, 0.65, color=Color_dict[Subset_name], bottom=Num_cumul, tick_label=Set_name)  
            Num_cumul+=Cont_num
            
        Set_name_base_list.append(Set_name)
        x_coords.append(j)
        j+=1        
    
    #ax.set_xlim([-0.5,0.5])
    ax.tick_params(axis='x', which='major', labelsize=12)
    ax.set_xticks(x_coords, labels=Set_name_base_list, minor=False, size=12, rotation=90)
    ax.set_ylabel('Number of contigs', size=12)
    ax.set_xlabel('Metagenomes', size=12)
    ax.set_title(f'After gclust with id {clust_identity_level}', size=14)
    ax.spines['top'].set_color(None)
    ax.spines['right'].set_color(None)
    
    #Create and place legend outside of a graph. Taken from: https://stackoverflow.com/questions/4700614/how-to-put-the-legend-out-of-the-plot
    box=ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    handles_ar=[]
    for Subset_name, subset_color in Color_dict.items():
        legend_patch=mpatches.Patch(color=subset_color, label=Subset_name)
        handles_ar.append(legend_patch)   
    ax.legend(handles=handles_ar, fontsize=10, ncol=1, loc='upper left', frameon=False, markerscale=1, handlelength=0.7, handletextpad=0.3, columnspacing=0.7, bbox_to_anchor=(1, 0.95))
    plt.tight_layout(rect=[0,0,0.8,1])
    plt.show()
    plt.savefig(f'{Barplot_name}.png', dpi=300, figsize=(5, 5))
    plt.savefig(f'{Barplot_name}.svg', dpi=300, figsize=(5, 5))    
    
    return

#######
#Count contigs shared between samples.
#######

def count_shared_cont(Clusters_dict, List_of_all_sets, Cont_num_dataframe):
    
    for clust_name, clust_ar in Clusters_dict.items():
        #Filling the dataframe.
        for i in range(len(clust_ar)):
            #Filling the diagonal.
            contig_name_1=clust_ar[i][0]
            contig_name_1=contig_name_1.lstrip('BCT_').lstrip('VRS_')
            if '_Il_' in contig_name_1:
                contig_name_base_1=contig_name_1.split('_Il_')[0]
            elif '_Na_' in contig_name_1:
                contig_name_base_1=contig_name_1.split('_Na_')[0]
            Cont_num_dataframe.at[contig_name_base_1, contig_name_base_1]+=1
            for j in range(len(clust_ar)):
                #Filling non-diagonal elements.
                if j>i:
                    contig_name_2=clust_ar[j][0]
                    contig_name_2=contig_name_2.lstrip('BCT_').lstrip('VRS_')
                    if '_Il_' in contig_name_2:
                        contig_name_base_2=contig_name_2.split('_Il_')[0]
                    elif '_Na_' in contig_name_2:
                        contig_name_base_2=contig_name_2.split('_Na_')[0]                    
                    #Make upper triangle matrix.
                    Index_1=List_of_all_sets.index(contig_name_base_1)
                    Index_2=List_of_all_sets.index(contig_name_base_2)
                    if Index_1>Index_2:
                        Cont_num_dataframe.at[contig_name_base_2, contig_name_base_1]+=1 
                    else:
                        Cont_num_dataframe.at[contig_name_base_1, contig_name_base_2]+=1                    
                    
    print(Cont_num_dataframe)
    
    return Cont_num_dataframe

#######
#Draw heatmap.
#######  

def draw_heatmap(Matrix_of_shared, keys_list, clust_identity_level, outpath, size):
    #Visualize data with heatmap.
    fig=plt.figure(figsize=(size,size), dpi=100)
    ax=fig.add_subplot(111)
    #Based on: https://matplotlib.org/gallery/images_contours_and_fields/image_annotated_heatmap.html#sphx-glr-gallery-images-contours-and-fields-image-annotated-heatmap-py
    Matrix_of_shared_log=np.log(Matrix_of_shared)
    im=ax.imshow(Matrix_of_shared_log, cmap='Wistia') #gist_yarg
    ax.set_xticks(np.arange(len(keys_list)))
    ax.set_yticks(np.arange(len(keys_list)))
    ax.set_xticklabels(keys_list, size=10, rotation=90)
    ax.set_yticklabels(keys_list, size=10)  
    ax.set_ylim(sorted(ax.get_xlim(), reverse=True)) #Solves a bug in matplotlib 3.1.1 discussed here: https://stackoverflow.com/questions/56942670/matplotlib-seaborn-first-and-last-row-cut-in-half-of-heatmap-plot
    
    # Minor ticks taken from https://stackoverflow.com/questions/38973868/adjusting-gridlines-and-ticks-in-matplotlib-imshow
    ax.set_xticks(np.arange(-.5, len(keys_list), 1), minor=True)
    ax.set_yticks(np.arange(-.5, len(keys_list), 1), minor=True)
    ax.tick_params(axis='both', which='minor', length=0)
    
    #Set the heatmap title
    ax.set_title(f'Shared contigs after gclust with id {clust_identity_level}', size=12)
    
    # Gridlines based on minor ticks
    ax.grid(which='minor', color='grey', linestyle='-', linewidth=1)    
    
    for i in range(len(keys_list)):
        for j in range(len(keys_list)):
            if Matrix_of_shared.at[keys_list[i], keys_list[j]]!=0:
                text = ax.text(j, i, Matrix_of_shared.at[keys_list[i], keys_list[j]], ha="center", va="center", color="k", size=7)
    
    fig.tight_layout()
    plt.show()
    plt.savefig(f'{outpath}All_rrepresentative_viral_contigs_gclust_id_{clust_identity_level}_heatmap.png', dpi=400, figsize=(size, size))
    plt.savefig(f'{outpath}All_rrepresentative_viral_contigs_gclust_id_{clust_identity_level}_heatmap.svg', dpi=400, figsize=(size, size))
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

def count_contig_set_overlap(Clusters_dict, clust_identity_level, outpath):
    
    #Grouping of samples.
    List_of_all_dict={'2016' : {'H_panicea' : ['H_pan_2016_Il'], 
                                'H_sitiens' : ['H_sit_2016_Il', 'H_sit_2016_Na'], 
                                'I_palmata' : ['I_pal_2016_Il', 'I_pal_2016_Na'], 
                                'M_water'   : ['MW_2016_Il']},
                      '2018' : {'H_panicea' : ['BCT_H_pan_2018_Il', 'VRS_H_pan_2018_Il'], 
                                'H_sitiens' : ['BCT_H_sit_2018_Il', 'VRS_H_sit_2018_Il'], 
                                'I_palmata' : ['BCT_I_pal_2018_Il', 'VRS_I_pal_2018_Il'], 
                                'M_water'   : ['BCT_MW_2018_Il', 'VRS_MW_2018_Il']},
                      '2022' : {'H_panicea' : ['H_pan_2022_Il', 'H_pan_2022_Na'], 
                                'H_sitiens' : ['H_sit_2022_Il', 'H_sit_2022_Na'], 
                                'I_palmata' : ['I_pal_2022_Il', 'I_pal_2022_Na'], 
                                'M_water'   : ['MW_2022_Il']},                      
                      }
    
    #Make a non-redundant sets of contigs within samples.
    for sampling_year, year_dict in List_of_all_dict.items():
        for sample_name, sample_set_list in year_dict.items():
            for sample_set_name in sample_set_list:
                #Select representative contig within a set - the longest one.
                Clusters_dict=select_repres_cont(Clusters_dict, sample_set_name)
                
    #Make a non-redundant sets of contigs within sample types (a particular sponge or a marine water).
    #E.g., merge illumina with nanopore and VRS with BCT.
    for sampling_year, year_dict in List_of_all_dict.items():
        for sample_name, sample_set_list in year_dict.items():  
            Clusters_dict=merg_meth_vars(Clusters_dict, sample_set_list)
            
    #Count numbers of contigs and visualize them with a barplot.
    Barplot_name=f'{outpath}\\All_viral_contigs_gclust_id_{clust_identity_level}'
    Color_dict={'Illum contigs' : '#3cd7eb', 'In both fractions (Illum+Nano)': '#e3f0a5', 'Nano contigs' : '#ed9892',
                'Bac fraction contigs': '#c390ed', 'In both fractions (Bac+Vir)': '#8a7ff0', 'Vir fraction contigs' : '#77f2ea'}
    barplot_viral_contigs_all(Clusters_dict, List_of_all_dict, Color_dict, Barplot_name, clust_identity_level)    
    
    #Create a mock dataframe to store the numbers of contigs shared between datasets.
    List_of_all_sets=[]
    for sampling_year, year_dict in List_of_all_dict.items():
        for sample_name, sample_set_list in year_dict.items():
            for sample_set_name in sample_set_list:
                sample_set_name_base=sample_set_name.rstrip('_Il').rstrip('_Na').lstrip('BCT_').lstrip('VRS_')
                if sample_set_name_base not in List_of_all_sets:
                    List_of_all_sets.append(sample_set_name_base)
    
    Cont_num_dataframe=pd.DataFrame(0, index=List_of_all_sets, columns=List_of_all_sets)
    print(Cont_num_dataframe)
    
    #Count contigs shared between samples.
    Cont_num_dataframe_filled=count_shared_cont(Clusters_dict, List_of_all_sets, Cont_num_dataframe)
    
    #Draw heatmap.
    draw_heatmap(Cont_num_dataframe_filled, List_of_all_sets, clust_identity_level, outpath, 5)
    
    
    
    #spacer_all_sets_overlap_matrix_non_clustered=pd.DataFrame(0, index=List_of_all_sets, columns=List_of_all_sets)
    #spacer_all_sets_overlap_matrix_clustered=pd.DataFrame(0,     index=List_of_all_sets, columns=List_of_all_sets)
    #
    #for cluster_id, spacer_set in Clusters_dict.items():
    #    sample_id_list=[]
    #    for spacer_id in spacer_set:
    #        if '_NODE_' in spacer_id:
    #            sample_id=spacer_id.split('_NODE_')[0]
    #        elif '_CUTOFF_' in spacer_id:
    #            sample_id=spacer_id.split('_CUTOFF_')[0]
    #        sample_id_list.append(sample_id)
    #    sample_id_set=list(set(sample_id_list))  
    #    
    #    for i in range(len(sample_id_set)):
    #        #Different spacers counting (with clustering).
    #        spacer_all_sets_overlap_matrix_clustered.at[sample_id_set[i], sample_id_set[i]]+=1
    #        #All spacers counting (without clustering).
    #        spacer_all_sets_overlap_matrix_non_clustered.at[sample_id_set[i], sample_id_set[i]]+=sample_id_list.count(sample_id_set[i])
    #        
    #        for j in range(len(sample_id_set)):
    #            if j>i:
    #                #Sort elements in the pair to match the order of elements in the matrix to prepare a triangle matrix.
    #                sample_id_pair=[[sample_id_set[i], List_of_all_sets.index(sample_id_set[i])], [sample_id_set[j], List_of_all_sets.index(sample_id_set[j])]]
    #                sample_id_pair_sorted=sorted(sample_id_pair, key=lambda x: x[1], reverse=True)
    #                #Different spacers counting (with clustering).
    #                spacer_all_sets_overlap_matrix_clustered.at[sample_id_pair_sorted[0][0], sample_id_pair_sorted[1][0]]+=1
    #                #All spacers counting (without clustering).
    #                spacer_all_sets_overlap_matrix_non_clustered.at[sample_id_pair_sorted[0][0], sample_id_pair_sorted[1][0]]+=(sample_id_list.count(sample_id_pair_sorted[0][0])+sample_id_list.count(sample_id_pair_sorted[1][0]))
    #                        
    #draw_heatmap(spacer_all_sets_overlap_matrix_clustered,     List_of_all_sets, 'All_datasets_clustered',     outpath, 5)
    #draw_heatmap(spacer_all_sets_overlap_matrix_non_clustered, List_of_all_sets, 'All_datasets_not_clustered', outpath, 5)
    #    
    #    
    ##Grouping by location.
    #List_of_loc_sets=['H_panicea', 'H_sitiens', 'I_palmata', 'MW']
    #
    #spacer_loc_sets_overlap_matrix_non_clustered=pd.DataFrame(0, index=List_of_loc_sets, columns=List_of_loc_sets)
    #spacer_loc_sets_overlap_matrix_clustered=pd.DataFrame(0,     index=List_of_loc_sets, columns=List_of_loc_sets)
    #
    #for cluster_id, spacer_set in Clusters_dict.items():
    #    sample_id_list=[]
    #    for spacer_id in spacer_set:
    #        sample_id=spacer_id.split('_20')[0]
    #        sample_id=sample_id.lstrip('VRS_').lstrip('BTC_')
    #        sample_id_list.append(sample_id)
    #    sample_id_set=list(set(sample_id_list))  
    #    
    #    for i in range(len(sample_id_set)):
    #        #Different spacers counting (with clustering).
    #        spacer_loc_sets_overlap_matrix_clustered.at[sample_id_set[i], sample_id_set[i]]+=1
    #        #All spacers counting (without clustering).
    #        spacer_loc_sets_overlap_matrix_non_clustered.at[sample_id_set[i], sample_id_set[i]]+=sample_id_list.count(sample_id_set[i])
    #        
    #        for j in range(len(sample_id_set)):
    #            if j>i:
    #                #Sort elements in the pair to match the order of elements in the matrix to prepare a triangle matrix.
    #                sample_id_pair=[[sample_id_set[i], List_of_loc_sets.index(sample_id_set[i])], [sample_id_set[j], List_of_loc_sets.index(sample_id_set[j])]]
    #                sample_id_pair_sorted=sorted(sample_id_pair, key=lambda x: x[1], reverse=True)
    #                #Different spacers counting (with clustering).
    #                spacer_loc_sets_overlap_matrix_clustered.at[sample_id_pair_sorted[0][0], sample_id_pair_sorted[1][0]]+=1
    #                #All spacers counting (without clustering).
    #                spacer_loc_sets_overlap_matrix_non_clustered.at[sample_id_pair_sorted[0][0], sample_id_pair_sorted[1][0]]+=(sample_id_list.count(sample_id_pair_sorted[0][0])+sample_id_list.count(sample_id_pair_sorted[1][0]))
    #                        
    #draw_heatmap(spacer_loc_sets_overlap_matrix_clustered,     List_of_loc_sets, 'Locations_clustered',     outpath, 3)
    #draw_heatmap(spacer_loc_sets_overlap_matrix_non_clustered, List_of_loc_sets, 'Locations_not_clustered', outpath, 3)
    #print('Here')
    #
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
    #List_of_BV_sets=['Bacterial', 'Viral']
    #
    #spacer_BV_sets_overlap_matrix_non_clustered=pd.DataFrame(0, index=List_of_BV_sets, columns=List_of_BV_sets)
    #spacer_BV_sets_overlap_matrix_non_clustered_venn=[0,0,[0,0]]
    #spacer_BV_sets_overlap_matrix_clustered=pd.DataFrame(0,     index=List_of_BV_sets, columns=List_of_BV_sets)
    #
    #for cluster_id, spacer_set in Clusters_dict.items():
    #    sample_id_list=[]
    #    for spacer_id in spacer_set:
    #        if 'VRS_' in spacer_id:
    #            sample_id_list.append('Viral')
    #        elif 'BCT_' in spacer_id:
    #            sample_id_list.append('Bacterial')
    #    sample_id_set=list(set(sample_id_list))  
    #    
    #    for i in range(len(sample_id_set)):
    #        #Different spacers counting (with clustering).
    #        spacer_BV_sets_overlap_matrix_clustered.at[sample_id_set[i], sample_id_set[i]]+=1
    #        #All spacers counting (without clustering).
    #        spacer_BV_sets_overlap_matrix_non_clustered.at[sample_id_set[i], sample_id_set[i]]+=sample_id_list.count(sample_id_set[i])
    #        spacer_BV_sets_overlap_matrix_non_clustered_venn[List_of_BV_sets.index(sample_id_set[i])]+=sample_id_list.count(sample_id_set[i])  
    #        
    #        for j in range(len(sample_id_set)):
    #            if j>i:
    #                #Sort elements in the pair to match the order of elements in the matrix to prepare a triangle matrix.
    #                sample_id_pair=[[sample_id_set[i], List_of_BV_sets.index(sample_id_set[i])], [sample_id_set[j], List_of_BV_sets.index(sample_id_set[j])]]
    #                sample_id_pair_sorted=sorted(sample_id_pair, key=lambda x: x[1], reverse=True)
    #                #Different spacers counting (with clustering).
    #                spacer_BV_sets_overlap_matrix_clustered.at[sample_id_pair_sorted[0][0], sample_id_pair_sorted[1][0]]+=1
    #                #All spacers counting (without clustering).
    #                spacer_BV_sets_overlap_matrix_non_clustered.at[sample_id_pair_sorted[0][0], sample_id_pair_sorted[1][0]]+=(sample_id_list.count(sample_id_pair_sorted[0][0])+sample_id_list.count(sample_id_pair_sorted[1][0]))
    #                spacer_BV_sets_overlap_matrix_non_clustered_venn[2][0]+=sample_id_list.count(sample_id_pair_sorted[1][0])
    #                spacer_BV_sets_overlap_matrix_non_clustered_venn[2][1]+=sample_id_list.count(sample_id_pair_sorted[0][0])
    #                        
    #draw_heatmap(spacer_BV_sets_overlap_matrix_clustered,          List_of_BV_sets, 'Bact_Vir_clustered',     outpath, 2)
    #draw_heatmap(spacer_BV_sets_overlap_matrix_non_clustered,      List_of_BV_sets, 'Bact_Vir_not_clustered', outpath, 2)     
    #draw_venn(spacer_BV_sets_overlap_matrix_clustered,             List_of_BV_sets, 'Bact_Vir_clustered',     outpath, 3)
    #draw_venn(spacer_BV_sets_overlap_matrix_non_clustered_venn,    List_of_BV_sets, 'Bact_Vir_not_clustered', outpath, 3)    
                
    
    return Clusters_dict

#######
#Write representative contigs after clustering.
####### 

def write_representative_contigs(Clusters_data_dict_repr, All_contigs_dict, output_path, clust_identity_level):
    
    fileout_rep=open(f'{output_path}All_representative_viral_contigs_gclust_id_{clust_identity_level}.fasta', 'w')
    fileout_nr=open(f'{output_path}All_non_redundant_viral_contigs_gclust_id_{clust_identity_level}.fasta', 'w')
    Fin_num_cont_rep=1
    Fin_num_cont_nr=1
    for cluster_name, cluster_ar in Clusters_data_dict_repr.items():
        for contig_data in cluster_ar:
            fileout_rep.write(f'>{contig_data[0]}\n{All_contigs_dict[contig_data[0]]}\n')
            Fin_num_cont_rep+=1
        fileout_nr.write(f'>{cluster_ar[0][0]}\n{All_contigs_dict[cluster_ar[0][0]]}\n')
        Fin_num_cont_nr+=1
            
    fileout_rep.close()
    fileout_nr.close()
    
    print(f'Final number of representative viral contigs (clustered with gclust with identity level {clust_identity_level}): {Fin_num_cont_rep}')
    print(f'Final number of non-redundant viral contigs (clustered with gclust with identity level {clust_identity_level}): {Fin_num_cont_nr}')
    
    return

#######
#Wrapper script for clustering data analysis.
#######  

def clustering_wrapper(clust_res_path, dataset_name, clust_identity_level, all_contigs_fasta, output_path):
    
    #Read output of clustering.
    Clust_data_dict=read_gclust_clusters_data(clust_res_path)
    
    #Analyse cluster size distribution.
    analyse_clusters_size(Clust_data_dict, dataset_name, clust_identity_level, output_path)
    
    #Get representative sequences and analyse clusters overlaps.
    Clusters_data_dict_repr=count_contig_set_overlap(Clust_data_dict, clust_identity_level, output_path)
    
    #Return representative contigs.
    All_contigs_dict=read_fasta(all_contigs_fasta)
    write_representative_contigs(Clusters_data_dict_repr, All_contigs_dict, output_path, clust_identity_level)    
    
    return

clustering_wrapper(Clust_res_path, Dataset_name, Clust_identity_level, All_contigs_fasta, Output_path)




