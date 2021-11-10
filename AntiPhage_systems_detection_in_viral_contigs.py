###############################################
##Dmitry Sutormin, 2021##
##Anti phage defense systems detection in viral contigs##

# 1) Script takes protein sequences from PADS Arsenal database and groups them by homology.
# 2) Analyses the output of hmm search performed in metagenomic ORFs with hmm profiles prepared for antiphage genes from PADS Arsenal.
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
import copy



###############
###############
##
## Part 1.
## Detect and count defense systems.
##
###############
###############

#######
#Data to be used.
#######

#PWD.
PWD="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\Spongy\Anti_phage_systems\\"

E_value_thr=10e-10


#######
#Get contigs coverage depth.
#######

def read_cov_data_file(pwd):
    
    pwd_in=pwd+"Defense_systems\Contigs_coverage_5000\\"
    
    #List of coverage depth files.
    input_cov_list=os.listdir(pwd_in)
    
    Sample_cont_cov_dict={}
    Sample_assembly_len={}
    Sample_read_num={}
    
    for filename in input_cov_list:
        
        Sample_name=filename.split('_bowtie_')[0]
        
        Sample_cont_cov_dict[Sample_name]={}
        Sample_read_num[Sample_name]=0
        
        filein=open(pwd_in+filename, 'r')
        for line in filein:
            if line[0]!='#':
                line=line.rstrip().split('\t')
                contig_name=line[0]
                contig_length=int(line[2])
                num_reads=int(line[3])
                cov_depth=float(line[6])
                
                #Store all data.
                Sample_cont_cov_dict[Sample_name][contig_name]=[contig_length, num_reads, cov_depth]
                
                #Store accumulative len data for step plot.
                if Sample_name not in Sample_assembly_len:
                    Sample_assembly_len[Sample_name]=[contig_length]
                else:
                    Sample_assembly_len[Sample_name].append(Sample_assembly_len[Sample_name][-1]+contig_length)
                               
                #Count total number of mapped reads.
                Sample_read_num[Sample_name]+=num_reads
        
    #Plot cumulative assembly length curve.
    plot_cumul_step_curve(Sample_assembly_len, pwd)
    
    return Sample_cont_cov_dict, Sample_read_num, Sample_assembly_len


#######
#Plot step plot.
#######

def plot_cumul_step_curve(Sample_assembly_len, pwd):
  
    plt.figure(figsize=(4, 3))
    
    for Sample_type, len_data in Sample_assembly_len.items():
    
        plt.step(range(len(len_data)), len_data, label=Sample_type)
    
    plt.ylabel("Cumulative length, nt", size=12)
    plt.xlabel("Contig #", size=12)
    plt.legend(frameon=False)
    plt.tight_layout()
    plt.savefig(pwd+"Cumulative_length_of_assemblies.png", dpi=300)   
    
    return



#######
#Read string containing a pythonic list.
#######

def read_text_list(string):
    annot_ar=[]
    
    string_split=string.lstrip('[').rstrip(']').split('], [')
    for annot in string_split:
        annot=annot.lstrip('[').rstrip(']').split(', ')
        gene_name=annot[0]
        gene_e_value=float(annot[1])
        gene_score=float(annot[2])
        gene_bias=float(annot[3])
        
        annot_ar.append([gene_name, gene_e_value, gene_score, gene_bias])  
        
    return annot_ar


#######
#Read the pivot tabel with all data on defense systems genes.
#######

def read_DS_table_group(pwd):
    
    #Read sample-based table.
    Sample_cont_based=open(pwd+'Sample_ViralContig_Gene_Def_Gene_Annotation.txt', 'r')
    Sample_cont_based_dict={}
    for line in Sample_cont_based:
        line=line.rstrip().split('\t')
        Sample_type=line[0]
        Contig_name=line[1]
        Gene_name=line[2]
        Gene_annot=read_text_list(line[3])
        
        if Sample_type not in Sample_cont_based_dict:
            Sample_cont_based_dict[Sample_type]={Contig_name : {Gene_name : Gene_annot}}
        else:
            if Contig_name not in Sample_cont_based_dict[Sample_type]:
                Sample_cont_based_dict[Sample_type][Contig_name]={Gene_name : Gene_annot}
            else:
                Sample_cont_based_dict[Sample_type][Contig_name][Gene_name]=Gene_annot
                
    
    Sample_cont_based.close()    
    
    print(Sample_cont_based_dict)
    
    return Sample_cont_based_dict

Sample_cont_based_dict=read_DS_table_group(PWD)


#######
#Group genes in a contig by proximity. Get putative defense systems.
#######

def group_genes_by_proximity(contig_dict):
    
    gene_num_ar=[]
    gene_name_ar=[]
    
    for gene_name, gene_annot in contig_dict.items():
        gene_num=int(gene_name.split('|')[0].split('_')[1])
        gene_num_ar.append(gene_num)
        gene_name_ar.append(gene_name)
    
    gene_num_ar_sorted=sorted(gene_num_ar, reverse=False)
    
    gene_name_ar_sorted=[]
    for gene_num in gene_num_ar_sorted:
        gene_index=gene_num_ar.index(gene_num)
        gene_name=gene_name_ar[gene_index]
        gene_name_ar_sorted.append(gene_name)
        
    #print(gene_name_ar)
    #print(gene_name_ar_sorted)
    #print(gene_num_ar)
    #print(gene_num_ar_sorted)
    
    #Only one gene detected as defense in a contig.
    if len(gene_num_ar_sorted)==1:
        return 0
    
    else:
        def_genes_groups=[[gene_num_ar_sorted[0]]]
        def_genes_groups_names=[[gene_name_ar_sorted[0]]]
        
        for i in range(len(gene_num_ar_sorted)-1):
            if abs(gene_num_ar_sorted[i]-gene_num_ar_sorted[i+1])<5:
                def_genes_groups[-1].append(gene_num_ar_sorted[i+1])
                def_genes_groups_names[-1].append(gene_name_ar_sorted[i+1])
            else:
                def_genes_groups.append([gene_num_ar_sorted[i+1]])
                def_genes_groups_names.append([gene_name_ar_sorted[i+1]])
          
    return def_genes_groups_names


#######
#Detect defense systems in a proximity groups by common name of genes.
#######

def detect_DS_system(gene_prox_group, contig_dict, Cont_cov_depth, Sample_reads_num, Assembly_len):
    
    #Scaling constants.
    Assembly_len_const=10e6
    Number_of_reads_constant=10e6
    
    DS_sys_dict_redundant={}
    DS_sys_dict_nonredundant={}  
    DS_sys_dict_redundant_norm={}
    DS_sys_dict_nonredundant_norm={}      
    
    #Only one defense gene in a proximity group. Raw and normalized.
    if len(gene_prox_group)==1:
        gene_name=gene_prox_group[0]
        gene_annot=contig_dict[gene_name]
        gene_annot_sorted_by_evalue=sorted(gene_annot, key = lambda x: float(x[1]))
        
        #Non redundant counting.
        DS_gene_name=gene_annot_sorted_by_evalue[0][0]
        DS_sys_type=DS_gene_name.lstrip("'").rstrip("'").lstrip('"').rstrip('"').split('_')[0]
            
        DS_sys_dict_redundant[DS_sys_type]=1  
        DS_sys_dict_redundant_norm[DS_sys_type]=(1*Cont_cov_depth)*(Assembly_len_const/Assembly_len)*(Number_of_reads_constant/Sample_reads_num)

        return DS_sys_dict_redundant, DS_sys_dict_nonredundant, DS_sys_dict_redundant_norm, DS_sys_dict_nonredundant_norm
    
    else:
        DS_gene_name_nonredundant=[]
        
        for gene_name in gene_prox_group:
            gene_annot=contig_dict[gene_name]
            gene_annot_sorted_by_evalue=sorted(gene_annot, key = lambda x: float(x[1]))
            
            #Non redundant counting.
            DS_gene_name=gene_annot_sorted_by_evalue[0][0]
            if DS_gene_name not in DS_gene_name_nonredundant:
                DS_gene_name_nonredundant.append(DS_gene_name)
                
                DS_sys_type=DS_gene_name.lstrip("'").rstrip("'").lstrip('"').rstrip('"').split('_')[0]
                if DS_sys_type in DS_sys_dict_nonredundant:
                    DS_sys_dict_nonredundant[DS_sys_type]+=1
                else:
                    DS_sys_dict_nonredundant[DS_sys_type]=1   
                    DS_sys_dict_nonredundant_norm[DS_sys_type]=(1*Cont_cov_depth)*(Assembly_len_const/Assembly_len)*(Number_of_reads_constant/Sample_reads_num)
            
            #Redundant counting.
            DS_sys_type=gene_annot_sorted_by_evalue[0][0].lstrip("'").rstrip("'").lstrip('"').rstrip('"').split('_')[0]
            if DS_sys_type in DS_sys_dict_redundant:
                DS_sys_dict_redundant[DS_sys_type]+=1
                DS_sys_dict_redundant_norm[DS_sys_type]+=(1*Cont_cov_depth)*(Assembly_len_const/Assembly_len)*(Number_of_reads_constant/Sample_reads_num)
            else:
                DS_sys_dict_redundant[DS_sys_type]=1
                DS_sys_dict_redundant_norm[DS_sys_type]=(1*Cont_cov_depth)*(Assembly_len_const/Assembly_len)*(Number_of_reads_constant/Sample_reads_num)
                
        #print(DS_sys_dict_redundant)
        #print(DS_sys_dict_nonredundant) 
    
    return DS_sys_dict_redundant, DS_sys_dict_nonredundant, DS_sys_dict_redundant_norm, DS_sys_dict_nonredundant_norm


#######
#Count defense systems.
#Count defense genes.
#######

def count_DSs(sample_name, Count_genes, Count_systems, Count_genes_norm, Count_systems_norm, DS_sys_dict_redundant, DS_sys_dict_nonredundant, DS_sys_dict_redundant_norm, DS_sys_dict_nonredundant_norm):
    
    #Count defense genes. Raw and normalized.
    if sample_name not in Count_genes:
        Count_genes[sample_name]={}
        Count_genes_norm[sample_name]={}
        
    for sys_name, gene_num in DS_sys_dict_redundant.items():
        if sys_name in Count_genes[sample_name]:
            Count_genes[sample_name][sys_name]+=gene_num
            Count_genes_norm[sample_name][sys_name]+=DS_sys_dict_redundant_norm[sys_name]
        else:
            Count_genes[sample_name][sys_name]=gene_num
            Count_genes_norm[sample_name][sys_name]=DS_sys_dict_redundant_norm[sys_name]
    
    #Count defense systems.
    if sample_name not in Count_systems:
        Count_systems[sample_name]={}
        Count_systems_norm[sample_name]={}
    
    for sys_name, dif_gene_num in DS_sys_dict_nonredundant.items():
        if dif_gene_num>=2:
            print(f'System detected in {sample_name}: {sys_name}')
            if sys_name in Count_systems[sample_name]:
                Count_systems[sample_name][sys_name]+=1
                Count_systems_norm[sample_name][sys_name]+=DS_sys_dict_nonredundant_norm[sys_name]
            else:
                Count_systems[sample_name][sys_name]=1
                Count_systems_norm[sample_name][sys_name]=DS_sys_dict_nonredundant_norm[sys_name]
                
    return Count_genes, Count_systems, Count_genes_norm, Count_systems_norm


#######
#Count defense systems.
#Count defense genes.
#######

def DS_combine_in_dataframe(Count_genes, Count_systems, Count_genes_norm, Count_systems_norm):
    
    #Prepare initiating dictionary.
    #Samples_list=['BCT_MW_2016', 'BCT_MW_2018', 'BCT_H_panicea_2016', 'BCT_H_panicea_2018', 'BCT_H_sitiens_2016', 'BCT_H_sitiens_2018', 'BCT_I_palmata_2016', 'BCT_I_palmata_2018', 'VRS_MW_2018', 'VRS_H_panicea_2018', 'VRS_H_sitiens_2018', 'VRS_I_palmata_2018']
    #Samples_list=['MW_2016', 'MW_2018', 'H_panicea_2016', 'H_panicea_2018', 'H_sitiens_2016', 'H_sitiens_2018', 'I_palmata_2016', 'I_palmata_2018']
    Samples_list=['MW_2016', 'H_panicea_2016', 'H_sitiens_2016', 'I_palmata_2016', 'MW_2018', 'H_panicea_2018', 'H_sitiens_2018', 'I_palmata_2018']    
    #Samples_list=['MW_2016', 'MW_2018', 'H_panicea_2016', 'H_panicea_2018', 'H_sitiens_2016', 'H_sitiens_2018', 'I_palmata_2016']
    DS_list=['ABI', 'BREX', 'CRISPR-CAS', 'DISARM', 'DND', 'DRUANTIA', 'GABIJA', 'HACHIMAN', 'KIWA', 'LAMASSU', 'PAGOS', 'RM', 'SEPTU', 'SHEDU', 'TA', 'THOERIS', 'WADJET', 'ZORYA']
    
    Samples_DS_dict_genes={}
    
    for sample_type in Samples_list:
        Samples_DS_dict_genes[sample_type]={}
        for DS_type in DS_list:
            Samples_DS_dict_genes[sample_type][DS_type]=0
    
    Samples_DS_dict_systems=copy.deepcopy(Samples_DS_dict_genes)
    Samples_DS_dict_genes_norm=copy.deepcopy(Samples_DS_dict_genes)
    Samples_DS_dict_systems_norm=copy.deepcopy(Samples_DS_dict_genes)
    
    #Fill the dictionary with raw genes number (redundant) and normalized values.
    for sample_type, DS_dict in Count_genes.items():
        for DS_type, genes_num in DS_dict.items():
            Samples_DS_dict_genes[sample_type][DS_type]=genes_num
            Samples_DS_dict_genes_norm[sample_type][DS_type]=Count_genes_norm[sample_type][DS_type]
            
    #Fill the dictionary with systems number (nonredundant) and normalized values.
    for sample_type, DS_dict in Count_systems.items():
        for DS_type, sys_num in DS_dict.items():
            Samples_DS_dict_systems[sample_type][DS_type]=sys_num
            Samples_DS_dict_systems_norm[sample_type][DS_type]=Count_systems_norm[sample_type][DS_type]
    
    print(Samples_DS_dict_genes)
    print(Samples_DS_dict_systems)
    print(Samples_DS_dict_genes_norm)
    print(Samples_DS_dict_systems_norm)    
    
    Samples_DS_genes_df=pd.DataFrame(Samples_DS_dict_genes)
    Samples_DS_sys_df=pd.DataFrame(Samples_DS_dict_systems)
    Samples_DS_genes_df_norm=pd.DataFrame(Samples_DS_dict_genes_norm)
    Samples_DS_sys_df_norm=pd.DataFrame(Samples_DS_dict_systems_norm)    
    
    return Samples_DS_genes_df, Samples_DS_sys_df, Samples_DS_genes_df_norm, Samples_DS_sys_df_norm


#######
#Visualize heatmap.
#######

def heatmap_viz(DS_df, outpath, Title, DS_what):
    
    # definitions for the axes
    left, width = 0.1, 0.65
    bottom, height = 0.20, 0.60
    spacing = 0.01
    
    rect_scatter=[left, bottom, width, height]
    rect_histx=[left+0.12, bottom + height + spacing, width-0.245, 0.15]
    rect_histy=[left + width + spacing, bottom, 0.2, height]
    
    plt.figure(figsize=(8, 8))
    
    ax_scatter=plt.axes(rect_scatter)
    ax_scatter.tick_params(direction='out', top=False, right=False)
    ax_histx=plt.axes(rect_histx)
    ax_histx.tick_params(direction='out', length=2, labelbottom=False)
    ax_histy=plt.axes(rect_histy)
    ax_histy.tick_params(direction='out', length=2, labelleft=False)    
    
    Min=min([DS_df.min().min()])
    Max=max([DS_df.max().max()])
    
    print(Min, Max)
    
    print(DS_df)
    
    #sns.heatmap(Control_df, annot=True)
    #ax_scatter.set_title(Title)
    ax=sns.heatmap(DS_df, ax=ax_scatter, annot=False, square=True, vmin=Min, vmax=Max, cmap='viridis')
    
    Sample_type_count=DS_df.sum(axis=0).tolist()
    DS_type_count=DS_df.sum(axis=1).tolist()
    
    print(Sample_type_count, DS_type_count)
    
    ax_histx.bar(range(len(Sample_type_count)), Sample_type_count, color='#a7f6ff')
    ax_histx.set_xlim(-0.5, len(Sample_type_count)-0.4)
    ax_histx.set_ylabel(f'Number of {DS_what}', size=11)
    ax_histx.set_xticks(range(len(Sample_type_count)), minor=False)
    ax_histx.spines['top'].set_visible(False)
    ax_histx.spines['right'].set_visible(False)
    ax_histx.spines['bottom'].set_visible(True)
    ax_histx.spines['left'].set_visible(True)    
    ax_histy.barh(range(len(DS_type_count)), DS_type_count[::-1], color='#ffd9d0')
    ax_histy.set_ylim(-0.5, len(DS_type_count)-0.4)
    ax_histy.set_xlabel(f'Number of {DS_what}', size=11)
    ax_histy.set_yticks(range(len(DS_type_count)), minor=False)
    ax_histy.spines['top'].set_visible(False)
    ax_histy.spines['right'].set_visible(False)
    ax_histy.spines['bottom'].set_visible(True)
    ax_histy.spines['left'].set_visible(True)    
    
    #plt.tight_layout()
    plt.savefig(outpath, dpi=300)   
    
    return


#######
#Detects defence systems.
#Wrapper function.
#######

def find_DS_in_contigs(Sample_cont_based_dict, pwd):
    
    #Retrive data for contigs and assembly: coverage depth, length, number of reads.
    #Sample_cont_cov_dict, Sample_read_num, Sample_assembly_len=read_cov_data_file(pwd)
    
    #print(Sample_read_num)
    
    Count_genes={}
    Count_systems={}
    Count_genes_norm={}
    Count_systems_norm={}
    
    for sample_name, sample_dict in Sample_cont_based_dict.items():
        for contig_name, contig_dict in sample_dict.items():
            
            #Get sample and contig features for normalization.
            Cont_cov_depth=1
            Sample_reads_num=100000
            Assembly_len=10000
            
            print(f'{contig_name}')
            
            #Detect groups of defense genes by proximity.
            Genes_grouped_by_proximity=group_genes_by_proximity(contig_dict)
            if Genes_grouped_by_proximity!=0:
                for gene_prox_group in Genes_grouped_by_proximity:
                    
                    #Count number of defense genes and defense systems. Raw numbers and normalized numbers.
                    DS_sys_dict_redundant, DS_sys_dict_nonredundant, DS_sys_dict_redundant_norm, DS_sys_dict_nonredundant_norm=detect_DS_system(gene_prox_group, contig_dict, Cont_cov_depth, Sample_reads_num, Assembly_len)
                    Count_genes, Count_systems, Count_genes_norm, Count_systems_norm=count_DSs(sample_name, Count_genes, Count_systems, Count_genes_norm, Count_systems_norm, DS_sys_dict_redundant, DS_sys_dict_nonredundant, DS_sys_dict_redundant_norm, DS_sys_dict_nonredundant_norm)
    
    print(Count_genes)
    print(Count_systems)
    #print(Count_genes_norm)
    #print(Count_systems_norm)    
    
    ##Combine counted DS genes and DSs in dataframe.
    #Samples_DS_genes_df, Samples_DS_sys_df, Samples_DS_genes_df_norm, Samples_DS_sys_df_norm=DS_combine_in_dataframe(Count_genes, Count_systems, Count_genes_norm, Count_systems_norm)
    #
    ##Visualize raw abundance of DS genes and DS.
    #heatmap_viz(Samples_DS_genes_df, pwd+'VIR_Number_of_defense_systems_genes_in_datasets.png', 'Number of defense systems genes in datasets', '\ndefense genes')
    #heatmap_viz(Samples_DS_sys_df, pwd+'VIR_Number_of_defense_systems_in_datasets.png', 'Number of defense systems in datasets', '\ndefense systems')
    ##Visualize normalized abundance of DS genes and DS.
    #heatmap_viz(Samples_DS_genes_df_norm, pwd+'VIR_Normalized_number_of_defense_systems_genes_in_datasets.png', 'Normalized number of defense systems genes in datasets', '\ndefense genes, norm')
    #heatmap_viz(Samples_DS_sys_df_norm, pwd+'VIR_Normalized_number_of_defense_systems_in_datasets.png', 'Normalized number of defense systems in datasets', '\ndefense systems, norm')
    
    return

find_DS_in_contigs(Sample_cont_based_dict, PWD)
