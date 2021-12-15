###############################################
##Dmitry Sutormin, 2021##
##Spacers blast analysis##

# Script reads spacers blast output and groups it by spacers source and target.
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

#########
##
## Part 1.
## Group spacers by the targeting genome; Group genomes by approaching spacers.
##
#########

#######
#Data to be used.
#######

#Blast output.
Blast_inpath="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\Spongy\Cctyper_Spacers_unclust_vs_Viral_contigs_clust\\all_viruses_vs_all_spacers_2018_evalue_filt"

#Spacers clustering data.
Spacers_clust_inpath="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\Spongy\CRISPR_spacers\\2021\Cctyper_Spacers\Spacers_clustering_cdhit\\all_spacers_cdhit_id95_c95.clstr"

#Output path.
Outpath="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\Spongy\Cctyper_Spacers_unclust_vs_Viral_contigs_clust\\"

#Blast output filering parameters: identity, number of allowed mismatches, evalue.
Ident_thr=85
Mm_thr=5
Evalue_thr=0.1

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
#Read spacers blast results.
#######

def read_spacers_blast_data(datapath, ident_thr, mm_thr, evalue_thr):
    filein=open(datapath, 'r')
    Spacer_Blast_dict={}
    Virus_Blast_dict={}
    for line in filein:
        line=line.rstrip('\r\n').split('\t')
        spacer_name=line[0]
        virus_name=line[1]
        ident=float(line[3])
        mm=int(line[5])
        evalue=float(line[11])
        if (ident>=ident_thr) and (mm<=mm_thr) and (evalue<=evalue_thr):
            if spacer_name not in Spacer_Blast_dict:
                Spacer_Blast_dict[spacer_name]=[virus_name]
            else:
                Spacer_Blast_dict[spacer_name].append(virus_name)
            
            if virus_name not in Virus_Blast_dict:
                Virus_Blast_dict[virus_name]=[spacer_name]
            else:
                Virus_Blast_dict[virus_name].append(spacer_name)
    
    filein.close()
    
    print(Spacer_Blast_dict)
    print('\n')
    print(Virus_Blast_dict)
    print('\n')
    
    return Spacer_Blast_dict, Virus_Blast_dict
    
    
#######
#Group blast data by spacers clustering.
#######

def clust_blast_data(Spacer_Blast_dict, Virus_Blast_dict, Clusters_dict, outpath, ident_thr, mm_thr, evalue_thr):
    
    #Group spacers-based dict.
    Spacer_Blast_dict_gr={}
    for cluster_name, spacers_ar in Clusters_dict.items():
        for spacer_name in spacers_ar:
            if spacer_name in Spacer_Blast_dict:
                if cluster_name not in Spacer_Blast_dict_gr:
                    Spacer_Blast_dict_gr[cluster_name]=[[spacer_name, Spacer_Blast_dict[spacer_name]]]
                else:
                    Spacer_Blast_dict_gr[cluster_name].append([spacer_name, Spacer_Blast_dict[spacer_name]])
    
    #Group viral-based dict.   
    Virus_Blast_dict_pre_gr={}
    for viral_name, spacers_name_ar in Virus_Blast_dict.items():
        for spacer_name in spacers_name_ar:
            for cluster_name, spacers_ar in Clusters_dict.items():
                if spacer_name in spacers_ar:
                    if viral_name not in Virus_Blast_dict_pre_gr:
                        Virus_Blast_dict_pre_gr[viral_name]=[cluster_name]
                    else:
                        Virus_Blast_dict_pre_gr[viral_name].append(cluster_name)
    
    Virus_Blast_dict_gr={}
    for viral_name, spacers_clust_ar in Virus_Blast_dict_pre_gr.items():
        spacers_gr_ar=[]
        for clust_name in list(set(spacers_clust_ar)):
            spacers_gr_ar.append(Clusters_dict[clust_name])
        Virus_Blast_dict_gr[viral_name]=spacers_gr_ar
    
    #for spacer_name, vir_contigs_at in Spacer_Blast_dict.items():
    
    print(Spacer_Blast_dict_gr)
    print('\n')
    print(Virus_Blast_dict_gr)
    
    spb_out=open(f'{outpath}\Cctyper_spacer_based_Vir_Ref_Seq_blast_results_grouping_ident{ident_thr}_mm{mm_thr}_eval{evalue_thr}.txt', 'w')
    for clust_name, spacer_vir_ar in Spacer_Blast_dict_gr.items():
        spb_out.write(f'{clust_name} : ')
        for spacer_vir_pair in spacer_vir_ar:
            spb_out.write(f'{spacer_vir_pair}\t')
        spb_out.write('\n')
    spb_out.close()
    
    virb_out=open(f'{outpath}\Vir_Ref_Seq_based_cctyper_spacer_blast_results_grouping_ident{ident_thr}_mm{mm_thr}_eval{evalue_thr}.txt', 'w')
    for vir_name, spacer_clust_ar in Virus_Blast_dict_gr.items():
        virb_out.write(f'{vir_name} : ')
        for spacer_ar in spacer_clust_ar:
            virb_out.write(f'{spacer_ar}\t')
        virb_out.write('\n')
    virb_out.close()    
    
    return
    
    
#Spacers_clusters_dict=read_cdhit_clusters_data(Spacers_clust_inpath)
#Spacer_Blast_dict, Virus_Blast_dict=read_spacers_blast_data(Blast_inpath, Ident_thr, Mm_thr, Evalue_thr)
#clust_blast_data(Spacer_Blast_dict, Virus_Blast_dict, Spacers_clusters_dict, Outpath, Ident_thr, Mm_thr, Evalue_thr)



#########
##
## Part 2.
## Attack-protect game.
##
#########

#######
#Data to be used.
#######

#Blast output.
Blast_inpath="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\Spongy\Cctyper_Spacers_unclust_vs_Viral_contigs_clust\\all_viruses_vs_all_spacers_2018_evalue_filt_no_FP_arrays"

#Output path.
Outpath="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\Spongy\Cctyper_Spacers_unclust_vs_Viral_contigs_clust\\"

#Blast output filering parameters: identity, number of allowed mismatches, evalue.
Ident_thr=95
Mm_thr=10
Evalue_thr=10


#######
#Read spacers blast results.
#######

def read_spacers_blast_data_attack_protect_game(datapath, ident_thr, mm_thr, evalue_thr):
    Samples_ar=['H_panicea_2018', 'H_sitiens_2018', 'I_palmata_2018', 'MW_2018']
    Vir_dict={}
    Sapcer_dict={}
    Inter_dict={}
    
    for i in Samples_ar:
        Vir_dict[i]={}
        Sapcer_dict[i]={}
        Inter_dict[i]={}        
        for j in Samples_ar:
            Vir_dict[i][j]=[]
            Sapcer_dict[i][j]=[]
            Inter_dict[i][j]=0           
    
    
    filein=open(datapath, 'r')

    for line in filein:
        line=line.rstrip('\r\n').split('\t')
        spacer_sample_name=line[0].lstrip('VRS_').lstrip('BCT_').split('_NODE')[0].split('_CUTOFF_')[0]
        virus_sample_name=line[1].lstrip('VRS_').lstrip('BCT_').split('_NODE')[0].split('_CUTOFF_')[0]
        virus_name=line[1]
        spacer_name=line[0]
        ident=float(line[3])
        mm=int(line[5])
        evalue=float(line[11])
        if (ident>=ident_thr) and (mm<=mm_thr) and (evalue<=evalue_thr):
            Inter_dict[spacer_sample_name][virus_sample_name]+=1
            if virus_name not in Vir_dict[spacer_sample_name][virus_sample_name]:
                Vir_dict[spacer_sample_name][virus_sample_name].append(virus_name)
            if spacer_name not in Sapcer_dict[spacer_sample_name][virus_sample_name]:
                Sapcer_dict[spacer_sample_name][virus_sample_name].append(spacer_name)

    
    filein.close()
    
    for sample_name_1, in_dict in Vir_dict.items():
        for sample_name_2, vir_ar in in_dict.items():
            Vir_dict[sample_name_1][sample_name_2]=len(vir_ar)
            
    for sample_name_1, in_dict in Sapcer_dict.items():
        for sample_name_2, spac_ar in in_dict.items():
            Sapcer_dict[sample_name_1][sample_name_2]=len(spac_ar)
    
    print(Inter_dict)
    print('\n')
    print(Vir_dict)
    print('\n')
    print(Sapcer_dict)
    print('\n')    
    
    return Inter_dict, Vir_dict, Sapcer_dict

Inter_dict, Vir_dict, Sapcer_dict=read_spacers_blast_data_attack_protect_game(Blast_inpath, Ident_thr, Mm_thr, Evalue_thr)