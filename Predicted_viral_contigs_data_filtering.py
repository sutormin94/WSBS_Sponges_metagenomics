###############################################
##Dmitry Sutormin, 2021##
##Prediction of viral contigs data##

#1) Script takes as input data from CheckV, ViralComplete, and VirSorter2
#2) Filters data from CheckV.
#3) Overlaps filtered contigs from CheckV, ViralComplete, and VirSorter2.
#4) Returns mfa file with viral contigs of a high confidence (detected by at least 2 pipelines).
###############################################


#######
#Packages to be imported.
#######

import numpy as np
import Bio
from Bio import SeqIO
import re
import matplotlib.pyplot as plt
from matplotlib_venn import venn2, venn3, venn3_circles
import csv
import pandas as pd
import seaborn as sns

#######
#Data to be used.
#######

#Dataset name.
Dataset_name="BCT_I_palmata_2018"

#CheckV input data.
CheckV_inpath=f'C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\Spongy\Viral_contigs\Viral_contigs_prediction_and_filtering\{Dataset_name}\CheckV_metaviralSPAdes_{Dataset_name}\quality_summary.tsv'

#VirSorter2 input data.
VirSorter2_inpath=f'C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\Spongy\Viral_contigs\Viral_contigs_prediction_and_filtering\{Dataset_name}\\virsorter2_metaviralSPades_{Dataset_name}\\final-viral-score.tsv'

#ViralComplete input data.
ViralComplete_inpath=f'C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\Spongy\Viral_contigs\Viral_contigs_prediction_and_filtering\{Dataset_name}\\viralComplete_after_viralVerify_metaviralSPAdes_{Dataset_name}\contigs_virus_result_table.csv'

#ViralVerify input contigs.
Input_contigs=f'C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\Spongy\Viral_contigs\Viral_contigs_prediction_and_filtering\{Dataset_name}\\viralVerify_metaviralSPAdes_{Dataset_name}\contigs_input_with_circ.fasta'

#Output path
Outpath=f'C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\Spongy\Viral_contigs\Viral_contigs_prediction_and_filtering\{Dataset_name}\\'



#######
#Read multy-fasta.
#######

def read_fasta_dict(fasta_path):
    fasta_dict={}
    fasta_input=open(fasta_path, 'r')
    for record in SeqIO.parse(fasta_input, "fasta"):
        fasta_id=str(record.name) 
        fasta_seq=str(record.seq) 
        fasta_dict[fasta_id]=fasta_seq
    fasta_input.close()
    
    print(f'Number of input contigs: {len(fasta_dict)}')
    
    return fasta_dict


#######
#Read CheckV data and filter.
#######

def read_checkV_and_filter(checkV_inpath):
    filein=open(checkV_inpath, 'r')
    
    Viral_contigs=[]
    
    for line in filein:
        line=line.split('\t')
        if line[0]!="contig_id":
            num_viral_genes=int(line[5])
            Checkv_quality=line[7]
            Warnings=line[13]
            if num_viral_genes>0:
                if Checkv_quality!="Not-determined":
                    Viral_contigs.append(line[0])
                    #print(line)
    
    filein.close()
    
    Contig_names_dict={}
    for contig_name in Viral_contigs:
        contig_name_base, cutoff=contig_name.split('_cutoff_')
        if contig_name_base in Contig_names_dict:
            Contig_names_dict[contig_name_base].append(cutoff)
        else:
            Contig_names_dict[contig_name_base]=[cutoff]
    
    #print(Contig_names_dict)
    
    Final_checkV_contigs_names=[]
    for contig_name_base, cutoff_ar in Contig_names_dict.items():
        cutoff_value_ar=[]
        ind=0
        for cutoff in cutoff_ar:
            cutoff_value, lin_type=cutoff.split('_type_')
            if cutoff_value=='metaplasmid':
                Final_checkV_contigs_names.append(contig_name_base+"_cutoff_metaplasmid_type_circular")
                ind=1
            else:
                cutoff_value_ar.append(int(cutoff_value))
        if ind==0:
            max_cutoff_value=max(cutoff_value_ar)
            Final_checkV_contigs_names.append(f'{contig_name_base}_cutoff_{max_cutoff_value}_type_{lin_type}')
        
    #print(Final_checkV_contigs_names)
    
    print(f'Number of Viral contigs detected by CheckV: {len(Final_checkV_contigs_names)}')
    
    return Final_checkV_contigs_names

Final_checkV_contigs_names=read_checkV_and_filter(CheckV_inpath)



#######
#Read VirSorter2 data.
#######

def read_virsorter2(virSorter2_inpath):
    filein=open(virSorter2_inpath, 'r')
    
    VirSorter2_contigs_names=[]
    
    for line in filein:
        line=line.split('\t')
        if line[0]!="seqname":
            contig_name=line[0].split('||')[0]
            VirSorter2_contigs_names.append(contig_name)
    
    filein.close()   
    
    #print(VirSorter2_contigs_names)
    
    Contig_names_dict={}
    for contig_name in VirSorter2_contigs_names:
        contig_name_base, cutoff=contig_name.split('_cutoff_')
        if contig_name_base in Contig_names_dict:
            Contig_names_dict[contig_name_base].append(cutoff)
        else:
            Contig_names_dict[contig_name_base]=[cutoff]    
            
    #print(Contig_names_dict)
            
    Final_virsorter2_contigs_names=[]
    for contig_name_base, cutoff_ar in Contig_names_dict.items():
        cutoff_value_ar=[]
        ind=0
        for cutoff in cutoff_ar:
            cutoff_value, lin_type=cutoff.split('_type_')
            if cutoff_value=='metaplasmid':
                Final_virsorter2_contigs_names.append(contig_name_base+"_cutoff_metaplasmid_type_circular")
                ind=1
            else:
                cutoff_value_ar.append(int(cutoff_value))
        if ind==0:
            max_cutoff_value=max(cutoff_value_ar)
            Final_virsorter2_contigs_names.append(f'{contig_name_base}_cutoff_{max_cutoff_value}_type_{lin_type}')
        
    #print(Final_virsorter2_contigs_names)
    
    print(f'Number of Viral contigs detected by VirSorter2: {len(Final_virsorter2_contigs_names)}')
    
    return Final_virsorter2_contigs_names

VirSorter2_contigs_names=read_virsorter2(VirSorter2_inpath)



#######
#Read ViralVerify -> ViralComplete data.
#######

def read_ViralComplete_data(viralComplete_inpath):
    filein=open(viralComplete_inpath, 'r')
    
    VirComp_contigs_names=[]
    
    for line in filein:
        line=line.split(',')
        contig_name=line[0]
        VirComp_contigs_names.append(contig_name)
    
    filein.close()   
    
    #print(VirComp_contigs_names)
    
    Contig_names_dict={}
    for contig_name in VirComp_contigs_names:
        contig_name_base, cutoff=contig_name.split('_cutoff_')
        if contig_name_base in Contig_names_dict:
            Contig_names_dict[contig_name_base].append(cutoff)
        else:
            Contig_names_dict[contig_name_base]=[cutoff]    
            
    #print(Contig_names_dict)
            
    Final_VirComp_contigs_names=[]
    for contig_name_base, cutoff_ar in Contig_names_dict.items():
        cutoff_value_ar=[]
        ind=0
        for cutoff in cutoff_ar:
            cutoff_value, lin_type=cutoff.split('_type_')
            if cutoff_value=='metaplasmid':
                Final_VirComp_contigs_names.append(contig_name_base+"_cutoff_metaplasmid_type_circular")
                ind=1
            else:
                cutoff_value_ar.append(int(cutoff_value))
        if ind==0:
            max_cutoff_value=max(cutoff_value_ar)
            Final_VirComp_contigs_names.append(f'{contig_name_base}_cutoff_{max_cutoff_value}_type_{lin_type}')
        
    #print(Final_VirComp_contigs_names)
    
    print(f'Number of Viral contigs detected by ViralVerify & ViralComplete: {len(Final_VirComp_contigs_names)}')
    
    return Final_VirComp_contigs_names
    

ViralComplete_contigs_names=read_ViralComplete_data(ViralComplete_inpath)


#######
#Overlap sets of viral contigs detected by different pipelines.
#######

def overlap_vir_conting_sets(Final_checkV_contigs_names, Final_virsorter2_contigs_names, Final_VirComp_contigs_names, fasta_inpath, dataset_name, outpath):
    
    venn3([set(Final_checkV_contigs_names), set(Final_virsorter2_contigs_names), set(Final_VirComp_contigs_names)], set_labels=('CheckV', 'VirSorter2', 'ViralVerify&\nViralComplete'))
    venn3_circles([set(Final_checkV_contigs_names), set(Final_virsorter2_contigs_names), set(Final_VirComp_contigs_names)])
    plt.title(f'{Dataset_name}')
    plt.show()
    
    plt.savefig(f'{outpath}\{dataset_name}_vir_contig_sets_overlap.png', dpi=400, figsize=(4, 4))
    plt.savefig(f'{outpath}\{dataset_name}_vir_contig_sets_overlap.svg', dpi=400, figsize=(4, 4))
    
    #Get contigs in the intersection.
    CheckV_and_VirSorter2=set(Final_checkV_contigs_names).intersection(set(Final_virsorter2_contigs_names))
    VirComp_and_Virsorter2=set(Final_VirComp_contigs_names).intersection(set(Final_virsorter2_contigs_names))
    CheckV_and_VirComp=set(Final_checkV_contigs_names).intersection(set(Final_VirComp_contigs_names))
    
    Double_overlap=set(list(CheckV_and_VirSorter2)+list(VirComp_and_Virsorter2)+list(CheckV_and_VirComp))
    
    #Read input contig sequences.
    Input_contigs_dict=read_fasta_dict(fasta_inpath)
    
    #Return sequences of filtered viral contigs.
    Output_viral_contigs=open(f'{outpath}\{dataset_name}_vir_contig_sets_overlap.fasta', 'w')
    
    for contig_name in list(Double_overlap):
        Contig_seq=Input_contigs_dict[contig_name]
        Output_viral_contigs.write(f'>{contig_name}\n{Contig_seq}\n')
    
    Output_viral_contigs.close()   
    
    return

overlap_vir_conting_sets(Final_checkV_contigs_names, VirSorter2_contigs_names, ViralComplete_contigs_names, Input_contigs, Dataset_name, Outpath)