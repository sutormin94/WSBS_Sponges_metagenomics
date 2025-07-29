###############################################
##Dmitry Sutormin, 2023##
##Prediction of viral contigs data##

## Part 1.
#1) Script takes as input data from CheckV, ViralComplete, and VirSorter2
#2) Filters data from CheckV.
#3) Overlaps filtered contigs from CheckV, ViralComplete, and VirSorter2.
#4) Returns mfa file with viral contigs of a high confidence (detected by at least 2 pipelines).

## Part 2.
#Script takes the output of DeepVirFinder and filters viral contigs with a high confidence.
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

###############
## Part 1. Initial filtering and overlapping (Virsorter2, CheckV, ViralVerify+ViralComplete)
###############

#######
#Data to be used.
#######

#Path to the working directory.
PWD="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\Spongy\Viral_contigs\Viral_contigs_prediction_and_filtering\Viruses_2016_predictions"
#PWD="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\MetaRUS\Coelentrazine_biosynthesis\Virus_prediction"

#Dataset name.
Dataset_name="I_palmata_2016_nanopore"

#Assembly method: metaviralspades, MetaSpades, MetFlye or any other.
Assembly_method='MetFlye'

#CheckV input data.
CheckV_inpath=f'{PWD}\\{Dataset_name}\\{Dataset_name}_quality_summary.tsv'

#VirSorter2 input data.
VirSorter2_inpath=f'{PWD}\\{Dataset_name}\\{Dataset_name}_final-viral-score.tsv'

#ViralComplete input data.
ViralComplete_inpath=f'{PWD}\\{Dataset_name}\\{Dataset_name}_viralverify_virus_result_table.csv'

#ViralVerify input contigs.
Input_contigs=f'{PWD}\\{Dataset_name}\\{Dataset_name}_input_with_circ.fasta'

#Output path
Outpath=f'{PWD}\\{Dataset_name}\\'



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

def read_checkV_and_filter(checkV_inpath, assembly_method):
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
    
    if assembly_method=='metaviralspades':
        
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
    
    elif assembly_method!='metaviralspades':
        Final_checkV_contigs_names=list(set(Viral_contigs))
        
    #print(Final_checkV_contigs_names)
    
    print(f'Number of Viral contigs detected by CheckV: {len(Final_checkV_contigs_names)}')
    
    return Final_checkV_contigs_names



#######
#Read VirSorter2 data.
#######

def read_virsorter2(virSorter2_inpath, assembly_method):
    filein=open(virSorter2_inpath, 'r')
    
    VirSorter2_contigs_names=[]
    
    for line in filein:
        line=line.split('\t')
        if line[0]!="seqname":
            contig_name=line[0].split('||')[0]
            VirSorter2_contigs_names.append(contig_name)
    
    filein.close()   
    
    #print(VirSorter2_contigs_names)
    
    if assembly_method=='metaviralspades':
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
    
    elif assembly_method!='metaviralspades':
        Final_virsorter2_contigs_names=list(set(VirSorter2_contigs_names))
        
    #print(Final_virsorter2_contigs_names)
    
    print(f'Number of Viral contigs detected by VirSorter2: {len(Final_virsorter2_contigs_names)}')
    
    return Final_virsorter2_contigs_names



#######
#Read ViralVerify -> ViralComplete data.
#######

def read_ViralComplete_data(viralComplete_inpath, assembly_method):
    filein=open(viralComplete_inpath, 'r')
    
    VirComp_contigs_names=[]
    
    for line in filein:
        line=line.split(',')
        contig_name=line[0]
        VirComp_contigs_names.append(contig_name)
    
    filein.close()   
    
    #print(VirComp_contigs_names)
    
    if assembly_method=='metaviralspades':
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
        
    elif assembly_method!='metaviralspades':
        Final_VirComp_contigs_names=list(set(VirComp_contigs_names))
        
    #print(Final_virsorter2_contigs_names)
    
    print(f'Number of Viral contigs detected by ViralVerify & ViralComplete: {len(Final_VirComp_contigs_names)}')
    
    return Final_VirComp_contigs_names
    


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


#######
#Initial filtering and overlapping wrapper function.
#######

def Init_filtering_wrapper(checkV_inpath, virSorter2_inpath, viralComplete_inpath, input_contigs, dataset_name, assembly_method, outpath):
    
    #Read CheckV data.
    Final_checkV_contigs_names=read_checkV_and_filter(checkV_inpath, assembly_method)
    
    #Read Virsorter2 data.
    VirSorter2_contigs_names=read_virsorter2(virSorter2_inpath, assembly_method)
    
    #Read ViralComplete data.
    ViralComplete_contigs_names=read_ViralComplete_data(viralComplete_inpath, assembly_method)
    
    #Combine and overlap.
    overlap_vir_conting_sets(Final_checkV_contigs_names, VirSorter2_contigs_names, ViralComplete_contigs_names, input_contigs, dataset_name, outpath)
    
    return

Init_filtering_wrapper(CheckV_inpath, VirSorter2_inpath, ViralComplete_inpath, Input_contigs, Dataset_name, Assembly_method, Outpath)



###############
## Part 2. Additional filtering using DeepVirFinder.
###############

#Path to the working directory.
PWD="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\Spongy\Viral_contigs\Viral_contigs_prediction_and_filtering\Viruses_2016_predictions"
#PWD="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\MetaRUS\Coelentrazine_biosynthesis\Virus_prediction"

#Dataset name.
Dataset_name="H_sitiens_2016_nanopore"

#Assembly method: metaviralspades, MetaSpades, MetFlye or any other.
Assembly_method='MetFlye'

#Score and P-value thresholds for data filtering.
Score_thr=0.6
Pvalue_thr=0.05

#DeepVirFinder input data.
DeepVirFinder_inpath=f'{PWD}\\{Dataset_name}\\{Dataset_name}_vir_contig_sets_overlap_gt4000bp_dvfpred.txt'

#Input contigs.
Input_contigs=f'{PWD}\\{Dataset_name}\\{Dataset_name}_vir_contig_sets_overlap.fasta'

#Output path
Outpath=f'{PWD}\{Dataset_name}\\'


#######
#Data filtering using DeepVirFinder.
#######

def DeepVirFinder_filter(deepvirfinder_inpath, dataset_name, score_thr, pvalue_thr, outpath):
    
    Filt_cont_list=[]
    
    i=0
    DVF_input=open(deepvirfinder_inpath, 'r')
    for line in DVF_input:
        line=line.rstrip().split('\t')
        if line[0] not in ['name']:
            Contig_name=line[0]
            Contig_len=int(line[1])
            Contig_score=float(line[2])
            Contig_pvalue=float(line[3])
            i+=1
            if (Contig_score>=score_thr) and (Contig_pvalue<pvalue_thr):
                Filt_cont_list.append(Contig_name)
            
    DVF_input.close()
    
    print(f'Number of {dataset_name} contigs survived with Score threshold {score_thr} and P-value threshold {pvalue_thr}: {len(Filt_cont_list)}/{i} ({np.round((float(len(Filt_cont_list))/i)*100, 1)}%)')
    
    return Filt_cont_list


#######
#DeepVirFinder data analysis wrapper.
#######

def DeepVirFinder_filter_set_wrapper(fasta_inpath, deepvirfinder_inpath, dataset_name, score_thr, pvalue_thr, outpath):
    
    #Read input contigs to be filtered.
    Input_contigs_dict=read_fasta_dict(fasta_inpath)
    
    #Filter contigs using DeepVirFinder predictions.
    Filt_cont_list=DeepVirFinder_filter(deepvirfinder_inpath, dataset_name, score_thr, pvalue_thr, outpath)
    
    #Return sequences of filtered viral contigs.
    Output_viral_contigs=open(f'{outpath}\{dataset_name}_sc_{score_thr}_pval_{pvalue_thr}_vir_contig_set_final.fasta', 'w')
    
    for contig_name in list(Filt_cont_list):
        Contig_seq=Input_contigs_dict[contig_name]
        Output_viral_contigs.write(f'>{contig_name}\n{Contig_seq}\n')
    
    Output_viral_contigs.close()       
    
    return

#DeepVirFinder_filter_set_wrapper(Input_contigs, DeepVirFinder_inpath, Dataset_name, Score_thr, Pvalue_thr, Outpath)