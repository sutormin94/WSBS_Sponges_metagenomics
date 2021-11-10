###############################################
##Dmitry Sutormin, 2021##
##Anti phage defense systems analysis##

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
## Prepare PADS Arsenal database for hmm profiles construction. Group genes by name.
##
###############
###############

#######
#Data to be used.
#######

#PWD.
PWD="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\Spongy\Anti_phage_systems\\"


#######
#Read fasta.
#######

def read_fasta_dict(fasta_path):
    fasta_dict={}
    fasta_in=open(fasta_path, 'r')
    for record in SeqIO.parse(fasta_in, "fasta"):
        prot_id=str(record.name) 
        prot_group=prot_id.split(',')[2]
        prot_seq=str(record.seq) 
        if prot_group in fasta_dict:
            fasta_dict[prot_group][prot_id]=prot_seq
        else:
            fasta_dict[prot_group]={prot_id : prot_seq}
    fasta_in.close()
    
    print(len(fasta_dict))
    
    return fasta_dict



#######
#Read PADS arsenal files, group proteins, output files with orthologs.
#######

def read_pads_arsenal_group(inpath):
    if os.path.isdir(f'{inpath}PADS_Arsenal_grouped')==False:
        os.mkdir(f'{inpath}PADS_Arsenal_grouped')
    for pd_system_name in os.listdir(f'{inpath}PADS_Arsenal\\'):
        infile_path=f'{inpath}PADS_Arsenal\\{pd_system_name}\\PADS_Arsenal_Bacteria_{pd_system_name}_v1_2019.09.09.faa'
        protein_dict=read_fasta_dict(infile_path)
        
        if os.path.isdir(f'{inpath}PADS_Arsenal_grouped\\{pd_system_name}')==False:
            os.mkdir(f'{inpath}PADS_Arsenal_grouped\\{pd_system_name}')
        
        for gene_group, prot_dict in protein_dict.items():
            gene_group=gene_group.replace('/', '_or_')
            fileout=open(f'{inpath}PADS_Arsenal_grouped\\{pd_system_name}\PADS_Arsenal_Bacteria_{pd_system_name}_{gene_group}.fasta', 'w')
            for prot_id, prot_seq in prot_dict.items():
                fileout.write(f'>{prot_id}\n{prot_seq}\n')
                
            fileout.close()
        
    return

#read_pads_arsenal_group(PWD)



###############
###############
##
## Part 2.
## Calculate the amount of defense genes relative to all annotated genes.
## Group annotation for one gene. Group defense genes by contigs, group defense genes by sample.
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
#Read hmmsearch output and collect hits in a dictionary.
#######

def read_hmm_file(filein_path, PADS_annotation, e_value_thr):
    
    filein=open(filein_path, 'r')
    for line in filein:
        if line[0]!='#':
            line=line.rstrip().split(' ')
            
            line_filt=[]
            for ele in line:
                if ele!='':
                    line_filt.append(ele)
                    
            full_orf_id=line_filt[0]
            query_name=line_filt[2]
            e_value=float(line_filt[4])
            score=float(line_filt[5])
            bias=float(line_filt[6])
            if e_value<e_value_thr:
                if full_orf_id in PADS_annotation:
                    PADS_annotation[full_orf_id].append([query_name, e_value, score, bias])
                else:
                    PADS_annotation[full_orf_id]=[[query_name, e_value, score, bias]]
    
    filein.close()
    return 


#######
#Group by sample.
#######

def group_annot_by_sample(pwd, PADS_annotation):
    
    prot_id_dict={}
    Sample_type_dict={}
    
    for input_file_name in os.listdir(f'{pwd}\\mgm_prot_all\\Full_annotation_files'):
        filein_path=f'{pwd}\\mgm_prot_all\\Full_annotation_files\\{input_file_name}'  
        print(filein_path)    
        
        fasta_id_ar=[]
        fasta_in=open(filein_path, 'r')
        for line in fasta_in:
            if line[0]=='>':
                prot_id=line.rstrip().lstrip('>')
                #print(prot_id)
                
                if prot_id in fasta_id_ar:
                    print('Intra-sample name duplication!')
                else:
                    fasta_id_ar.append(prot_id)     
        
        prot_id_dict[input_file_name.rstrip('_prot.faa')]=fasta_id_ar
        Sample_type_dict[input_file_name.rstrip('_prot.faa')]={}

        fasta_in.close()
               
    #Count defense genes by sample type.
    for gene_id, gene_annot in PADS_annotation.items():
        for sample_id, sample_gene_list in prot_id_dict.items():
            if gene_id in sample_gene_list:
                Sample_type_dict[sample_id][gene_id]=gene_annot
    
    fileout=open(pwd+"Sample_Gene_Def_Gene_Annotation.txt", 'w')           
    for sample_id, sample_gene_list in Sample_type_dict.items():
        print(f'{sample_id} : {len(sample_gene_list)}/{len(prot_id_dict[sample_id])} defense genes found')
        for gene_id, gene_annot in sample_gene_list.items():
            fileout.write(f'{sample_id}\t{gene_id}\t{gene_annot}\n')
        
    fileout.close()
        
    return Sample_type_dict


#######
#Group by contig.
#######

def group_annot_by_contig(pwd, PADS_annotation):
    
    Contig_dict={}
    
    #Group defense genes by contig id.
    for gene_id, gene_annot in PADS_annotation.items():  
        
        contig_id=gene_id.split('_>')[1]
        
        if contig_id in Contig_dict:
            Contig_dict[contig_id][gene_id]=gene_annot
        else:
            Contig_dict[contig_id]={}
            Contig_dict[contig_id][gene_id]=gene_annot
     
    print(f'Total number of contigs with defense genes detected : {len(Contig_dict)}')   
    
    #Count defense genes in contigs.
    fileout=open(pwd+"Contig_Gene_Def_Gene_Annotation.txt", 'w')   
    Cont_def_g_num_ar=[]
    for contig_id, contig_annot in Contig_dict.items():
        Cont_def_g_num_ar.append(len(contig_annot))
        for gene_id, gene_annot in contig_annot.items():
            fileout.write(f'{contig_id}\t{gene_id}\t{gene_annot}\n')
        
    fileout.close()
        
    plt.hist(Cont_def_g_num_ar)
    plt.yscale('log')
    plt.xlabel('Number of defense genes in a contig')
    plt.ylabel('Number of contigs')
    plt.show()
    plt.savefig(f'{pwd}Defense_genes_in_contigs_distribution.png') 
        
    
    return Contig_dict


#######
#Read all hmmsearch output files and assemle the data.
#######

def take_all_data(pwd, e_value_thr):
    
    #Read defense systems annotation.
    PADS_annotation={}
    
    for input_file_name in os.listdir(f'{pwd}\\PADs_search_in_mgm_prot_filtered'):
        filein_path=f'{pwd}\\PADs_search_in_mgm_prot_filtered\\{input_file_name}'  
        print(filein_path)
        
        read_hmm_file(filein_path, PADS_annotation, e_value_thr)
        
    print(PADS_annotation)
    print(f'Total number of defense genes detected : {len(PADS_annotation)}')
    
    #Read the full annotation of all input metagenomes and count defense genes by sample.
    Sample_type_dict=group_annot_by_sample(pwd, PADS_annotation)
    
    #Group defense genes by contigs.
    Contig_dict=group_annot_by_contig(pwd, PADS_annotation)
     
    return

#take_all_data(PWD, E_value_thr)


###############
###############
##
## Part 3.
## Assemble {Sample}\t{contig}\t{gene}\t{gene annotation} table.
##
###############
###############

#PWD.
PWD="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\Spongy\Anti_phage_systems\\"


def combine_data_in_one_table(pwd):
    
    #Read sample-based table.
    Sample_based=open(pwd+'Sample_Gene_Def_Gene_Annotation.txt', 'r')
    Sample_based_dict={}
    for line in Sample_based:
        line=line.rstrip().split('\t')
        Sample_type=line[0]
        Gene_name=line[1]
        Gene_annot=line[2]
        
        if Gene_name not in Sample_based_dict:
            Sample_based_dict[Gene_name]=[Sample_type, Gene_annot]
        else:
            print(f'Ambigous gene name! {Gene_name}')
    
    Sample_based.close()
    
    #Read contig-based table.
    Contig_based=open(pwd+'Contig_Gene_Def_Gene_Annotation.txt', 'r')
    Contig_based_dict={}
    for line in Contig_based:
        line=line.rstrip().split('\t')  
        Contig_name=line[0]
        Gene_name=line[1]
        Gene_annot=line[2]  
        
        if Gene_name not in Contig_based_dict:
            Contig_based_dict[Gene_name]=[Contig_name, Gene_annot]
        else:
            print(f'Ambigous gene name! {Gene_name}')        
    
    Contig_based.close()
    
    #Combine and write resultant table.
    Sample_contig_based=open(pwd+'Sample_Contig_Gene_Def_Gene_Annotation.txt', 'w')
    for Gene_name, cont_gene_data in Contig_based_dict.items():
        sample_gene_data=Sample_based_dict[Gene_name]
        Sample_contig_based.write(f'{sample_gene_data[0]}\t{cont_gene_data[0]}\t{Gene_name}\t{cont_gene_data[1]}\n')
    
    Sample_contig_based.close()
    
    return
    
#combine_data_in_one_table(PWD)


###############
###############
##
## Part 4.
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
    
    pwd_in=pwd+"Defense_systems\Contigs_coverage_5000_2018\\"
    
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
    plt.legend(frameon=False, markerscale=10, handlelength=1.7, handletextpad=1, fontsize=8)
    plt.tight_layout()
    plt.savefig(pwd+"Cumulative_length_of_assemblies_2018.png", dpi=300) 
    plt.savefig(pwd+"Cumulative_length_of_assemblies_2018.svg", dpi=300)   
    
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
    Sample_cont_based=open(pwd+'PADS_based_DS_detection\Sample_Contig_Gene_Def_Gene_Annotation_filtered_2018.txt', 'r')
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
    
    #print(Sample_cont_based_dict)
    
    return Sample_cont_based_dict

Sample_cont_based_dict=read_DS_table_group(PWD)


#######
#Group genes in a contig by proximity. Get putative defense systems.
#######

def group_genes_by_proximity(contig_dict):
    
    Proximity_distance_thr=5
    
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
            if abs(gene_num_ar_sorted[i]-gene_num_ar_sorted[i+1])<Proximity_distance_thr:
                def_genes_groups[-1].append(gene_num_ar_sorted[i+1])
                def_genes_groups_names[-1].append(gene_name_ar_sorted[i+1])
            else:
                def_genes_groups.append([gene_num_ar_sorted[i+1]])
                def_genes_groups_names.append([gene_name_ar_sorted[i+1]])
          
    return def_genes_groups_names


#######
#Detect defense systems in a proximity groups by common name of genes.
#######

def detect_DS_system(gene_prox_group, contig_dict, Cont_cov_depth, Sample_reads_num, Assembly_len, Putative_systems_output, Putative_systems_output_long, sample_name, contig_name):
    
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
        
        #Output defense system.
        Putative_systems_output.write(f'{sample_name}\t{contig_name}\t{gene_name}\t{DS_gene_name}\n')

        return DS_sys_dict_redundant, DS_sys_dict_nonredundant, DS_sys_dict_redundant_norm, DS_sys_dict_nonredundant_norm
    
    else:
        DS_gene_name_nonredundant=[]
        
        #Output defense system.
        Putative_systems_output.write(f'{sample_name}\t{contig_name}')
        Putative_systems_output_long.write(f'{sample_name}\t{contig_name}')
        
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
            Putative_systems_output.write(f'\t{gene_name}\t{DS_gene_name}')
            Putative_systems_output_long.write(f'\t{gene_name}\t{DS_gene_name}')
            
            DS_sys_type=gene_annot_sorted_by_evalue[0][0].lstrip("'").rstrip("'").lstrip('"').rstrip('"').split('_')[0]
            if DS_sys_type in DS_sys_dict_redundant:
                DS_sys_dict_redundant[DS_sys_type]+=1
                DS_sys_dict_redundant_norm[DS_sys_type]+=(1*Cont_cov_depth)*(Assembly_len_const/Assembly_len)*(Number_of_reads_constant/Sample_reads_num)
            else:
                DS_sys_dict_redundant[DS_sys_type]=1
                DS_sys_dict_redundant_norm[DS_sys_type]=(1*Cont_cov_depth)*(Assembly_len_const/Assembly_len)*(Number_of_reads_constant/Sample_reads_num)
                
        #print(DS_sys_dict_redundant)
        #print(DS_sys_dict_nonredundant) 
        Putative_systems_output.write('\n')
        Putative_systems_output_long.write('\n')         
    
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
    #Samples_list=['MW_2016', 'H_panicea_2016', 'H_sitiens_2016', 'I_palmata_2016', 'MW_2018', 'H_panicea_2018', 'H_sitiens_2018', 'I_palmata_2018']    
    #Samples_list=['MW_2016', 'MW_2018', 'H_panicea_2016', 'H_panicea_2018', 'H_sitiens_2016', 'H_sitiens_2018', 'I_palmata_2016']
    Samples_list=['MW_2018', 'H_panicea_2018', 'H_sitiens_2018', 'I_palmata_2018']
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

def heatmap_viz(DS_df, outpath, Title, DS_what, Adj_params):
    
    # definitions for the axes
    left, width = 0.1, 0.65
    bottom, height = 0.20, 0.60
    spacing = 0.01
    
    rect_scatter=[left, bottom, width, height]
    rect_histx=[left+Adj_params[0], bottom + height + spacing, width+Adj_params[1], 0.15]
    rect_histy=[left + width + spacing, bottom, 0.2, height]
    
    plt.figure(figsize=(7, 8))
    
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
    
    Sample_type_count=DS_df.sum(axis=0).tolist()
    DS_type_count=DS_df.sum(axis=1).tolist()
    
    print(Sample_type_count, DS_type_count) 
    
    DS_types_list=DS_df.index.tolist()
    print(DS_types_list)
    
    #sns.heatmap(Control_df, annot=True)
    #ax_scatter.set_title(Title)
    ax=sns.heatmap(DS_df, ax=ax_scatter, annot=False, square=True, vmin=Min, vmax=Max, cmap='viridis')
    ax.set_yticks(np.array(range(len(DS_types_list)))+0.5, minor=False)
    ax.set_yticklabels(DS_types_list)
    
    
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
    Sample_cont_cov_dict, Sample_read_num, Sample_assembly_len=read_cov_data_file(pwd)
    
    #Output putative defense systems.
    Putative_systems_output=open(f'{pwd}\PADS_based_DS_detection\Sample_contig_system_genes_annotations.txt', 'w')
    Putative_systems_output_long=open(f'{pwd}\PADS_based_DS_detection\Sample_contig_system_genes_annotations_long.txt', 'w')
    
    print(Sample_read_num)
    
    Count_genes={}
    Count_systems={}
    Count_genes_norm={}
    Count_systems_norm={}
    
    
    for sample_name, sample_dict in Sample_cont_based_dict.items():
        for contig_name, contig_dict in sample_dict.items():
            
            #Get sample and contig features for normalization.
            Cont_cov_depth=Sample_cont_cov_dict[sample_name][contig_name][2]
            Sample_reads_num=Sample_read_num[sample_name]
            Assembly_len=Sample_assembly_len[sample_name][-1]
            
            #Detect groups of defense genes by proximity.
            Genes_grouped_by_proximity=group_genes_by_proximity(contig_dict)
            
            if Genes_grouped_by_proximity!=0:
                for gene_prox_group in Genes_grouped_by_proximity:
                    
                    #Count number of defense genes and defense systems. Raw numbers and normalized numbers.
                    DS_sys_dict_redundant, DS_sys_dict_nonredundant, DS_sys_dict_redundant_norm, DS_sys_dict_nonredundant_norm=detect_DS_system(gene_prox_group, contig_dict, Cont_cov_depth, Sample_reads_num, Assembly_len, Putative_systems_output, Putative_systems_output_long, sample_name, contig_name)
                    Count_genes, Count_systems, Count_genes_norm, Count_systems_norm=count_DSs(sample_name, Count_genes, Count_systems, Count_genes_norm, Count_systems_norm, DS_sys_dict_redundant, DS_sys_dict_nonredundant, DS_sys_dict_redundant_norm, DS_sys_dict_nonredundant_norm)
    
    Putative_systems_output.close()
    Putative_systems_output_long.close()
    
    print(Count_genes)
    print(Count_systems)
    print(Count_genes_norm)
    print(Count_systems_norm)    
    
    #Combine counted DS genes and DSs in dataframe.
    Samples_DS_genes_df, Samples_DS_sys_df, Samples_DS_genes_df_norm, Samples_DS_sys_df_norm=DS_combine_in_dataframe(Count_genes, Count_systems, Count_genes_norm, Count_systems_norm)
    
    #Visualize raw abundance of DS genes and DS.
    #Polishing parameters.
    Adj_params=[0.365, -0.49]
    heatmap_viz(Samples_DS_genes_df, pwd+'PADS_based_DS_detection\\2018_Number_of_defense_systems_genes_in_datasets.svg', 'Number of defense systems genes in datasets', '\ndefense genes', Adj_params)
    heatmap_viz(Samples_DS_sys_df, pwd+'PADS_based_DS_detection\\2018_Number_of_defense_systems_in_datasets.svg', 'Number of defense systems in datasets', '\ndefense systems', Adj_params)
    #Visualize normalized abundance of DS genes and DS.
    heatmap_viz(Samples_DS_genes_df_norm, pwd+'PADS_based_DS_detection\\2018_Normalized_number_of_defense_systems_genes_in_datasets.svg', 'Normalized number of defense systems genes in datasets', '\ndefense genes, norm', Adj_params)
    heatmap_viz(Samples_DS_sys_df_norm, pwd+'PADS_based_DS_detection\\2018_Normalized_number_of_defense_systems_in_datasets.svg', 'Normalized number of defense systems in datasets', '\ndefense systems, norm', Adj_params)
    
    return

#find_DS_in_contigs(Sample_cont_based_dict, PWD)




###############
###############
##
## Part 5.
## Read and analyse DefenseFinder output.
##
###############
###############

#######
#Data to be used.
#######

#PWD.
PWD="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\Spongy\Anti_phage_systems\DefenseFinder\\"


#######
#Read DefenseFinder output, group by sample, by defense system type.
#######

def Read_DefenseFinder_output_group_norm_data(pwd):
    
    #Retrive data for contigs and assembly: coverage depth, length, number of reads.
    Sample_cont_cov_dict, Sample_read_num, Sample_assembly_len=read_cov_data_file('C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\Spongy\Anti_phage_systems\\')    
    
    #Prepare dataframe.
    Samples_list=['MW_2018', 'H_panicea_2018', 'H_sitiens_2018', 'I_palmata_2018']
    DS_list=['Abi', 'BREX', 'Cas', 'DISARM', 'Dnd', 'Druantia', 'Gabija', 'Hachiman', 'Kiwa', 'Lamassu', 'PAGOS', 'RM', 'Septu', 'Shedu', 'TA', 'Thoeris', 'Wadjet', 'Zorya', 'AVAST', 'CBASS', 'dGTPase', 'Nhi', 'Gao', 'Lit', 'Retron', 'Rst', 'DRT', 'Dsr', 'Viperin', 'PrrC', 'PARIS', 'BtsA', 'dCTPdeaminase']
    
    Samples_DS_dict_genes={}
    
    for sample_type in Samples_list:
        Samples_DS_dict_genes[sample_type]={}
        for DS_type in DS_list:
            Samples_DS_dict_genes[sample_type][DS_type]=0
    
    Samples_DS_dict_systems=copy.deepcopy(Samples_DS_dict_genes)
    Samples_DS_dict_genes_norm=copy.deepcopy(Samples_DS_dict_genes)
    Samples_DS_dict_systems_norm=copy.deepcopy(Samples_DS_dict_genes) 
    
    #Scaling constants.
    Assembly_len_const=10e6
    Number_of_reads_constant=10e6    
    
    #Retrive DefenseFinder output.
    for DF_dir in os.listdir(f'{pwd}Bacterial_fraction\\'):
        sample_name=DF_dir.split('__')[0]
        for Job_ID in os.listdir(f'{pwd}Bacterial_fraction\\{DF_dir}'):
            DF_output=open(f'{pwd}Bacterial_fraction\\{DF_dir}\\{Job_ID}\defense_finder_systems.tsv', 'r')
            for line in DF_output:
                line=line.rstrip().split('\t')
                if line[0]!='sys_id':
                    System_type=line[1]
                    if System_type[:3]=='Abi':
                        System_type='Abi'
                    elif System_type[:3]=='Gao':
                        System_type='Gao'
                    elif System_type[:3]=='Rst':
                        System_type='Rst'                    
                    Contig_name=line[3].split('__>')[1]
                    Num_genes=int(line[6])
                    
                    #Raw counting of DS and DS genes.
                    Samples_DS_dict_genes[sample_name][System_type]+=Num_genes
                    Samples_DS_dict_systems[sample_name][System_type]+=1
                    
                    #Normalized counting of DS and DS genes.
                    Cont_cov_depth=Sample_cont_cov_dict[sample_name][Contig_name][2]
                    Assembly_len=Sample_assembly_len[sample_name][-1]
                    Sample_reads_num=Sample_read_num[sample_name]
                    
                    Samples_DS_dict_genes_norm[sample_name][System_type]+=(Num_genes*Cont_cov_depth)*(Assembly_len_const/Assembly_len)*(Number_of_reads_constant/Sample_reads_num)
                    Samples_DS_dict_systems_norm[sample_name][System_type]+=(1*Cont_cov_depth)*(Assembly_len_const/Assembly_len)*(Number_of_reads_constant/Sample_reads_num)
        
            DF_output.close()
            
    #Prepare dataframes for plotting.
    Samples_DS_genes_df=pd.DataFrame(Samples_DS_dict_genes)
    Samples_DS_sys_df=pd.DataFrame(Samples_DS_dict_systems)
    Samples_DS_genes_df_norm=pd.DataFrame(Samples_DS_dict_genes_norm)
    Samples_DS_sys_df_norm=pd.DataFrame(Samples_DS_dict_systems_norm)  
    
    #Visualize raw abundance of DS genes and DS.
    #Polishing parameters.
    Adj_params=[0.437, -0.565]    
    heatmap_viz(Samples_DS_genes_df, pwd+'\\2018_DefenseFinder_Number_of_defense_systems_genes_in_datasets.svg', 'Number of defense systems genes in datasets', '\ndefense genes', Adj_params)
    heatmap_viz(Samples_DS_sys_df, pwd+'\\2018_DefenseFinder_Number_of_defense_systems_in_datasets.svg', 'Number of defense systems in datasets', '\ndefense systems', Adj_params)
    #Visualize normalized abundance of DS genes and DS.
    heatmap_viz(Samples_DS_genes_df_norm, pwd+'\\2018_DefenseFinder_Normalized_number_of_defense_systems_genes_in_datasets.svg', 'Normalized number of defense systems genes in datasets', '\ndefense genes, norm', Adj_params)
    heatmap_viz(Samples_DS_sys_df_norm, pwd+'\\2018_DefenseFinder_Normalized_number_of_defense_systems_in_datasets.svg', 'Normalized number of defense systems in datasets', '\ndefense systems, norm', Adj_params)
          
    return

Read_DefenseFinder_output_group_norm_data(PWD)



###############
###############
##
## Part 6.
## Read and analyse PADLOC output.
##
###############
###############

#######
#Data to be used.
#######

#PWD.
PWD="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\Spongy\Anti_phage_systems\PADLOC\\"


#######
#Read PADLOC output, group by sample, by defense system type.
#######

def Read_PADLOC_output_group_norm_data(pwd):
    
    #Retrive data for contigs and assembly: coverage depth, length, number of reads.
    Sample_cont_cov_dict, Sample_read_num, Sample_assembly_len=read_cov_data_file('C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\Spongy\Anti_phage_systems\\')    
    
    #Prepare dataframe.
    Samples_list=['MW_2018', 'H_panicea_2018', 'H_sitiens_2018', 'I_palmata_2018']
    DS_list=['Abi', 'BREX', 'Cas', 'DISARM', 'Dnd', 'Druantia', 'gabija', 'Hachiman', 'kiwa', 'Lamassu', 'PAGOS', 'RM', 'Septu', 'shedu', 'TA', 'Thoeris', 'Wadjet', 'Zorya', 'AVAST', 'CBASS', 'dGTPase', 'Nhi', 'Gao', 'Lit', 'Retron', 'Rst', 'DRT', 'Dsr', 'Viperin', 'PrrC', 'PARIS', 'BtsA', 'dCTPdeaminase', 'ietAS', 'mza', 'qatABCD', 'upx']
    
    Samples_DS_dict_genes={}
    
    for sample_type in Samples_list:
        Samples_DS_dict_genes[sample_type]={}
        for DS_type in DS_list:
            Samples_DS_dict_genes[sample_type][DS_type]=0
    
    Samples_DS_dict_systems=copy.deepcopy(Samples_DS_dict_genes)
    Samples_DS_dict_genes_norm=copy.deepcopy(Samples_DS_dict_genes)
    Samples_DS_dict_systems_norm=copy.deepcopy(Samples_DS_dict_genes) 
    
    #Scaling constants.
    Assembly_len_const=10e6
    Number_of_reads_constant=10e6    
    
    #Retrive PADLOC output.
    for PADLOC_file in os.listdir(f'{pwd}Bacterial_fraction_csv\\'):
        sample_name=PADLOC_file.split('__')[0]
        PADLOC_output=open(f'{pwd}Bacterial_fraction_csv\\{PADLOC_file}', 'r')
        print(sample_name)
        if 'System_ID' in locals():
            del System_ID
        for line in PADLOC_output:
            line=line.rstrip().split(',')
            if line[0]!='system.number':
                if 'System_ID' in locals():
                    if System_ID==line[0]:
                        Num_genes+=1
                    else:
                        #Raw counting of DS and DS genes.
                        Samples_DS_dict_genes[sample_name][System_type]+=Num_genes
                        Samples_DS_dict_systems[sample_name][System_type]+=1
                    
                        #Normalized counting of DS and DS genes.
                        Cont_cov_depth=Sample_cont_cov_dict[sample_name][Contig_name][2]
                        Assembly_len=Sample_assembly_len[sample_name][-1]
                        Sample_reads_num=Sample_read_num[sample_name]
                    
                        Samples_DS_dict_genes_norm[sample_name][System_type]+=(Num_genes*Cont_cov_depth)*(Assembly_len_const/Assembly_len)*(Number_of_reads_constant/Sample_reads_num)
                        Samples_DS_dict_systems_norm[sample_name][System_type]+=(1*Cont_cov_depth)*(Assembly_len_const/Assembly_len)*(Number_of_reads_constant/Sample_reads_num)   
                        
                        System_ID=line[0]
                        Num_genes=1
                        Contig_name=line[1]
                        System_type=line[2] 
                        if System_type.split('_')[0]=='zorya':
                            System_type='Zorya'
                        elif System_type.split('_')[0]=='GAO':
                            System_type='Gao'
                        elif System_type.split('_')[0]=='cas':
                            System_type='Cas'  
                        elif System_type.split('_')[0]=='hachiman':
                            System_type='Hachiman' 
                        elif System_type.split('_')[0]=='cbass':
                            System_type='CBASS'  
                        elif System_type.split('_')[0]=='AVAST':
                            System_type='AVAST'   
                        elif System_type.split('_')[0]=='druantia':
                            System_type='Druantia'  
                        elif System_type.split('_')[0]=='lamassu':
                            System_type='Lamassu'   
                        elif System_type.split('_')[0]=='DRT':
                            System_type='DRT'  
                        elif System_type.split('_')[0]=='qatABCD':
                            System_type='qatABCD'   
                        elif System_type[:3]=='dsr':
                            System_type='Dsr' 
                        elif System_type.split('_')[0]=='wadjet':
                            System_type='Wadjet'     
                        elif System_type.split('_')[0]=='thoeris':
                            System_type='Thoeris'  
                        elif System_type.split('_')[0]=='septu':
                            System_type='Septu'  
                        elif System_type.split('_')[0]=='mza':
                            System_type='mza'                        
                        elif System_type[:3]=='Rst':
                            System_type='Rst'                        
                        
                else:
                    System_ID=line[0]
                    Num_genes=1
                    Contig_name=line[1]
                    System_type=line[2]
                    if System_type.split('_')[0]=='zorya':
                        System_type='Zorya'
                    elif System_type.split('_')[0]=='GAO':
                        System_type='Gao'
                    elif System_type.split('_')[0]=='cas':
                        System_type='Cas' 
                    elif System_type.split('_')[0]=='hachiman':
                        System_type='Hachiman'  
                    elif System_type.split('_')[0]=='cbass':
                        System_type='CBASS'   
                    elif System_type.split('_')[0]=='AVAST':
                        System_type='AVAST'  
                    elif System_type.split('_')[0]=='druantia':
                        System_type='Druantia'    
                    elif System_type.split('_')[0]=='lamassu':
                        System_type='Lamassu'  
                    elif System_type.split('_')[0]=='DRT':
                        System_type='DRT'
                    elif System_type.split('_')[0]=='qatABCD':
                        System_type='qatABCD'   
                    elif System_type[:3]=='dsr':
                        System_type='Dsr'  
                    elif System_type.split('_')[0]=='wadjet':
                        System_type='Wadjet'  
                    elif System_type.split('_')[0]=='thoeris':
                        System_type='Thoeris'   
                    elif System_type.split('_')[0]=='septu':
                        System_type='Septu'  
                    elif System_type.split('_')[0]=='mza':
                        System_type='mza'                    
                    elif System_type[:3]=='Rst':
                        System_type='Rst'                    
        PADLOC_output.close()
    
    
    #Prepare dataframes for plotting.
    Samples_DS_genes_df=pd.DataFrame(Samples_DS_dict_genes)
    Samples_DS_sys_df=pd.DataFrame(Samples_DS_dict_systems)
    Samples_DS_genes_df_norm=pd.DataFrame(Samples_DS_dict_genes_norm)
    Samples_DS_sys_df_norm=pd.DataFrame(Samples_DS_dict_systems_norm)  
    
    #Visualize raw abundance of DS genes and DS.
    #Polishing parameters.
    Adj_params=[0.45, -0.58]    
    heatmap_viz(Samples_DS_genes_df, pwd+'\\2018_PADLOC_Number_of_defense_systems_genes_in_datasets.svg', 'Number of defense systems genes in datasets', '\ndefense genes', Adj_params)
    heatmap_viz(Samples_DS_sys_df, pwd+'\\2018_PADLOC_Number_of_defense_systems_in_datasets.svg', 'Number of defense systems in datasets', '\ndefense systems', Adj_params)
    #Visualize normalized abundance of DS genes and DS.
    heatmap_viz(Samples_DS_genes_df_norm, pwd+'\\2018_PADLOC_Normalized_number_of_defense_systems_genes_in_datasets.svg', 'Normalized number of defense systems genes in datasets', '\ndefense genes, norm', Adj_params)
    heatmap_viz(Samples_DS_sys_df_norm, pwd+'\\2018_PADLOC_Normalized_number_of_defense_systems_in_datasets.svg', 'Normalized number of defense systems in datasets', '\ndefense systems, norm', Adj_params)        
    
    return
    

Read_PADLOC_output_group_norm_data(PWD)