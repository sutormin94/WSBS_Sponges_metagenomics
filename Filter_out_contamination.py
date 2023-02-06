###############################################
##Dmitry Sutormin, 2021##
##Filters out contamination from assemblies, as reported by NCBI##

# Filters out contigs reported as contamination by NCBI.
###############################################

import os
import Bio
from Bio import SeqIO

# Path to_folder with assemblies.
Assemblies_path="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\Spongy\Metagenome_assembly\Assemblies\Contigs_2018_200bp_unfilt\\"

# Contamination files path.
Contam_path="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\Spongy\Metagenome_assembly\Assemblies\Contamination_reports_NCBI\\"

# Filtered assemblies path out.
Filt_assemblies_path="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\Spongy\Metagenome_assembly\Assemblies\Contigs_2018_200bp_filt\\"

if not os.path.isdir(Filt_assemblies_path):
    os.mkdir(Filt_assemblies_path)
    

def read_contam(inpath):
    filein=open(inpath, 'r')
    Contigs_to_remove=[]
    for line in filein:
        if line[:4]=="NODE":
            Contig_ID=line.split('\t')[0]
            print(Contig_ID)
            Contigs_to_remove.append(Contig_ID)
        elif line[:4]=="lcl|":
            Contig_ID=line.split(' ')[0].replace('lcl|', '')
            print('Here!', Contig_ID)
            Contigs_to_remove.append(Contig_ID)
    
    filein.close()
    return Contigs_to_remove


def read_assemblies_filt_write(inpath, outpath, Contigs_to_remove):
    filein=open(inpath, 'r')
    fileout=open(outpath, "w+")
        
    for record in SeqIO.parse(filein, "fasta"):
        Contig_ID=record.id
        if Contig_ID not in Contigs_to_remove:
            fileout.write(f'>{Contig_ID}\n{str(record.seq)}\n')
        else:
            print(Contig_ID)
            
    filein.close()
    fileout.close()
    return
    
   
    
Contamfiles_list=os.listdir(Contam_path)

for Contamfile in Contamfiles_list:
    Basename=Contamfile.replace('.txt', '')
    print(Basename)
    Contigs_to_remove=read_contam(os.path.join(Contam_path, Contamfile))
    read_assemblies_filt_write(os.path.join(Assemblies_path,Basename+'.fasta'), 
                               os.path.join(Filt_assemblies_path,Basename+'.fasta'), 
                               Contigs_to_remove)
    