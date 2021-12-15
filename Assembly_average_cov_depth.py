###############################################
##Dmitry Sutormin, 2021##
##Average coverage depth calculation for assemblies##

# Calculates average coverage depth for de novo metagenome assemblies.
###############################################

import os


#######
#Packages to be imported.
#######

# File with average coverage for contigs (output of samtools depth).
Depth_for_contigs_path="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\Spongy\Anti_phage_systems\Defense_systems\Contigs_coverage_5000_2018\\"

Depthfiles_list=os.listdir(Depth_for_contigs_path)

print('File_name\tAssembly_total_len\tTotal_reads_len\tAverage_coverage_depth')

for file_name in Depthfiles_list:
    Depthfile_full_path=os.path.join(Depth_for_contigs_path, file_name)
    Depthfile_full=open(Depthfile_full_path, 'r')
    Full_assembly_len=0
    Full_reads_len=0
    for line in Depthfile_full:
        line=line.rstrip().split('\t')
        if line[0] not in ['#rname']:
            contig_len=int(line[2])
            contig_cov_depth=float(line[6])
            contig_sum_reads_len=contig_len*contig_cov_depth
            Full_assembly_len+=contig_len
            Full_reads_len+=contig_sum_reads_len
    Average_cov_depth=float(Full_reads_len)/Full_assembly_len
    print(f'{file_name}\t{Full_assembly_len}\t{Full_reads_len}\t{Average_cov_depth}')