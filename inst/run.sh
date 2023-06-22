######################################
#whole workkflow of PhaDPClust
#Haoqi, 5th Nov, 2022
#R Version 3.6.0
######################################


# Setp 1
# Recount SNVs from BAM file by RSamtools
# Input: SNV loci file; BAM file
# Script: Recount.r 
# Loci.file: four columns: CHR, POS, WT, MUT 
#            can be generated from VCF file
# Chr ranges from 1:22; ranges from chr1:chr22 if using hg38
# Run under R/3.6.0
# For some samples, there may be ERRORs:
#      ERROR : subscript contains out-of-bounds indices 
#      This script still works

mkdir -p 00_Count_Rsamtools
mkdir -p 00_Count_Rsamtools/${samplename}
R --vanilla --slave -q -f Rsamtools_reads.R --args ${samplename} ${bam} ${loci_file} ${chr} ${outdir}

# Example(hg37): 
# 		R --vanilla --slave -q -f ../Rsamtools_reads.R --args A22-I ./A22-I.bam ./A22_loci.txt 1 ./00_Count_Rsamtools/A22_I
# Example(hg38): 
# 		R --vanilla --slave -q -f ../Rsamtools_reads.R --args A22-I ./A22-I.bam ./A22_loci.txt chr1 ./00_Count_Rsamtools/A22_I


# Step 2 
# Gnenerating DPClust input
# Input: RSamtoolsCount from Step 1
# count_path: output from Step 1 
# Battenberg_Path: path for ${subclonal_file} and  ${rho_psi_file} 


# Single sample 
R --vanilla --slave -q -f Generate_singlesample.R --args ${samplename} \
															${count_path} \
															${Battenberg_Path} \ \
															${outdir} 

# Example
# R --vanilla --slave -q -f Generate_singlesample.R --args A22-I 00_Count_Rsamtools/ \
#    Battenberg_Output \
#    ./ \


# Multiple sample 
R --vanilla --slave -q -f Generate_multisample.R --args ${samplenames} \
															${count_path} \
															${Battenberg_Path} \ \
															${outdir} 

# Example
# R --vanilla --slave -q -f Generate_singlesample.R --args c("A22-I","A22-H") 00_Count_Rsamtools/ \
#    Battenberg_Output \
#    ./ \

###################
# Step 3
# Find all phased pairs of mutations from BAM files
# Required: fai_file; ign_file (file path should be changed in Phasing.r)
# Input: SNV loci file (or VCF file); BAM file; BAI file
# Script: Phasing_locifile.r / Phasing_vcf.r
# Default setting: threshold of distance is 700 (line 20)
# hg19.new.fa.fai: for chromosome in 1:22
# hg19.fa.fai: for chromosome in chr1:chr22
mkdir -p 03_Phasing
R --vanilla --slave -q -f ./Phasing_locifile.R --args ${samplename} ${bam} ${loci_file} ${outdir}
R --vanilla --slave -q -f ./Phasing_vcf.R --args ${samplename} ${bam} ${vcf} ${outdir}

# Example
# R --vanilla --slave -q -f ./Phasing_locifile.R --args A22-I \
# ./A22-I.bam \
# ./A22_loci.txt \
# ./03_Phasing \
# /hg19.new.fa.fai \  #new because chr is 1:22 rather than chr1:chr22
# /ignore_dpclust_preprocessing.txt 


###################
# Step 4
# Extract phased pairs from specific copy number states
# Currently only use pairs from diplodi and haploidy region 
# Required: dpclust input (end with "_ssDPI.txt")
# Input: SNV loci file (or VCF file); BAM file;
# Script: Generate_Phasing.r
# Output: clphasing file (clphasing represent clonal phasing, that there is no subclonal gain or loss)

# The last three number represent the copy.major, copy.minor, copy.frac
# For muts from diploid region, they shoud be 1,1,1
# For muts for haploid region, they should be 1,0,1
# For muts from both diploid and haploid, they should be 1, 1

mkdir -p 04_cl_phasing
R --vanilla --slave -q -f generate_clphasing.r --args ${samplename} \
${Phasing_dir} \
${Battenberg_Path} \
${DPClust_input_path} \
${outdir} ${copy_major} ${copy_frac} ${copy_minor}

# Example
# R --vanilla --slave -q -f generate_clphasing.R --args A22-I ./ ./ ./ ./04_cl_phasing 1 1
# R --vanilla --slave -q -f generate_clphasing.R --args A22-I ./ ./ ./ ./04_cl_phasing 1 1 1
# R --vanilla --slave -q -f generate_clphasing.R --args A22-I ./ ./ ./ ./04_cl_phasing 1 1 0


###################
# Step 5
# Run PhaDPClust
# Required file: dataset.RData,clphasing.RData,DPClust input file
#               function.R
# Script: workflow_single.R
# Num_subsamples: number of samples used, 10 at default
# killclu: if cluster contains less than this threshold of clusters, it would be killed; 0.05 at default

# move cl_phasing.RData and ssDPI.txt into one folder 
mkdir -p 05_Data
cp 02_Input_all_SNVs/${samplename}_ssDPI.txt ./05_Data/${samplename}_ssDPI.txt
cp 04_cl_phasing/${samplename}_cl_phasing_pass.RData ./05_Data/${samplename}_cl_phasing_pass.RData

# run DPClust 
mkdir -p 05_PhaDPClust
mkdir -p 05_DPClust
# This line for normal DPclust
R --vanilla --slave -q -f pipeline.R --args -r 1 -d ../02_Input_all_SNVs -o ./05_DPClust -i ./dp.txt
# This line for phasing DPclust
R --vanilla --slave -q -f Pipeline_phasingversion.R --args -r 1 -d ./05_Data -o ./05_PhaDPClust -i ./${samplename}_phadp.txt -a phasing_ass


































