#'
#' Simple DPClust preprocessing pipeline that takes a VCF file with SNV calls. Then fetches
#' allele counts from the specified bam file. With the counts and the copy number fit it
#' creates a DPClust input file with all the required columns.
#' 
#' Dependencies:
#'  * The alleleCounter C utility must be in $PATH
#' 
#' v1.0 - 2015-11-04 - sd11@sanger.ac.uk

args = commandArgs(T)
samplename = toString(args[1]) # Name of the sample, used to name output files
bam_file = toString(args[2]) # Full path to the bam file, there must be a bai file in the same directory
vcf_file = toString(args[3]) # Full path to the vcf file with SNV calls. All calls in this will be used
samplestatistics.file = toString(args[4]) # Full path to a samplestatistics output file from ASCAT NGS
copynumber.caveman.file = toString(args[5]) # Full path to the copynumber.cabeman output file from ASCAT NGS
sex = toString(args[6]) # Specify male or female
output_dir = toString(args[7]) # Full path to where the output should be written
fai_file = toString(args[8]) # Full path to the reference genome index file used for this sample
ign_file = toString(args[9]) # Full path to simple list of chromosome names to ignore (must contain at least Y and MT)

library(dpclust3p)

# Define the final output file
dpoutput_file = paste(output_dir, "/", samplename, "_allDirichletProcessInfo.txt", sep="")

# Define various temp files
loci_file = paste(output_dir, "/", samplename, "_loci.txt", sep="")
allelecounts_file = paste(output_dir, "/", samplename, "_alleleFrequencies.txt", sep="")
battenberg_output_prefix = paste(output_dir, "/", samplename, sep="")
subclones_file = paste(battenberg_output_prefix, "_subclones.txt", sep="")
rho_and_psi_file = paste(battenberg_output_prefix, "_rho_and_psi.txt", sep="")

# Dump loci - this function can take multiple vcf files when multiple samples from same donor
vcf2loci(vcf_file=vcf_file, fai_file=fai_file, ign_file=ign_file, outfile=loci_file)

# Fetch allele counts
alleleCount(locifile=loci_file, bam=bam_file, outfile=allelecounts_file, min_baq=20, min_maq=35)

# Transform ASCAT NGS output into Battenberg output
ascatNgsToBattenberg(outfile.prefix=battenberg_output_prefix, copynumber.caveman.file=copynumber.caveman.file, samplestatistics.file=samplestatistics.file)

# Create dpIn file
runGetDirichletProcessInfo(loci_file=loci_file, 
                           allele_frequencies_file=allelecounts_file, 
                           cellularity_file=rho_and_psi_file, 
                           subclone_file=subclones_file, 
                           gender=sex, 
                           SNP.phase.file="NA", 
                           mut.phase.file="NA", 
                           output_file=dpoutput_file)
