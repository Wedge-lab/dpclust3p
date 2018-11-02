#'
#' Simple DPClust preprocessing pipeline that takes a VCF file with SNV calls. It fetches
#' allele counts from the specified VCF file. With the counts and the copy number fit it
#' creates a DPClust input file with all the required columns.
#' 
#' Dependencies:
#'  * The alleleCounter C utility must be in $PATH
#' 
#' v1.1 - 2018-11-01 - sd11 [at] sanger.ac.uk

library(dpclust3p)
library(optparse)

option_list = list(
  make_option(c("-s", "--samplename"), type="character", default=NULL, help="Samplename of the sample to run", metavar="character"),
  make_option(c("-v", "--vcf"), type="character", default=NULL, help="VCF file with mutation data", metavar="character"),
  make_option(c("-r", "--rho_and_psi"), type="character", default=NULL, help="Battenberg rho and psi output file", metavar="character"),
  make_option(c("-c", "--copynumber"), type="character", default=NULL, help="Battenberg subclones output file with copy number data", metavar="character"),
  make_option(c("--sex"), type="character", default=NULL, help="Sex of the sample", metavar="character"),
  make_option(c("-o", "--output"), type="character", default=NULL, help="Output directory", metavar="character"),
  make_option(c("--fai"), type="character", default=NULL, help="Reference genome index", metavar="character"),
  make_option(c("--ign"), type="character", default=NULL, help="File with a list of contigs to ignore", metavar="character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

samplename = opt$samplename
vcf_file = opt$vcf
rho_and_psi_file = opt$rho_and_psi
subclones_file = opt$copynumber
sex = opt$sex
output_dir = opt$output
fai_file = opt$fai
ign_file = opt$ign

.checkfile = function(infile) {
  if (!file.exists(infile)) {
    stop(paste("File", infile, "does not exist", sep=""))
  }
}

.checkfile(vcf_file)
.checkfile(rho_and_psi_file)
.checkfile(subclones_file)
.checkfile(fai_file)
.checkfile(ign_file)

if (!sex %in% c("male", "female")) {
  stop("Provide male or female as sex")
}

# Define the final output file
dpoutput_file = file.path(output_dir, paste(samplename, "_allDirichletProcessInfo.txt", sep=""))

# Define various temp files
loci_file = file.path(output_dir, paste(samplename, "_loci.txt", sep=""))
allelecounts_file = file.path(output_dir, paste(samplename, "_alleleFrequencies.txt", sep=""))

# Dump loci - this function can take multiple vcf files when multiple samples from same donor
vcf2loci(vcf_file=vcf_file, fai_file=fai_file, ign_file=ign_file, outfile=loci_file)

# Fetch allele counts
dumpCounts.Sanger(vcf_infile=vcf_file, tumour_outfile=allelecounts_file)

# Create dpIn file
runGetDirichletProcessInfo(loci_file=loci_file, 
                           allele_frequencies_file=allelecounts_file, 
                           cellularity_file=rho_and_psi_file, 
                           subclone_file=subclones_file, 
                           gender=sex, 
                           SNP.phase.file="NA", 
                           mut.phase.file="NA", 
                           output_file=dpoutput_file)

# Cleanup
file.remove(loci_file)
file.remove(allelecounts_file)
