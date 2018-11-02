#'
#' Simple pipeline that produces a DPClust master file
#' 
#' v1.0 - 2018-11-02 - sd11 [at] sanger.ac.uk

library(dpclust3p)
library(optparse)

option_list = list(
  make_option(c("-s", "--samplenames"), type="character", default=NULL, help="Comma separated list of samplenames", metavar="character"),
  make_option(c("-d", "--donornames"), type="character", default=NULL, help="Comma separated list of donornames (use the same donor name to match multiple samples)", metavar="character"),
  make_option(c("-p", "--purities"), type="character", default=NULL, help="Comma separated list of purities (supply -p or -r)", metavar="character"),
  make_option(c("-r", "--rho_and_psi"), type="character", default=NULL, help="List of Battenberg rho and psi output files (supply -p or -r)", metavar="character"),
  make_option(c("--sex"), type="character", default=NULL, help="Comma separated list of sex of the sample", metavar="character"),
  make_option(c("-o", "--output"), type="character", default=NULL, help="Output file", metavar="character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

donornames = unlist(strsplit(opt$donornames, ","))
samplenames = unlist(strsplit(opt$samplenames, ","))
sex = unlist(strsplit(opt$sex, ","))
outputfile = opt$output

if (!is.null(opt$rho_and_psi)) {
  rho_and_psi_files = strsplit(opt$rho_and_psi, ",")
} else {
  rho_and_psi_files = NULL 
}
if (!is.null(opt$purities)) {
  purities = strsplit(opt$purities, ",")
} else {
  purities = NULL
}

if (!all(sex %in% c("male", "female"))) {
  stop("Provide male or female as sex")
}

createProjectFile(outputfile=outputfile, 
                  donornames=donornames, 
                  samplenames=samplenames, 
                  sex=sex, 
                  purities=purities, 
                  rho_and_psi_files=rho_and_psi_files)