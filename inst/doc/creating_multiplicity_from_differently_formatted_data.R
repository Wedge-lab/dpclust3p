## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
devtools::install_github("Wedge-lab/dpclust3p", ref = "hotfix_if_is_na", dependencies = TRUE)
library(dpclust3p)

## -----------------------------------------------------------------------------
pcawg_muts <- read.delim("https://pcawg-hub.s3.us-east-1.amazonaws.com/download/October_2016_all_patients_2778.snv_mnv_indel.maf.xena.nonUS", sep = "\t", header = TRUE, nrows = 20000)

## -----------------------------------------------------------------------------
pcawg_cn <- read.delim("https://pcawg-hub.s3.us-east-1.amazonaws.com/download/20170119_final_consensus_copynumber_sp", sep = "\t", header = TRUE)

## -----------------------------------------------------------------------------
pcawg_purity <- read.delim("https://pcawg-hub.s3.us-east-1.amazonaws.com/download/consensus.20170217.purity.ploidy_sp", sep = "\t", header = TRUE)

## -----------------------------------------------------------------------------
pcawg_clinical <- read.delim("https://pcawg-hub.s3.us-east-1.amazonaws.com/download/pcawg_donor_clinical_August2016_v9_sp", sep = "\t", header = TRUE)

## -----------------------------------------------------------------------------
if(file.exists("./pseudo_dpclust_files")){
} else {
  dir.create("./pseudo_dpclust_files")
}
if(file.exists("./pseudo_dpclust_files/pseudo_allele_frequency/")){
} else {
  dir.create("./pseudo_dpclust_files/pseudo_allele_frequency/")
}
if(file.exists("./pseudo_dpclust_files/pseudo_cellularity/")){
} else {
  dir.create("./pseudo_dpclust_files/pseudo_cellularity/")
}
if(file.exists("./pseudo_dpclust_files/pseudo_loci_file/")){
} else {
  dir.create("./pseudo_dpclust_files/pseudo_loci_file/")
}
if(file.exists("./pseudo_dpclust_files/pseudo_subclones/")){
} else {
  dir.create("./pseudo_dpclust_files/pseudo_subclones/")
}
if(file.exists("./dpclust_info")){
} else {
  dir.create("./dpclust_info")
}

## -----------------------------------------------------------------------------
sample_1_muts <- subset(pcawg_muts, Sample == pcawg_muts$Sample[1])
sample_1_muts <- subset(sample_1_muts, reference %in% c("A","T","G","C") & alt %in% c("A","T","G","C"))

## -----------------------------------------------------------------------------
sample_1_cn_pcawg <- subset(pcawg_cn, sampleID == pcawg_muts$Sample[1])
sample_1_cn <- sample_1_cn_pcawg[,c("chr","start","end")]
colnames(sample_1_cn) <- c("chr","startpos","endpos")
sample_1_cn$BAF <- NA
sample_1_cn$pval <- NA
sample_1_cn$LogR	 <- NA
sample_1_cn$ntot	 <- sample_1_cn_pcawg$total_cn
sample_1_cn$nMaj1_A	 <- sample_1_cn_pcawg$major_cn
sample_1_cn$nMin1_A	 <- sample_1_cn_pcawg$minor_cn
sample_1_cn$frac1_A	 <- 1
sample_1_cn$nMaj2_A	 <- NA
sample_1_cn$nMin2_A	 <- NA
sample_1_cn$frac2_A	 <- NA
sample_1_cn$SDfrac_A	 <- NA
sample_1_cn$SDfrac_A_BS	 <- NA
sample_1_cn$frac1_A_0_025	 <- NA
sample_1_cn$frac1_A_0_975	 <- NA
sample_1_cn$p_maj	 <- NA
sample_1_cn$p_min	 <- NA
sample_1_cn$flag <- NA
# chr	startpos	endpos	BAF	pval	LogR	ntot	nMaj1_A	nMin1_A	frac1_A	nMaj2_A	nMin2_A	frac2_A	SDfrac_A	SDfrac_A_BS	frac1_A_0_025	frac1_A_0_975	p_maj	p_min	flag
head(sample_1_cn)

## -----------------------------------------------------------------------------
write.table(sample_1_cn, "./pseudo_dpclust_files/pseudo_subclones/sample_1_subclones.txt", col.names = TRUE, quote = FALSE, row.names = FALSE, sep = "\t")

## -----------------------------------------------------------------------------
sample_1_purity <- subset(pcawg_purity, samplename == pcawg_muts$Sample[1])

## -----------------------------------------------------------------------------
sample_1_clinical <- subset(pcawg_clinical, xena_sample == pcawg_muts$Sample[1])

## -----------------------------------------------------------------------------
sample_1_muts$Count_A[sample_1_muts$reference == "A"] <- sample_1_muts$t_ref_count[sample_1_muts$reference == "A"]
sample_1_muts$Count_C[sample_1_muts$reference == "C"] <- sample_1_muts$t_ref_count[sample_1_muts$reference == "C"]
sample_1_muts$Count_G[sample_1_muts$reference == "G"] <- sample_1_muts$t_ref_count[sample_1_muts$reference == "G"]
sample_1_muts$Count_T[sample_1_muts$reference == "T"] <- sample_1_muts$t_ref_count[sample_1_muts$reference == "T"]
sample_1_muts$Count_A[sample_1_muts$alt == "A"] <- sample_1_muts$t_alt_count[sample_1_muts$alt == "A"]
sample_1_muts$Count_C[sample_1_muts$alt == "C"] <- sample_1_muts$t_alt_count[sample_1_muts$alt == "C"]
sample_1_muts$Count_G[sample_1_muts$alt == "G"] <- sample_1_muts$t_alt_count[sample_1_muts$alt == "G"]
sample_1_muts$Count_T[sample_1_muts$alt == "T"] <- sample_1_muts$t_alt_count[sample_1_muts$alt == "T"]
sample_1_muts$Total_depth <- rowSums(sample_1_muts[,c("Count_A","Count_C","Count_G","Count_T")], na.rm = TRUE)

sample_1_alleles <- (sample_1_muts[,c("chr","start","Count_A","Count_C","Count_G","Count_T","Total_depth")])
sample_1_alleles[is.na(sample_1_alleles)] <- 0
head(sample_1_alleles)

## -----------------------------------------------------------------------------
colnames(sample_1_alleles) <- c("CHR","START","Count_A","Count_C","Count_G","Count_T","Total_depth")

## -----------------------------------------------------------------------------
write.table(sample_1_alleles, "./pseudo_dpclust_files/pseudo_allele_frequency/sample_1_pseudo_allele_frequency.tsv", col.names = TRUE, quote = FALSE, row.names = FALSE, sep = "\t")

## -----------------------------------------------------------------------------
sample_1_loci <- sample_1_muts[,c("chr","start","reference","alt")]
colnames(sample_1_loci) <- c("CHR","START","REF","ALT")

## -----------------------------------------------------------------------------
write.table(sample_1_loci, "./pseudo_dpclust_files/pseudo_loci_file/sample_1_pseudo_loci.txt", col.names = FALSE, quote = FALSE, row.names = FALSE, sep = "\t")

## -----------------------------------------------------------------------------
sample_1_pseudo_cellularity <- matrix(nrow = 3, ncol = 4, dimnames = list(c("ASCAT","FRAC_GENOME","REF_SEG"),c("rho","psi","distance","is.best")),
                             data = c(NA, NA, NA, NA, sample_1_purity$purity, sample_1_purity$ploidy, NA, "TRUE", NA, NA, 0, "FALSE"), byrow = TRUE)

## -----------------------------------------------------------------------------
write.table(sample_1_pseudo_cellularity, "./pseudo_dpclust_files/pseudo_cellularity/sample_1_pseudo_cellularity.tsv", col.names = TRUE, quote = FALSE, row.names = TRUE, sep = "\t")

## -----------------------------------------------------------------------------
runGetDirichletProcessInfo(loci_file = "./pseudo_dpclust_files/pseudo_loci_file/sample_1_pseudo_loci.txt",
                           allele_frequencies_file = "./pseudo_dpclust_files/pseudo_allele_frequency/sample_1_pseudo_allele_frequency.tsv",
                           cellularity_file = "./pseudo_dpclust_files/pseudo_cellularity/sample_1_pseudo_cellularity.tsv",
                           subclone_file = "./pseudo_dpclust_files/pseudo_subclones/sample_1_subclones.txt",
                           gender = sample_1_clinical$donor_sex,
                           SNP.phase.file = "NA", mut.phase.file = "NA",
                           output_file = "./dpclust_info/sample_1_dpclust_info")

## -----------------------------------------------------------------------------
sample_1_dpclust_info <- read.delim("./dpclust_info/sample_1_dpclust_info", sep = "\t", header = TRUE)

## -----------------------------------------------------------------------------
head(sample_1_dpclust_info)

