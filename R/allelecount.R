################################################
# Allele counting functions - SNV, indel and SV
################################################
#' Run alleleCount
#' 
#' Count the alleles for specified locations in the loci file. Expects alleleCount binary in $PATH
#' @param locifile A file with at least chromsome and position columns of the locations to be counted
#' @param bam A bam file
#' @param outfile Where to store the output
#' @param min_baq The minimum base quality required for a read to be counted
#' @param min_maq The minimum mapping quality required for a read to be counted
#' @author sd11
#' @export
alleleCount = function(locifile, bam, outfile, min_baq=20, min_maq=35) {
  cmd = paste(ALLELECOUNTER,
              "-b", bam,
              "-o", outfile,
              "-l", locifile,
              "-m", min_baq,
              "-q", min_maq, sep=" ")
  system(cmd, wait=T)
}

#' Obtain alt allele for all variants as a vector
#' 
#' Dumps the alt allele in a vector. It takes only the first element if multiple alt alleles are listed
#' @param vcf A VCF object
#' @return A vector with the alt allele
#' @author sd11
#' @export
getAltAllele = function(vcf) {
  clist = CharacterList(VariantAnnotation::alt(vcf))
  mult = elementNROWS(clist) > 1L
  
  if (any(mult)) {
    warning(paste("Took first alt allele only for these variants: ", paste(seqnames(vcf), start(vcf))))
  }
  
  clist[mult] = lapply(clist[mult], paste0, collapse=",")
  return(unlist(clist))
}

#' Dump allele counts from vcf - Sanger ICGC pancancer pipeline
#'
#' Dump allele counts stored in the sample columns of the VCF file. Output will go into a file
#' supplied as tumour_outfile and optionally normal_outfile. It will be a fully formatted
#' allele counts file as returned by alleleCounter.
#' @param vcf_infile The vcf file to read in
#' @param tumour_outfile File to save the tumour counts to
#' @param normal_outfile Optional parameter specifying where the normal output should go
#' @param refence_genome Optional parameter specifying the reference genome build used
#' @param samplename Optional parameter specifying the samplename to be used for matching the right column in the VCF
#' @author sd11
#' @export
dumpCounts.Sanger = function(vcf_infile, tumour_outfile, normal_outfile=NA, refence_genome="hg19", samplename=NA) {
  dumpCountsFromVcf(vcf_infile, tumour_outfile, centre="sanger", normal_outfile=normal_outfile, refence_genome=refence_genome, samplename=samplename)
}

#' Dump allele counts from vcf - DKFZ ICGC pancancer pipeline
#'
#' Dump allele counts stored in the sample columns of the VCF file. Output will go into a file
#' supplied as tumour_outfile and optionally normal_outfile. It will be a fully formatted
#' allele counts file as returned by alleleCounter.
#' @param vcf_infile The vcf file to read in
#' @param tumour_outfile File to save the tumour counts to
#' @param normal_outfile Optional parameter specifying where the normal output should go
#' @param refence_genome Optional parameter specifying the reference genome build used
#' @param samplename Optional parameter specifying the samplename to be used for matching the right column in the VCF
#' @author sd11
#' @export
dumpCounts.DKFZ = function(vcf_infile, tumour_outfile, normal_outfile=NA, refence_genome="hg19", samplename=NA) {
  dumpCountsFromVcf(vcf_infile, tumour_outfile, centre="dkfz", normal_outfile=normal_outfile, refence_genome=refence_genome, samplename=samplename)
}

#' Dump allele counts from vcf - Broad ICGC pancancer pipeline
#'
#' Dump allele counts stored in the sample columns of the VCF file. Output will go into a file
#' supplied as tumour_outfile and optionally normal_outfile. It will be a fully formatted
#' allele counts file as returned by alleleCounter.
#' @param vcf_infile The vcf file to read in
#' @param tumour_outfile File to save the tumour counts to
#' @param normal_outfile Optional parameter specifying where the normal output should go
#' @param refence_genome Optional parameter specifying the reference genome build used
#' @param samplename Optional parameter specifying the samplename to be used for matching the right column in the VCF
#' @author sd11
#' @export
dumpCounts.Broad = function(vcf_infile, tumour_outfile, normal_outfile=NA, refence_genome="hg19", samplename=NA) {
  dumpCountsFromVcf(vcf_infile, tumour_outfile, centre="broad", normal_outfile=normal_outfile, refence_genome=refence_genome, samplename=samplename)
}

#' Dump allele counts from vcf - MuSE ICGC pancancer pipeline
#'
#' Dump allele counts stored in the sample columns of the VCF file. Output will go into a file
#' supplied as tumour_outfile and optionally normal_outfile. It will be a fully formatted
#' allele counts file as returned by alleleCounter.
#' @param vcf_infile The vcf file to read in
#' @param tumour_outfile File to save the tumour counts to
#' @param normal_outfile Optional parameter specifying where the normal output should go
#' @param refence_genome Optional parameter specifying the reference genome build used
#' @param samplename Optional parameter specifying the samplename to be used for matching the right column in the VCF
#' @author sd11
#' @export
dumpCounts.Muse = function(vcf_infile, tumour_outfile, normal_outfile=NA, refence_genome="hg19", samplename=NA) {
  dumpCountsFromVcf(vcf_infile, tumour_outfile, centre="muse", normal_outfile=normal_outfile, refence_genome=refence_genome, samplename=samplename)
}

#' Dump allele counts from vcf - ICGC pancancer consensus SNV pipeline
#'
#' Dump allele counts stored in the info column of the VCF file. Output will go into a file
#' supplied as tumour_outfile. It will be a fully formatted allele counts file as returned 
#' by alleleCounter. There are no counts for the matched normal. 
#' @param vcf_infile The vcf file to read in
#' @param tumour_outfile File to save the tumour counts to
#' @param normal_outfile Optional parameter specifying where the normal output should go
#' @param refence_genome Optional parameter specifying the reference genome build used
#' @param samplename Optional parameter specifying the samplename to be used for matching the right column in the VCF
#' @author sd11
#' @export
dumpCounts.ICGCconsensusSNV = function(vcf_infile, tumour_outfile, normal_outfile=NA, refence_genome="hg19", samplename=NA) {
  dumpCountsFromVcf(vcf_infile, tumour_outfile, centre="icgc_consensus_snv", normal_outfile=normal_outfile, refence_genome=refence_genome, samplename=samplename)
}

#' Dump allele counts from vcf - ICGC pancancer consensus indel pipeline
#'
#' Dump allele counts stored in the info column of the VCF file. Output will go into a file
#' supplied as tumour_outfile. It will be a fully formatted allele counts file as returned 
#' by alleleCounter. There are no counts for the matched normal. 
#' @param vcf_infile The vcf file to read in
#' @param tumour_outfile File to save the tumour counts to
#' @param refence_genome Optional parameter specifying the reference genome build used
#' @param dummy_alt_allele Specify allele to be used to encode the alt counts (the indel can be multiple alleles)
#' @param dummy_ref_allele Specify allele to be used to encode the ref counts (the indel can be multiple alleles)
#' @author sd11
#' @export
dumpCounts.ICGCconsensusIndel = function(vcf_infile, tumour_outfile, refence_genome="hg19", dummy_alt_allele=NA, dummy_ref_allele=NA) {
  dumpCountsFromVcf(vcf_infile, tumour_outfile, centre="icgc_consensus_indel", normal_outfile=NA, refence_genome=refence_genome, samplename=NA, dummy_alt_allele=dummy_alt_allele, dummy_ref_allele=dummy_ref_allele)
}

#' Dump allele counts from vcf - Strelka indel pipeline
#'
#' Dump allele counts stored in the info column of the VCF file. Output will go into a file
#' supplied as tumour_outfile. It will be a fully formatted allele counts file as returned
#' by alleleCounter. There are no counts for the matched normal.
#' @param vcf_infile The vcf file to read in
#' @param tumour_outfile File to save the tumour counts to
#' @param refence_genome Optional parameter specifying the reference genome build used
#' @param dummy_alt_allele Specify allele to be used to encode the alt counts (the indel can be multiple alleles)
#' @param dummy_ref_allele Specify allele to be used to encode the ref counts (the indel can be multiple alleles)
#' @author sd11
#' @export
dumpCounts.StrelkaIndel = function(vcf_infile, tumour_outfile, refence_genome="hg19", dummy_alt_allele=NA, dummy_ref_allele=NA) {
	dumpCountsFromVcf(vcf_infile, tumour_outfile, centre="strelka_indel", normal_outfile=NA, refence_genome=refence_genome, samplename=NA, dummy_alt_allele=dummy_alt_allele, dummy_ref_allele=dummy_ref_allele)
}

#' Dump allele counts from vcf - CGP Pindel indel pipeline
#'
#' Dump allele counts stored in the info column of the VCF file. Output will go into a file
#' supplied as tumour_outfile. It will be a fully formatted allele counts file as returned
#' by alleleCounter. There are no counts for the matched normal.
#' @param vcf_infile The vcf file to read in
#' @param tumour_outfile File to save the tumour counts to
#' @param refence_genome Optional parameter specifying the reference genome build used
#' @param dummy_alt_allele Specify allele to be used to encode the alt counts (the indel can be multiple alleles)
#' @param dummy_ref_allele Specify allele to be used to encode the ref counts (the indel can be multiple alleles)
#' @author sd11
#' @export
dumpCounts.cgpPindel = function(vcf_infile, tumour_outfile, refence_genome="hg19", dummy_alt_allele=NA, dummy_ref_allele=NA) {
  dumpCountsFromVcf(vcf_infile, tumour_outfile, centre="cgppindel_indel", normal_outfile=NA, refence_genome=refence_genome, samplename=NA, dummy_alt_allele=dummy_alt_allele, dummy_ref_allele=dummy_ref_allele)
}

#' Dump allele counts from vcf - Mutect
#'
#' Dump allele counts stored in the info column of the VCF file. Output will go into a file
#' supplied as tumour_outfile. It will be a fully formatted allele counts file as returned 
#' by alleleCounter.
#' @param vcf_infile The vcf file to read in
#' @param tumour_outfile File to save the tumour counts to
#' @param samplename Parameter specifying the samplename to be used for matching the right column in the VCF
#' @param normal_outfile Optional parameter specifying where the normal output should go
#' @param refence_genome Optional parameter specifying the reference genome build used
#' @author sd11
#' @export
dumpCounts.mutect = function(vcf_infile, tumour_outfile, samplename, normal_outfile=NA, refence_genome="hg19") {
  dumpCountsFromVcf(vcf_infile, tumour_outfile, centre="mutect", normal_outfile=normal_outfile, refence_genome=refence_genome, samplename=samplename)
}

#' Dump allele counts from VCF
#' 
#' This function implements all the steps required for dumping counts from VCF
#' as supplied by the ICGC pancancer pipelines. 
#' @noRd
dumpCountsFromVcf = function(vcf_infile, tumour_outfile, centre, normal_outfile=NA, refence_genome="hg19", samplename=NA, dummy_alt_allele=NA, dummy_ref_allele=NA) {
  # Helper function for writing the output  
  write.output = function(output, output_file) {
    write.table(output, file=output_file, col.names=T, quote=F, row.names=F, sep="\t")
  }
  
  # Read in the vcf and dump the tumour counts in the right format
  v = VariantAnnotation::readVcf(vcf_infile, refence_genome)
  if (nrow(v) > 0) {
    tumour = getCountsTumour(v, centre=centre, samplename=samplename, dummy_alt_allele=dummy_alt_allele, dummy_ref_allele=dummy_ref_allele)
  } else {
    tumour = NA
  }
  tumour = formatOutput(tumour, v)
  write.output(tumour, tumour_outfile)
  
  # Optionally dump the normal counts in the right format
  if (!is.na(normal_outfile)) {
    if (nrow(v) > 0) {
      normal = getCountsNormal(v, centre=centre, samplename=samplename)
    } else {
      normal = NA
    }
    normal = formatOutput(normal, v)
    write.output(normal, normal_outfile)
  }
}

#' Format a 4 column counts table into the alleleCounter format. This function assumes A, C, G, T format.
#' @noRd
formatOutput = function(counts_table, v) {
  # Check if the vcf object is empty, if so, generate dummy data
  if (nrow(v)==0) {
    output = data.frame(matrix(ncol = 7, nrow = 0))
  } else {
    output = data.frame(as.character(seqnames(v)), start(ranges(v)), counts_table, rowSums(counts_table))
  }
  colnames(output) = c("#CHR","POS","Count_A","Count_C","Count_G","Count_T","Good_depth")
  return(output)
}

#' Dump allele counts from vcf for normal
#' 
#' Returns an allele counts table for the normal sample
#' @param v The vcf file
#' @param centre The sequencing centre of which pipeline the vcf file originates
#' @return An array with 4 columns: Counts for A, C, G, T
#' @author sd11
#' @noRd
getCountsNormal = function(v, centre="sanger", samplename=NA) {
  if (centre=="sanger") {
    return(getAlleleCounts.Sanger(v, 1))
  } else if(centre=="dkfz") {
    sample_col = which(colnames(VariantAnnotation::geno(v)$DP4) == "CONTROL")
    return(getAlleleCounts.DKFZ(v, sample_col))
  } else if (centre=="muse") {
    # Assuming the tumour name is provided
    sample_col = which(colnames(VariantAnnotation::geno(v)$AD) != samplename)
    return (getAlleleCounts.MuSE(v, sample_col))
  } else if (centre=="broad") {
    print("The Broad ICGC pipeline does not report allele counts for the matched normal")
    q(save="no", status=1)
  } else if (centre=="icgc_consensus_snv") {
    print("The ICGC consensus pipeline does not report allele counts for the matched normal")
    q(save="no", status=1)
  } else if (centre=="icgc_consensus_indel") {
    print("The ICGC consensus indel pipeline does not report allele counts for the matched normal")
    q(save="no", status=1)
  } else if (centre=="mutect") {
    sample_col = which(colnames(VariantAnnotation::geno(v)$AD) == samplename)
    return(getAlleleCounts.mutect(v, sample_col))
  } else {
    print(paste("Supplied centre not supported:", centre))
    q(save="no", status=1)
  }
}

#' Dump allele counts from vcf for tumour
#' 
#' Returns an allele counts table for the tumour sample
#' @param v The vcf file
#' @param centre The sequencing centre of which pipeline the vcf file originates
#' @param samplename Samplename to be used when matching sample columns (Default: NA)
#' @param dummy_alt_allele Supply if the alt allele is to be overwritten. The counts of the alt will be put in the column corresponding to this allele. (Default: NA)
#' @param dummy_ref_allele Supply if the ref allele is to be overwritten. The counts of the ref will be put in the column corresponding to this allele. (Default: NA)
#' @return An array with 4 columns: Counts for A, C, G, T
#' @author sd11
#' @noRd
getCountsTumour = function(v, centre="sanger", samplename=NA, dummy_alt_allele=NA, dummy_ref_allele=NA) {
  if (centre=="sanger") {
    return(getAlleleCounts.Sanger(v, 2))
  } else if (centre=="dkfz") {
    sample_col = which(colnames(VariantAnnotation::geno(v)$DP4) == "TUMOR")
    return(getAlleleCounts.DKFZ(v, sample_col))
  } else if (centre=="muse") {
    if (is.na(samplename)) { stop("When dumping allele counts from the Muse samplename must be supplied") }
    sample_col = which(colnames(VariantAnnotation::geno(v)$AD) == samplename)
    return (getAlleleCounts.MuSE(v, sample_col))
  } else if (centre=="broad") {
    if (is.na(samplename)) { stop("When dumping allele counts from the Broad ICGC pipeline samplename must be supplied") }
    return(getAlleleCounts.Broad(v, 1))
  } else if (centre=="icgc_consensus_snv") {
    return(getAlleleCounts.ICGC_consensus_snv(v))
  } else if (centre=="icgc_consensus_indel") {
    if (is.na(dummy_alt_allele) | is.na(dummy_ref_allele)) { stop("When dumping allele counts from the ICGC consensus indels dummy_alt_allele and dummy_ref_allele must be supplied") }
    return(getAlleleCounts.ICGC_consensus_indel(v, dummy_alt_allele=dummy_alt_allele, dummy_ref_allele=dummy_ref_allele))
  } else if (centre=="mutect") {
    if (is.na(samplename)) { stop("When dumping allele counts from the Mutect samplename must be supplied") }
    sample_col = which(colnames(VariantAnnotation::geno(v)$AD) == samplename)
    return(getAlleleCounts.mutect(v, sample_col))
  } else if (centre=="strelka_indel") {
    if (is.na(dummy_alt_allele) | is.na(dummy_ref_allele)) { stop("When dumping allele counts from the Strelka indels dummy_alt_allele and dummy_ref_allele must be supplied") }
    return(getAlleleCounts.Strelka_indel(v, dummy_alt_allele=dummy_alt_allele, dummy_ref_allele=dummy_ref_allele))
  } else if (centre=="scalpel_indel") {
    if (is.na(dummy_alt_allele) | is.na(dummy_ref_allele)) { stop("When dumping allele counts from the Strelka indels dummy_alt_allele and dummy_ref_allele must be supplied") }
    return(getAlleleCounts.Scalpel_indel(v, dummy_alt_allele=dummy_alt_allele, dummy_ref_allele=dummy_ref_allele))
  } else if (centre=="cgppindel_indel") {
    sample_col = which(colnames(geno(v)$WTR)==samplename)
    return(getAlleleCounts.cgpPindel_indel(v, sample_col=sample_col, dummy_alt_allele=dummy_alt_allele, dummy_ref_allele=dummy_ref_allele))
  } else {
    print(paste("Supplied centre not supported:", centre))
    q(save="no", status=1)
  }
}

#' Dump allele counts from Sanger pipeline vcf
#' 
#' Helper function that dumps the allele counts from a Sanger pipeline VCF file
#' @param v The vcf file
#' @param sample_col The column in which the counts are. If it's the first sample mentioned in the vcf this would be sample_col 1
#' @return An array with 4 columns: Counts for A, C, G, T
#' @author sd11
#' @noRd
getAlleleCounts.Sanger = function(v, sample_col) {
  return(cbind(VariantAnnotation::geno(v)$FAZ[,sample_col]+VariantAnnotation::geno(v)$RAZ[,sample_col], 
               VariantAnnotation::geno(v)$FCZ[,sample_col]+VariantAnnotation::geno(v)$RCZ[,sample_col], 
               VariantAnnotation::geno(v)$FGZ[,sample_col]+VariantAnnotation::geno(v)$RGZ[,sample_col], 
               VariantAnnotation::geno(v)$FTZ[,sample_col]+VariantAnnotation::geno(v)$RTZ[,sample_col]))
}

#' Dump allele counts from the DKFZ pipeline
#' 
#' Helper function that takes a sample column and fetches allele counts. As the DKFZ pipeline does not
#' provide counts for each base, but just the alt and reference, we will provide just those and the
#' other bases with 0s.
#' 
#' Note: If there are multiple ALT alleles this function will only take the first mentioned!
#' @param v The vcf file
#' @param sample_col The column in which the counts are. If it's the first sample mentioned in the vcf this would be sample_col 1
#' @return An array with 4 columns: Counts for A, C, G, T
#' @author sd11
#' @noRd
getAlleleCounts.DKFZ = function(v, sample_col) {
  # Fetch counts for both forward and reverse ref/alt
  counts = VariantAnnotation::geno(v)$DP4[,sample_col,]
  counts.ref = counts[,1] + counts[,2] # ref forward/reverse
  counts.alt = counts[,3] + counts[,4] # alt forward/reverse
  allele.ref = as.character(VariantAnnotation::ref(v))
  # allele.alt = unlist(lapply(VariantAnnotation::alt(v), function(x) { as.character(x[[1]]) }))
  allele.alt = getAltAllele(v)
  
  output = construct_allelecounter_table(count.ref, count.alt, allele.ref, allele.alt)
  return(output)
}

#' Dump allele counts from the MuSE pipeline
#' 
#' Helper function that takes a sample column and fetches allele counts. As the MuSE pipeline does not
#' provide counts for each base, but just the alt and reference, we will provide just those and the
#' other bases with 0s.
#' 
#' Note: If there are multiple ALT alleles this function will only take the first mentioned! 
#' Note2: This function assumes there are only two columns with read counts. The format allows for more
#' @param v The vcf file
#' @param sample_col The column in which the counts are. If it's the first sample mentioned in the vcf this would be sample_col 1
#' @return An array with 4 columns: Counts for A, C, G, T
#' @author sd11
#' @noRd
getAlleleCounts.MuSE = function(v, sample_col) {
  if (length(colnames(VariantAnnotation::geno(v)$AD)) > 2) {
    print("In getAlleleCounts.MuSE: Assuming 2 columns with read counts, but found more. This is not supported")
    q(save="no", status=1)
  }
  
  # An SNV can be represented by more than 1 alt alleles, here we pick the alt allele with the highest read count
  num.snvs = nrow(VariantAnnotation::geno(v)$AD)
  counts = array(NA, c(num.snvs, 2))
  allele.alt = array(NA, num.snvs)
  for (i in 1:num.snvs) {
    snv.counts = unlist(VariantAnnotation::geno(v)$AD[i,sample_col])
    counts[i,1] = snv.counts[1] # The reference is the first base for which read counts are mentioned
    select_base = which.max(snv.counts[2:length(snv.counts)])
    allele.alt[i] = as.character(VariantAnnotation::alt(v)[[i]][select_base])
    select_base = select_base+1 # The reference is the first base for which read counts are mentioned
    counts[i,2] = snv.counts[select_base] 
  }
  
  allele.ref = as.character(VariantAnnotation::ref(v))
  
  output = construct_allelecounter_table(counts[i,1], counts[i,2], allele.ref, allele.alt)
  return(output)
}

#' Dump allele counts from the Broad pipeline
#' 
#' Helper function that takes a sample column and fetches allele counts. As the Broad pipeline does not
#' provide counts for each base, but just the alt and reference, we will provide just those and the
#' other bases with 0s.
#' 
#' Note: If there are multiple ALT alleles this function will only take the first mentioned! 
#' @param v The vcf file
#' @param sample_col The column in which the counts are. If it's the first sample mentioned in the vcf this would be sample_col 1
#' @return An array with 4 columns: Counts for A, C, G, T
#' @author sd11
#' @noRd
getAlleleCounts.Broad = function(v, sample_col) {
  count.ref = as.numeric(unlist(VariantAnnotation::geno(v)$ref_count))
  count.alt = as.numeric(unlist(VariantAnnotation::geno(v)$alt_count))
  allele.ref = as.character(VariantAnnotation::ref(v))
  # allele.alt = unlist(lapply(VariantAnnotation::alt(v), function(x) { as.character(x[[1]]) }))
  allele.alt = getAltAllele(v)
  
  output = construct_allelecounter_table(count.ref, count.alt, allele.ref, allele.alt)
  return(output)
}

#' Dump allele counts in ICGC consensus SNV format
#' 
#' This function fetches allele counts from the info field in the VCF file.
#' Note: If there are multiple ALT alleles this function will only take the first mentioned! 
#' @param v The vcf file
#' @return An array with 4 columns: Counts for A, C, G, T
#' @author sd11
#' @noRd
getAlleleCounts.ICGC_consensus_snv = function(v) {
  count.alt = info(v)$t_alt_count
  count.ref = info(v)$t_ref_count
  allele.ref = as.character(VariantAnnotation::ref(v))
  # allele.alt = unlist(lapply(VariantAnnotation::alt(v), function(x) { as.character(x[[1]]) }))
  allele.alt = getAltAllele(v)
  
  output = construct_allelecounter_table(count.ref, count.alt, allele.ref, allele.alt)
  return(output)
}

#' Dump allele counts in ICGC consensus indel format
#' 
#' This function fetches allele counts from the info field in the VCF file.
#' Note: If there are multiple ALT alleles this function will only take the first mentioned! 
#' @param v The vcf file
#' @param dummy_alt_allele Dummy base to use to encode the alt allele, the counts will be saved in the corresponding column
#' @param dummy_ref_allele Dummy base to use to encode the ref allele, the counts will be saved in the corresponding column
#' @return An array with 4 columns: Counts for A, C, G, T
#' @author sd11
#' @noRd
getAlleleCounts.ICGC_consensus_indel = function(v, dummy_alt_allele="A", dummy_ref_allele="C") {
  c = construct_allelecounter_table(count.ref=info(v)$t_ref_count, 
                                    count.alt=info(v)$t_alt_count, 
                                    allele.ref=rep(dummy_ref_allele, nrow(v)),
                                    allele.alt=rep(dummy_alt_allele, nrow(v)))
  return(c)
}

#' Dump allele counts in Strelka indel format
#'
#' This function fetches allele counts from the info field in the VCF file.
#'
#' @param v The vcf file
#' @param dummy_alt_allele Dummy base to use to encode the alt allele, the counts will be saved in the corresponding column
#' @param dummy_ref_allele Dummy base to use to encode the ref allele, the counts will be saved in the corresponding column
#' @return An array with 4 columns: Counts for A, C, G, T
#' @author sd11
#' @noRd
getAlleleCounts.Strelka_indel = function(v, dummy_alt_allele="A", dummy_ref_allele="C") {
  c = construct_allelecounter_table(count.ref=geno(v)$TAR[,2,1],
				    count.alt=geno(v)$TIR[,2,1],
				    allele.ref=rep(dummy_ref_allele, nrow(v)),
				    allele.alt=rep(dummy_alt_allele, nrow(v)))
  return(c)
}

#' Dump allele counts in cgpPindel indel format
#'
#' This function fetches allele counts from the info field in the VCF file. Note that this requires running vafCorrect after cgpPindel.
#'
#' @param v The vcf file
#' @param sample_col The column in which the counts are. If it's the first sample mentioned in the vcf this would be sample_col 1
#' @param dummy_alt_allele Dummy base to use to encode the alt allele, the counts will be saved in the corresponding column
#' @param dummy_ref_allele Dummy base to use to encode the ref allele, the counts will be saved in the corresponding column
#' @return An array with 4 columns: Counts for A, C, G, T
#' @author sd11
#' @noRd
getAlleleCounts.cgpPindel_indel = function(v, sample_col, dummy_alt_allele="A", dummy_ref_allele="C") {
  c = construct_allelecounter_table(count.ref=geno(v)$WT[,sample_col],
                                    count.alt=geno(v)$MTR[,sample_col],
                                    allele.ref=rep(dummy_ref_allele, nrow(v)),
                                    allele.alt=rep(dummy_alt_allele, nrow(v)))
  return(c)
}

#' Dump allele counts in Mutect format
#' 
#' This function fetches allele counts from the info field in the VCF file.
#' Note: If there are multiple ALT alleles this function will only take the first mentioned! 
#' @param v The vcf file
#' @param sample_col The column in which the counts are. If it's the first sample mentioned in the vcf this would be sample_col 1
#' @return An array with 4 columns: Counts for A, C, G, T
#' @author sd11
#' @noRd
getAlleleCounts.mutect = function(v, sample_col) {
  counts = do.call(rbind, geno(v)$AD[,sample_col])
  count.ref = counts[,1]
  count.alt = counts[,2]
  allele.ref = as.character(VariantAnnotation::ref(v))
  # allele.alt = unlist(lapply(VariantAnnotation::alt(v), function(x) { as.character(x[[1]]) }))
  allele.alt = getAltAllele(v)
  output = construct_allelecounter_table(count.ref, count.alt, allele.ref, allele.alt)
  return(output)
}

#' Function that constructs a table in the format of the allele counter
#' @param count.ref Number of reads supporting the reference allele
#' @param count.alt Number of reads supporting the variant allele
#' @param allele.ref The reference allele
#' @param allele.alt The variant allele
#' @return A data.frame consisting of four columns: Reads reporting A, C, G and T
#' @author sd11
#' @noRd
construct_allelecounter_table = function(count.ref, count.alt, allele.ref, allele.alt) {
  output = array(0, c(length(allele.ref), 4))
  nucleotides = c("A", "C", "G", "T")
  # Propagate the alt allele counts
  nucleo.index = match(allele.alt, nucleotides)
  for (i in 1:nrow(output)) {
    output[i,nucleo.index[i]] = count.alt[i]
  }
  
  # Propagate the reference allele counts
  nucleo.index = match(allele.ref, nucleotides)
  for (i in 1:nrow(output)) {
    output[i,nucleo.index[i]] = count.ref[i]
  }
  return(output)
}
