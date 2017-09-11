# New phasing script to replace the old that used a separate perl script - sd11@sanger.ac.uk
# now requires Rsamtools v >= 2.14 

args = commandArgs(T)
samplename = args[1]
bam_file = args[2]
vcf_file = args[3]
outdir = args[4]
if (length(args) > 4) {
  fai_file = args[5]
  ign_file = args[6]
} else {
  fai_file = "/lustre/scratch116/casm/cgp/pancancer/reference/genome.fa.fai"
  ign_file = "/lustre/scratch116/casm/team113/sd11/reference/GenomeFiles/battenberg_ignore/ignore.txt"
}

max_distance = 700
minMapQ = 20

library(Rsamtools)
library(VariantAnnotation)

############################################
# VCF 2 LOCI
############################################
#' Parse genome index
#' 
#' Convenience function that parses a reference genome index as generated
#' by samtools index
#' @param fai_file The index
#' @return A data frame with columns "chromosome", "length", "offset", "fasta_line_length", "line_blen"
#' @author sd11
#' @export
parseFai = function(fai_file) {
  fai = read.table(fai_file, header=F, stringsAsFactors=F)
  colnames(fai) = c("chromosome", "length", "offset", "fasta_line_length", "line_blen")
  return(fai)
}

#' Parse chromosomes to ignore file
#' 
#' Convenience function that parses an ignore file. This file
#' is expected to have a single column with just chromosome names
#' @param ignore_file The file specifying to be ignored chromosomes
#' @return A data frame with a single column named "chromosome"
#' @author sd11
#' @export
parseIgnore = function(ignore_file) {
  ign = read.table(ignore_file, header=F, stringsAsFactors=F)
  colnames(ign) = c("chromosome")
  return(ign)
}

#' Transform vcf to loci file
#' 
#' Function that dumps the loci of snvs from a series of vcf files into a single loci file
#' @param vcf_files A vector of vcf files to be considered
#' @param fai_file Reference genome index
#' @param ign_file A file with chromosomes to be excluded from consideration
#' @param outfile Where to store the output
#' @param dummy_alt_allele The alt allele to store, supply when you want to override what is in the VCF (Default: NA)
#' @param dummy_ref_allele The reference allele to store, supply when you want to override what is in the VCF (Default: NA)
#' @author sd11
#' @export
vcf2loci = function(vcf_files, fai_file, ign_file, outfile, dummy_alt_allele=NA, dummy_ref_allele=NA) {
  fai = parseFai(fai_file)
  ign = parseIgnore(ign_file)
  allowed_chroms = which(!(fai$chromosome %in% ign$chromosome))
  
  # Run through each supplied vcf file, collect the loci from each
  combined.loci = data.frame()
  for (vcf_file in vcf_files) {
    vcf.cols = ncol(read.delim(vcf_file, comment.char="#", header=F, stringsAsFactor=F, nrows=1))
    vcf.cols.default = 10 # vcf file standard contains 10 columns
    vcf.colClasses = c(NA, NA, "NULL", NA, NA, rep("NULL", 5+(vcf.cols-vcf.cols.default)))
    vcf.loci = read.delim(vcf_file, comment.char="#", header=F, stringsAsFactor=F, colClasses=vcf.colClasses)
    colnames(vcf.loci) = c("chromosome", "pos", "ref","alt")
    vcf.loci.sel = subset(vcf.loci, chromosome %in% fai$chromosome[allowed_chroms])
    combined.loci = rbind(combined.loci, vcf.loci.sel)
  }
  
  # Override the ref and alt alleles if desired
  if (!is.na(dummy_alt_allele) & dummy_alt_allele!="NA") {
    combined.loci$alt = dummy_alt_allele
  }
  if (!is.na(dummy_ref_allele) & dummy_ref_allele!="NA") {
    combined.loci$ref = dummy_ref_allele
  }
  
  # Remove duplicate entries
  combined.loci = unique(combined.loci)
  write.table(combined.loci, col.names=F, quote=F, row.names=F, file=outfile, sep="\t")
}

############################################
# phasing
############################################

mut_mut_phasing = function(loci_file, phased_file, bam_file, bai_file, max_distance, minMapQ=20) {
  # TODO: this must be removed
  #chr.names = c(1:22,"X","Y")
  print("START")
  muts <- read.delim(loci_file, header=F, row.names=NULL, stringsAsFactors=F)
  names(muts) = c("CHR","POSITION","WT","MT")
  chr.names = unique(muts$CHR)
  muts = muts[order(match(muts$CHR,chr.names),muts$POSITION),]
  print("AFTER READING")
  # Pairwise comparison of all muts, only take those that are close to eachother
  output <- data.frame(Chr = vector(mode="character",length=0), Pos1 = vector(mode="numeric",length=0), Ref1 = vector(mode="character",length=0), Var1 = vector(mode="character",length=0), Pos2 = vector(mode="numeric",length=0), Ref2 = vector(mode="character",length=0), Var2 = vector(mode="character",length=0))
  for(chr in chr.names){
    chr.muts = muts[muts$CHR==chr,]
    no.muts = nrow(chr.muts)
    for(i in 1:(no.muts-1)){
      dist = chr.muts$POSITION[(i+1):no.muts] - chr.muts$POSITION[i]
      inds = which(dist <= max_distance) # 700
      if(length(inds)>0){
        output <- rbind(output,data.frame(Chr = chr, Pos1 = chr.muts$POSITION[i], Ref1 = chr.muts$WT[i], Var1 = chr.muts$MT[i], Pos2 = chr.muts$POSITION[i+inds], Ref2 = chr.muts$WT[i+inds], Var2 = chr.muts$MT[i+inds]))
      }
    }
  }
  print("AFTER FIRST FOR")
  
  if (nrow(output) > 0) {
    print("IN IF")
    count.data = data.frame()
    
    for (i in 1:nrow(output)) {
      save(file="test.RData", output, i)
      chrom = output$Chr[i]
      dat_pair = data.frame(chrom=rep(output$Chr[i],2), start=c(output$Pos1[i], output$Pos2[i]), end=c(output$Pos1[i], output$Pos2[i]))
      dat_pair = makeGRangesFromDataFrame(dat_pair)
      
      wt_wt = paste(output$Ref1[i], output$Ref2[i], sep="_")
      mut_wt = paste(output$Var1[i], output$Ref2[i], sep="_")
      wt_mut = paste(output$Ref1[i], output$Var2[i], sep="_")
      mut_mut = paste(output$Var1[i], output$Var2[i], sep="_")
    
      flag = scanBamFlag(isPaired=T, hasUnmappedMate=F, isDuplicate=F, isUnmappedQuery=F) #, isProperPair=T
      param <- ScanBamParam(which=dat_pair[1,], what=scanBamWhat(), flag=flag, mapqFilter=minMapQ) #, mapqFilter=minMapQ
      bamfile <- BamFile(bam_file, asMates=TRUE)
      bam <- scanBam(bamfile, param=param)[[1]]
    
      alleles = get_allele_combination_counts(bam, dat_pair)
      
      #' SNV to SNV
      allele_pairs = table(factor(paste(alleles$snv1, alleles$snv2, sep="_"), levels=c(wt_wt, mut_mut, wt_mut, mut_wt)))
      count.data = rbind(count.data, data.frame(output[i,], 
                                                Num_WT_WT=allele_pairs[[wt_wt]], 
                                                Num_Mut_Mut=allele_pairs[[mut_mut]], 
                                                Num_Mut_WT=allele_pairs[[mut_wt]], 
                                                Num_WT_Mut=allele_pairs[[wt_mut]]))
    }
    print("AFTER FIRST FOR in IF")
    # Categorise pairs of mutations
    count.data$phasing = NA
    for(h in 1:nrow(count.data)){
      counts = count.data[h,8:11]
      print(counts) 
      if(counts[2]>0 & counts[3]+counts[4] == 0){
        count.data$phasing[h]="phased"
        #}else if(counts[2]==0 & counts[3]+counts[4] > 0){
        #require BOTH WT-mut and mut-WT pairs
      }else if(counts[2]==0 & counts[3]>0 & counts[4]>0){  
        count.data$phasing[h]="anti-phased"
      }else if(counts[2]>0 && counts[3]>0 && counts[4]==0){
        count.data$phasing[h]="clone-subclone"
      }else if(counts[2]>0 && counts[3]==0 && counts[4]>0){
        count.data$phasing[h]="subclone-clone"
      }  
    }
    print("AFTER second FOR in IF")
    write.table(count.data[!is.na(count.data$phasing),],phased_file,sep="\t",quote=F,row.names=F)
  }
}


mut_cn_phasing = function(loci_file, phased_file, hap_file, bam_file, bai_file, outfile, max_distance) {
  
  if (file.info(loci_file)$size == 0) {
    linked.muts = data.frame(matrix(rep(NA, 13), nrow=1))
    colnames(linked.muts) = c("Chr","Pos1","Ref1","Var1","Pos2","Ref2","Var2","AF","AFphased","Num_linked_to_A","Num_linked_to_C","Num_linked_to_G","Num_linked_to_T")
    linked.muts = na.omit(linked.muts)
  } else if(nrow(read.delim(loci_file, header=F, stringsAsFactors=F, fill=T))==0) {
    linked.muts = data.frame(matrix(rep(NA, 13), nrow=1))
    colnames(linked.muts) = c("Chr","Pos1","Ref1","Var1","Pos2","Ref2","Var2","AF","AFphased","Num_linked_to_A","Num_linked_to_C","Num_linked_to_G","Num_linked_to_T")
    linked.muts = na.omit(linked.muts)
  } else {
    # TODO: is this filtering required when just supplying loci files from a single chromosome?
    chr.muts = read.delim(loci_file, header=F, stringsAsFactors=F, fill=T)
    names(chr.muts) = c("CHR","POSITION","REF_BASE","MUT_BASE")
    
    # Match phased SNPs and their haplotypes together
    phased = read.delim(phased_file, header=T, stringsAsFactors=F, quote="\"")
    # Compatible with both BB setups that have row.names and those that don't
    if (ncol(phased) == 6) {
      colnames(phased) = c("SNP", "Chr","Pos", "AF", "AFphased", "AFsegmented")
    } else {
      colnames(phased) = c("Chr","Pos", "AF", "AFphased", "AFsegmented")
    }
    
    # TODO: check that chromosomes are using the same names between loci and phased files
    phased = phased[phased$Chr %in% chr.muts$CHR,]
    
    hap.info = read.delim(hap_file, sep=" ", header=F, row.names=NULL, stringsAsFactors=F)
    # Compatible with both BB setups that have row.names and those that don't
    if (ncol(hap.info) == 7) {
      colnames(hap.info) = c("SNP","dbSNP","pos","ref","mut","ref_count","mut_count")
    } else {
      colnames(hap.info) = c("dbSNP","pos","ref","mut","ref_count","mut_count")
    }
    # get haplotypes that match phased heterozygous SNPs
    hap.info = hap.info[match(phased$Pos,hap.info$pos),]
    
    # Synchronise dfs in case some SNPs are not in hap.info
    selection = !is.na(hap.info$pos)
    hap.info = hap.info[selection,]
    phased = phased[selection,]
    
    #220212
    #phased$AF[hap.info$ref_count==1] = 1-phased$AF[hap.info$ref_count==1]
    phased$Ref = hap.info$ref
    phased$Var = hap.info$mut
    
    # Annotate the chr.muts df with the min abs distance to a phased SNP
    chr.muts$dist = sapply(1:dim(chr.muts)[1], function(i,p,m) min(abs(p$Pos - m$POSITION[i])), p=phased, m=chr.muts)
    chr.muts$snp.index = sapply(1:dim(chr.muts)[1], function(i,p,m) which.min(abs(p$Pos - m$POSITION[i])), p=phased, m=chr.muts)
    # Use only those to a SNP
    muts <- chr.muts[chr.muts$dist < max_distance,] #700
    snps <- phased[muts$snp.index,]
    
    # linked.muts = run_linkage_pull_snp(loci_file, bam_file, bai_file, muts$CHR, muts$POSITION, muts$REF_BASE, muts$MUT_BASE, snps$Pos, snps$Ref, snps$Var, snps$AF, snps$AFphased)
    
    output = data.frame(Chr=muts$CHR, 
                        Pos1=muts$POSITION, 
                        Ref1=muts$REF_BASE, 
                        Var1=muts$MUT_BASE, 
                        Pos2=snps$Pos, 
                        Ref2=snps$Ref, 
                        Var2=snps$Var, 
                        AF=snps$AF, 
                        AFphased=snps$AFphased)
    
    linked.muts = data.frame()
    for (i in 1:nrow(output)) {
    
      chrom = output$Chr[i]
      dat_pair = data.frame(chrom=rep(output$Chr[i],2), start=c(output$Pos1[i], output$Pos2[i]), end=c(output$Pos1[i], output$Pos2[i]))
      dat_pair = makeGRangesFromDataFrame(dat_pair)
      
      flag = scanBamFlag(isPaired=T, hasUnmappedMate=F, isDuplicate=F, isUnmappedQuery=F) #, isProperPair=T
      param <- ScanBamParam(which=dat_pair[1,], what=scanBamWhat(), flag=flag)
      bamfile <- BamFile(bam_file, asMates=TRUE)
      bam <- scanBam(bamfile, param=param)[[1]]
      
      count.data = get_allele_combination_counts(bam, dat_pair)
      #' SNV to SNP
      alleles_snv_ref = count.data[count.data$snv1==as.character(output$Ref1[i]),]
      nucl_counts = sapply(c("A", "C", "G", "T"), function(nucl) { sum(alleles_snv_ref$snv2==nucl, na.rm=T) })
      linked.muts = rbind(linked.muts, data.frame(output[i,],
                                                  Num_linked_to_A=nucl_counts[["A"]],
                                                  Num_linked_to_C=nucl_counts[["C"]],
                                                  Num_linked_to_G=nucl_counts[["G"]],
                                                  Num_linked_to_T=nucl_counts[["T"]]))
    }      
    
    # Categorise where the mutation is with respect to the CN event
    ACGT = 10:13
    names(ACGT) <- c("A", "C", "G", "T")
    linked.muts$Parental <- rep(NA, dim(linked.muts)[1])
    if (nrow(linked.muts) > 0) {
      for (i in 1:nrow(linked.muts)) {
        
        # Fetch allele frequency
        af = linked.muts$AF[i]
        # Get number of reads covering the ref mutation allele
        ref_count = hap.info[hap.info$pos==linked.muts$Pos2[i],]$ref_count
        # Get number of reads covering the alt mutation allele
        alt_count = hap.info[hap.info$pos==linked.muts$Pos2[i],]$mut_count
        # Number of reads covering SNP allele A, that also cover mutation alt
        linked_to_A = linked.muts[i,ACGT[linked.muts$Ref2[i]]]
        # Number of reads covering SNP allele B, that also cover mutation alt
        linked_to_B = linked.muts[i,ACGT[linked.muts$Var2[i]]]
        
        if (af < 0.5 & alt_count==1 & linked_to_A > 0 & linked_to_B == 0) {
          linked.muts$Parental[i] = "MUT_ON_DELETED"
        } else if (af < 0.5 & alt_count==1 & linked_to_A == 0 & linked_to_B > 0) {
          linked.muts$Parental[i] = "MUT_ON_RETAINED"
        } else if (af > 0.5 & alt_count==1 & linked_to_A > 0 & linked_to_B == 0) {
          linked.muts$Parental[i] = "MUT_ON_RETAINED"
        } else if (af > 0.5 & alt_count==1 & linked_to_A == 0 & linked_to_B > 0) {
          linked.muts$Parental[i] = "MUT_ON_DELETED"
        } else if (af > 0.5 & ref_count==1 & linked_to_A > 0 & linked_to_B == 0) {
          linked.muts$Parental[i] = "MUT_ON_DELETED"
        } else if (af > 0.5 & ref_count==1 & linked_to_A == 0 & linked_to_B > 0) {
          linked.muts$Parental[i] = "MUT_ON_RETAINED"
        } else if (af < 0.5 & ref_count==1 & linked_to_A > 0 & linked_to_B == 0) {
          linked.muts$Parental[i] = "MUT_ON_RETAINED"
        } else if (af < 0.5 & ref_count==1 & linked_to_A == 0 & linked_to_B > 0) {
          linked.muts$Parental[i] = "MUT_ON_DELETED"
        }
      }
    }
  }
  
  write.table(linked.muts,outfile, sep="\t", quote=F, row.names=F)
}



get_allele_combination_counts = function(bam, dat_pair) {
  alleles = data.frame()
  
  # alleles = lapply(unique(bam$groupid), function(i) {
  #for (j in seq(1, (length(bam$seq)), 2)) {
  for (i in unique(bam$groupid)) {
    if (sum(bam$groupid==i) != 2) {
      # No two reads in this pair - possibly removed due to mapping quality constraints
      next
    }
    
    # Take index of first read in pair
    j = which(bam$groupid==i)[1]
    
    rel_pos_snv1 = NA
    read_snv1 = NA
    rel_pos_snv2 = NA
    read_snv2 = NA
    read_1 = NA
    read_2 = NA
    
    rel_pos_read1_snv1 = (start(dat_pair[1])-bam$pos[j]) + 1
    rel_pos_read2_snv1 = (start(dat_pair[1])-bam$mpos[j]) + 1
    if (rel_pos_read1_snv1 < bam$qwidth[j] & rel_pos_read1_snv1 > 0) {
      rel_pos_snv1 = rel_pos_read1_snv1
      read_snv1 = j
      read_1 = "first"
    } else if (rel_pos_read2_snv1 < bam$qwidth[j] & rel_pos_read2_snv1 > 0) {
      rel_pos_snv1 = rel_pos_read2_snv1
      read_snv1 = j+1
      read_1 = "second"
    }
    
    rel_pos_read1_snv2 = (start(dat_pair[2])-bam$pos[j]) + 1
    rel_pos_read2_snv2 = (start(dat_pair[2])-bam$mpos[j]) + 1
    if (rel_pos_read1_snv2 < bam$qwidth[j] & rel_pos_read1_snv2 > 0) {
      rel_pos_snv2 = rel_pos_read1_snv2
      read_snv2 = j
      read_2 = "first"
    } else if (rel_pos_read2_snv2 < bam$qwidth[j] & rel_pos_read2_snv2 > 0) {
      rel_pos_snv2 = rel_pos_read2_snv2
      read_snv2 = j+1
      read_2 = "second"
    }
    
    if (!is.na(rel_pos_snv1)) {
      base_snv1 = substr(as.character(bam$seq[read_snv1]), rel_pos_snv1, rel_pos_snv1)
    } else {  
      base_snv1 = NA
    }
    
    if (!is.na(rel_pos_snv2)) {
      base_snv2 = substr(as.character(bam$seq[read_snv2]), rel_pos_snv2, rel_pos_snv2)
    } else {
      base_snv2 = NA
    }
    
    alleles = rbind(alleles, data.frame(snv1=base_snv1, snv2=base_snv2, read_1=read_1, read_2=read_2))
  }#)
  
  # alleles = do.call(rbind, alleles)
  alleles_m = na.omit(alleles)
  return(alleles_m)
}

############################################
# script
############################################
loci_file = file.path(outdir, paste0(samplename, "_loci.txt"))
phased_file = file.path(outdir, paste0(samplename, "_mutmut_phasing.txt"))
vcf2loci(list(vcf_file), fai_file, ign_file, loci_file)

mut_mut_phasing(loci_file=loci_file,
                phased_file=phased_file,
                bam_file=bam_file,
                bai_file=paste0(bam_file, ".bai"),
                max_distance=max_distance,
                minMapQ=minMapQ)



