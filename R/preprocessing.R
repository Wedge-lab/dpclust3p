ALLELECOUNTER = "alleleCounter"
LINKAGEPULL = "Linkage_pull.pl"

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
    read_data = tryCatch(read.delim(vcf_file, comment.char="#", header=F, stringsAsFactor=F, nrows=1),
                         error = function(e) NULL)
    
    # check that there was data in the file
    if (!is.null(read_data)) {
      vcf.cols = ncol(read_data)
      vcf.cols.default = 10 # vcf file standard contains 10 columns
      vcf.colClasses = c(NA, NA, "NULL", NA, NA, rep("NULL", 5+(vcf.cols-vcf.cols.default)))
      vcf.loci = read.delim(vcf_file, comment.char="#", header=F, stringsAsFactor=F, colClasses=vcf.colClasses)
      colnames(vcf.loci) = c("chromosome", "pos", "ref","alt")
      vcf.loci.sel = subset(vcf.loci, chromosome %in% fai$chromosome[allowed_chroms])
      combined.loci = rbind(combined.loci, vcf.loci.sel)
    }
  }
  
  if (nrow(combined.loci) > 0) {
    # Override the ref and alt alleles if desired
    if (!is.na(dummy_alt_allele) & dummy_alt_allele!="NA") {
      combined.loci$alt = dummy_alt_allele
    }
    if (!is.na(dummy_ref_allele) & dummy_ref_allele!="NA") {
      combined.loci$ref = dummy_ref_allele
    }
  } else {
    # no data found, generate an empty data.frame
    combined.loci = data.frame(matrix(ncol = 4, nrow = 0))
    colnames(combined.loci) = c("chromosome", "pos", "ref","alt")
  }
  
  # Remove duplicate entries
  chrpos = paste0(combined.loci$chromosome, "_", combined.loci$pos)
  chrpos_dup = chrpos[duplicated(chrpos)]
  combined.loci = combined.loci[!chrpos %in% chrpos_dup,]
  write.table(combined.loci, col.names=F, quote=F, row.names=F, file=outfile, sep="\t")
}

############################################
# Mutation signatures
############################################

#' Obtain tri-nucleotide context. It reads in the supplied loci file and querries
#' the provided reference genome fasta for the context. Finally it writes out the
#' loci with annotated context.
#' @param loci_file A 4 column loci file with chrom, position, reference allele and alt allele
#' @param outfile Where to write the output
#' @param ref_genome Full path to a reference genome Fasta file
#' @author sd11
#' @export
getTrinucleotideContext = function(loci_file, outfile, ref_genome) {
  loci = read.table(loci_file, sep='\t', header=F, stringsAsFactors=F)
  loci.g = GenomicRanges::GRanges(seqnames=loci[,1], ranges=IRanges::IRanges(loci[,2], loci[,2]))
  r = Rsamtools::scanFa(file=ref_genome, resize(GenomicRanges::granges(loci.g), 3, fix="center"))
  write.table(cbind(loci, as.data.frame(r)), file=outfile, row.names=F, col.names=F, quote=F, sep="\t")
}

# #' Checks each tri-nucleotide context whether its CAG or CTG
# #' @param vector_of_trinucleotides A vector containing the base to check and its immediate neighbors
# #' @author sd11
# isDeaminase = function(vector_of_trinucleotides) {
#   return(grepl("(CAG)|(CTG)|(GAC)|(GTC)", vector_of_trinucleotides))
# }

#' Filter a with tri-nucleotide context annotated list of loci for a particular signature
#' supplied as a regex like so: (CAG)|(CTG)|(GAC)|(GTC).
#' @param signature_anno_loci_file Filepath to a with tri-nucleotide context annotated loci
#' @param signature_regex A regex that captures the signature to keep.
#' @param outfile Filepath to where to store the output.
#' @param trinucleotide_column Integer representing the column within the input file that contains the context
#' @param alt_alleles An optional vector with reference alleles to allow.
#' @param alt_allele_column The column in the annotated loci file that contains the reference base.
#' @author sd11
#' @export
filterForSignature = function(signature_anno_loci_file, signature_regex, outfile, trinucleotide_column=5, alt_alleles=NA, alt_allele_column=4) {
  signature_anno_loci = read.table(signature_anno_loci_file, sep='\t', header=F, stringsAsFactors=F)
  signature_anno_loci_filt = signature_anno_loci[grepl(signature_regex, signature_anno_loci[,trinucleotide_column]), ]

  if (!is.na(alt_alleles)) {
    regex_split = gsub("(", "", signature_regex, fixed=T)
    regex_split = gsub(")", "", regex_split, fixed=T)
    regex_split = unlist(strsplit(regex_split, "|", fixed=T))
    selection = rep(F, nrow(signature_anno_loci_filt))
    
    # Go through all signatures and their specific reference context
    for (i in 1:length(regex_split)) {
      alt_i = alt_alleles[i]
      sign_i = regex_split[i]
      # Select only those mutations that have this particular context and reference
      selection[signature_anno_loci_filt[,trinucleotide_column]==sign_i & signature_anno_loci_filt[,alt_allele_column]==alt_i] = T
    }
    signature_anno_loci_filt = signature_anno_loci_filt[selection,]
  }
  write.table(signature_anno_loci_filt, file=outfile, row.names=F, col.names=F, quote=F, sep="\t")
}

#' Convenience function that filters a with tri-nucleotide context annotated list of loci for 
#' cytosine deaminase signature, or C>T at CpG.
#' @param loci_file Filepath to a with tri-nucleotide context annotated loci
#' @param outfile Filepath to where to store the output.
#' @param ref_genome Full path to an indexed reference genome fasta file
#' @param trinucleotide_column Integer representing the column within the input file that contains the context
#' @param alt_allele_column The column in the annotated loci file that contains the reference base.
#' @author sd11
#' @export
filterForDeaminase = function(loci_file, outfile, ref_genome, trinucleotide_column=5, alt_allele_column=4) {
  signature_anno_loci_file = gsub(".txt", "_signature_anno.txt", loci_file)
  getTrinucleotideContext(loci_file, signature_anno_loci_file, ref_genome)
  filterForSignature(signature_anno_loci_file=signature_anno_loci_file, 
                     signature_regex="(CCG)|(GCC)|(CGG)|(GGC)", 
                     outfile=outfile, 
                     trinucleotide_column=trinucleotide_column, 
                     alt_alleles=c("T", "T", "A", "A"), 
                     alt_allele_column=alt_allele_column)
}

############################################
# MUT 2 MUT phasing
############################################
#' Run linkage pull on SNVs vs SNVs
#' @noRd
run_linkage_pull_mut = function(output, loci_file, bam_file, bai_file) {
  #
  # Runs Linkage_pull 
  #
  nearbyMuts_file = paste(loci_file, "_nearbyMuts.txt", sep="")
  write.table(output, nearbyMuts_file, sep="\t", quote=F, col.names=F, row.names=F)
  
  zoomed_phase_file = paste(loci_file, "_zoomed_phase_info.txt", sep="")
  cmd=paste(LINKAGEPULL," ",nearbyMuts_file," ",bam_file," ",bai_file," ",zoomed_phase_file," Mut",sep="")
  print(paste("command=",cmd,sep=""))
  system(cmd,wait=T)
  
  count.data = read.delim(zoomed_phase_file,sep="\t",header=T, stringsAsFactors=F)
  
  return(count.data)
}

#' Phase mutation to mutation
#' 
#' Run mutation to mutation phasing. This function requires the Linkage_pull.pl script in $PATH.
#' @param loci_file A list of loci
#' @param phased_file File to save the output
#' @param bam_file Full path to the bam file
#' @param bai_file Full path to the bai file
#' @param max_distance The max distance of a mutation and SNP can be apart to be considered for phasing
#' @author sd11, dw9
#' @export
mut_mut_phasing = function(loci_file, phased_file, bam_file, bai_file, max_distance) {
  # Check if there are lines in the file, otherwise it will crash this script
  if (file.info(loci_file)$size == 0) {
    print("No lines in loci file")
    q(save="no")
  }
  
  # TODO: this must be removed
  chr.names = c(1:22,"X","Y")
  
  muts <- read.delim(loci_file, header=F, row.names=NULL, stringsAsFactors=F)
  names(muts) = c("CHR","POSITION","WT","MT")
  muts = muts[order(match(muts$CHR,chr.names),muts$POSITION),]
  
  # # TODO: Does this work with X and Y ??
  # if (!is.null(chrom)) {
  #   muts = muts[muts$CHR==chrom,]
  # }
  
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
  
  if (nrow(output) > 0) {
    # Run linkage pull on the chromosome locations mentioned in the data.frame output
    count.data = run_linkage_pull_mut(output, loci_file, bam_file, bai_file)
    
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
    
    write.table(count.data[!is.na(count.data$phasing),],phased_file,sep="\t",quote=F,row.names=F)
  }
}

############################################
# MUT 2 CN phasing
############################################
#' Run linkage pull on SNVs vs SNPs
#' @noRd
run_linkage_pull_snp = function(loci_file, bam_file, bai_file, chr, pos1, ref1, var1, pos2, ref2, var2, af, af_phased) {
  #
  # Runs the Linkage_pull.pl script with the given columns as its input. 
  # Returns a dataframe with the given columns, plus the linkage information that was pulled from the BAM file
  #
  
  output = data.frame(Chr=chr, Pos1=pos1, Ref1=ref1, Var1=var1, Pos2=pos2, Ref2=ref2, Var2=var2, AF=af, AFphased=af_phased)
  
  linkedFile = paste(loci_file, "_muts_linkedSNPs.txt",sep="")
  write.table(output[,1:7], linkedFile, sep="\t", quote=F, col.names=F, row.names=F)
  
  outfile = paste(loci_file, "_zoomed_mutcn_phase.txt", sep="")
  cmd=paste(LINKAGEPULL," ",linkedFile," ",bam_file," ",bai_file," ",outfile," SNP",sep="")
  print(paste("command=",cmd,sep=""))
  system(cmd,wait=T)
  
  input <- read.delim(outfile, header=T, sep="\t", stringsAsFactors=F)
  linked.muts <- cbind(output, input[,8:11])
  
  return(linked.muts)
}

#' Phase mutation to SNP/copy number
#' 
#' Run mutation to copy number phasing. This function requires the Linkage_pull.pl script in $PATH.
#' Note: This function should either be run separately per chromosome and then combined with \code{\link{concat_files}}
#' or on all chromsomes in one go, but then the _allHaplotypeInfo.txt Battenberg files need to be concatenated first.
#' @param loci_file A list of loci
#' @param phased_file The .BAFsegmented.txt output file from Battenberg
#' @param hap_file Path to the _allHaplotypeInfo.txt Battenberg output file to be used
#' @param bam_file Full path to the bam file
#' @param bai_file Full path to the bai file
#' @param outfile File to save the output
#' @param max_distance The max distance of a mutation and SNP can be apart to be considered for phasing
#' @author sd11, dw9
#' @export
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
    
    linked.muts = run_linkage_pull_snp(loci_file, bam_file, bai_file, muts$CHR, muts$POSITION, muts$REF_BASE, muts$MUT_BASE, snps$Pos, snps$Ref, snps$Var, snps$AF, snps$AFphased)
    
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

############################################
# Combine all the steps into a DP input file
############################################

#' Main function that creates the DP input file. A higher level function should be called by users
#' @noRd
GetDirichletProcessInfo<-function(outputfile, cellularity, info, subclone.file, is.male = F, out.dir = NULL, SNP.phase.file = NULL, mut.phase.file = NULL, adjust_male_y_chrom=F){

  write_output = function(info, outputfile) {
    if (!is.null(info)) {
  	  # convert GenomicRanges object to df
  	  df = data.frame(chr=as.data.frame(seqnames(info)),
  			  start=start(info)-1,
  			  end=end(info))
            df = cbind(df, as.data.frame(elementMetadata(info)))
            colnames(df)[1] = "chr"
            df = df[with(df, order(chr)),]
            print(head(df))
            
    } else {
      df = data.frame(matrix(ncol=16, nrow=0))
      colnames(df) = c("chr", "start", "end", "WT.count", "mut.count", "subclonal.CN", "nMaj1","nMin1", "frac1", "nMaj2", "nMin2", "frac2", "phase", "mutation.copy.number", "subclonal.fraction", "no.chrs.bearing.mut")
    }
    write.table(df, outputfile, sep="\t", row.names=F, quote=F)
  }
  
  if (is.null(info)) {
    # No entries were found in the supplied VCF file, therefore generate an empty output file
    write_output(info, outputfile)
    return()
  }
  
  subclone.data = read.table(subclone.file,sep="\t",header=T,stringsAsFactors=F)
  # Add in the Y chrom if donor is male and Battenberg hasn't supplied it (BB returns X/Y ad multiple copies of X for men)
  if (is.male & (! "Y" %in% subclone.data$chr) & adjust_male_y_chrom) {
    subclone.data = addYchromToBattenberg(subclone.data)
  }
  subclone.data.gr = GenomicRanges::GRanges(subclone.data$chr, IRanges::IRanges(subclone.data$startpos, subclone.data$endpos), rep('*', nrow(subclone.data)))
  elementMetadata(subclone.data.gr) = subclone.data[,3:ncol(subclone.data)]
  subclone.data.gr = sortSeqlevels(subclone.data.gr)
  
  info_anno = as.data.frame(cbind(array(NA, c(length(info), 7)))) 
  colnames(info_anno) = c('subclonal.CN','nMaj1','nMin1','frac1','nMaj2','nMin2','frac2')
  inds = findOverlaps(info, subclone.data.gr)  
  info_anno[queryHits(inds),2:7] = subclone.data[subjectHits(inds),][,c("nMaj1_A", "nMin1_A", "frac1_A", "nMaj2_A", "nMin2_A", "frac2_A")]
  
  CN1 = (info_anno[queryHits(inds),]$nMaj1 + info_anno[queryHits(inds),]$nMin1) * info_anno[queryHits(inds),]$frac1
  # If frac is not one for allele 1 (i.e. not only CN data for allele 1), add the CN contribution of allele 2 as well
  CN2 = (info_anno[queryHits(inds),]$nMaj2 + info_anno[queryHits(inds),]$nMin2) * info_anno[queryHits(inds),]$frac2 * ifelse(info_anno[queryHits(inds),]$frac1 != 1, 1, 0)
  CN2[is.na(CN2)] = 0
  info_anno[queryHits(inds),]$subclonal.CN = CN1 + CN2
  elementMetadata(info) = cbind(as.data.frame(elementMetadata(info)), info_anno)
  
  info$phase="unphased"
  if (!is.null(SNP.phase.file) & SNP.phase.file!="NA") {
    phasing = read.table(SNP.phase.file, header=T, stringsAsFactors=F) #header=T, skip=1, 
    phasing.gr = GenomicRanges::GRanges(phasing$Chr, IRanges::IRanges(phasing$Pos1, phasing$Pos1))
    phasing.gr$phasing = phasing$Parental
    inds = findOverlaps(info, phasing.gr)  
    info$phase[queryHits(inds)] = phasing.gr$phasing[subjectHits(inds)]
    
    info$phase[is.na(info$phase)]="unphased"
  }
  
  if(is.male & "chr" %in% names(info)){
    normal.CN = rep(2,nrow(info))
    normal.CN[info$chr=="X"| info$chr=="Y"] = 1
    info$mutation.copy.number = mutationBurdenToMutationCopyNumber(info$mut.count/ (info$mut.count + info$WT.count) , info$subclonal.CN, cellularity, normal.CN)
  }else{
    info$mutation.copy.number = mutationBurdenToMutationCopyNumber(info$mut.count/ (info$mut.count + info$WT.count) , info$subclonal.CN, cellularity)
  }
  
  # convert MCN to subclonal fraction - tricky for amplified mutations
  info$subclonal.fraction = info$mutation.copy.number
  expected.burden.for.MCN = dpclust3p:::mutationCopyNumberToMutationBurden(rep(1,length(info)),info$subclonal.CN,cellularity)
  info$no.chrs.bearing.mut = NA

  # Check if any SNVs have a copy number state, if there are none we can stop here
  if (all(is.na(info$nMaj1))) {
	  write_output(info, outputfile)
	  return()
  }

  non.zero.indices = which(info$mut.count>0 & !is.na(expected.burden.for.MCN))
  #test for mutations in more than 1 copy
  p.vals = sapply(1:length(non.zero.indices),function(v,e,i){
    prop.test(v$mut.count[i],v$mut.count[i] + v$WT.count[i],e[i],alternative="greater")$p.value 
  },v=info[non.zero.indices,], e=expected.burden.for.MCN[non.zero.indices])
  amplified.muts = non.zero.indices[p.vals<=0.05]
  
  info$no.chrs.bearing.mut = 1	
  
  #copy numbers of subclones can only differ by 1 or 0 (as assumed when calling subclones)
  if(length(amplified.muts)>0){		
    for(a in 1:length(amplified.muts)){
      max.CN2=0
      #use phasing info - if on 'deleted' (lower CN chromosome), use the minor copy number
      if(info$phase[amplified.muts[a]]=="MUT_ON_DELETED"){
        print("mut on minor chromosome")
        max.CN1 = info$nMin1[amplified.muts[a]]
        frac1 = info$frac1[amplified.muts[a]]
        frac2=0
        if(!is.na(info$nMin2[amplified.muts[a]])){
          #swap subclones, so that the one with the higher CN is first
          if(info$nMin2[amplified.muts[a]]>max.CN1){
            max.CN2 = max.CN1
            max.CN1 = info$nMin2[amplified.muts[a]]
            frac2 = frac1
            frac1 = info$frac2[amplified.muts[a]]
          }else{
            max.CN2 = info$nMin2[amplified.muts[a]]
            frac2 = info$frac2[amplified.muts[a]]
          }
        }					
      }else{
        max.CN1 = info$nMaj1[amplified.muts[a]]
        frac1 = info$frac1[amplified.muts[a]]
        frac2=0
        if(!is.na(info$nMaj2[amplified.muts[a]])){
          #swap subclones, so that the one with the higher CN is first
          if(info$nMaj2[amplified.muts[a]]>max.CN1){
            max.CN2 = max.CN1
            max.CN1 = info$nMaj2[amplified.muts[a]]
            frac2 = frac1
            frac1 = info$frac2[amplified.muts[a]]						
          }else{
            max.CN2 = info$nMaj2[amplified.muts[a]]
            frac2 = info$frac2[amplified.muts[a]]
          }
        }	
      }
      best.err = info$mutation.copy.number[amplified.muts[a]] - 1
      best.CN=1
      for(j in 1:max.CN1){
        for(k in (j-1):min(j,max.CN2)){
          potential.CN = j * frac1 + k * frac2
          err = abs(info$mutation.copy.number[amplified.muts[a]]/potential.CN-1)
          if(err<best.err){
            info$no.chrs.bearing.mut[amplified.muts[a]] = potential.CN
            best.err=err
            best.CN = potential.CN
          }
        }
      }
      info$subclonal.fraction[amplified.muts[a]] = info$mutation.copy.number[amplified.muts[a]] / best.CN
    }
  }
  
  ##########################################################################
  #test for subclonal mutations
  
  #test whether mut burden is less than expected value for MCN = 1
  p.vals1 = sapply(1:length(non.zero.indices),function(v,e,i){
    prop.test(v$mut.count[i],v$mut.count[i] + v$WT.count[i], e[i], alternative="less")$p.value
  },v=info[non.zero.indices,], e=expected.burden.for.MCN[non.zero.indices])
  #test whether mut burden is above error rate (assumed to be 1 in 200)
  p.vals2 = sapply(1:length(non.zero.indices),function(v,i){
    prop.test(v$mut.count[i],v$mut.count[i] + v$WT.count[i],0.005,alternative="greater")$p.value
  },v=info[non.zero.indices,])
  
  subclonal.muts = non.zero.indices[p.vals1<=0.05 & p.vals2<=0.05]
  
  # use subclonal CN that minimises the difference in subclonal fraction from 1
  if(length(subclonal.muts)>0){
    for(a in 1:length(subclonal.muts)){
      #if there are no subclonal CNVs, don't adjust subclonal fraction
      if(is.na(info$frac2[subclonal.muts[a]])){next}
      #assume subclonal muts are on one chromosome copy, therefore mutation copy number must be subclonal fraction of the higher CN subclone (i.e. lost in lower CN subclone) or 1 (i.e. present in both subclones)
      if(info$nMaj1[subclonal.muts[a]]+info$nMin1[subclonal.muts[a]] > info$nMaj2[subclonal.muts[a]]+info$nMin2[subclonal.muts[a]]){	
        possible.subclonal.fractions = c(info$frac1[subclonal.muts[a]],1)
      }else{
        possible.subclonal.fractions = c(info$frac2[subclonal.muts[a]],1)
      }
      best.CN = possible.subclonal.fractions[which.min(abs(info$mutation.copy.number[subclonal.muts[a]]/possible.subclonal.fractions - 1))]
      #extra test 200313 - check whether subclonal CN results in clonal mutation, otherwise subclonal CN doesn't explain subclonal MCN
      if(best.CN != 1 & prop.test(info$mut.count[subclonal.muts[a]],info$mut.count[subclonal.muts[a]]+info$WT.count[subclonal.muts[a]],expected.burden.for.MCN[subclonal.muts[a]] * best.CN)$p.value > 0.05){
        info$subclonal.fraction[subclonal.muts[a]] = info$mutation.copy.number[subclonal.muts[a]] / best.CN
        info$no.chrs.bearing.mut[subclonal.muts[a]] = best.CN
      }
    }
  }	
   
  #possible.zero.muts = intersect((1:length(info))[-non.zero.indices],which(!is.na(info$nMin1)))
  #possible.zero.muts = c(possible.zero.muts,non.zero.indices[p.vals2>0.05])
  #if(length(possible.zero.muts)>0){
  #  del.indices = which(info$nMin1[possible.zero.muts]==0 & !info$phase[possible.zero.muts]=="MUT_ON_RETAINED")
  #  info$subclonal.fraction[possible.zero.muts[del.indices]] = NA
  #  info$no.chrs.bearing.mut[possible.zero.muts[del.indices]] = 0
  #}
  
  write_output(info, outputfile)
}

#' Convenience function to load the cellularity from a rho_and_psi file
#' @noRd
GetCellularity <- function(rho_and_psi_file) {
  d = read.table(rho_and_psi_file, header=T, stringsAsFactors=F)
  return(d['FRAC_GENOME','rho'])
}

#' Convenience function to fetch WTCount and mutCount
#'@noRd
GetWTandMutCount <- function(loci_file, allele_frequencies_file) {
  subs.data = tryCatch(read.table(loci_file, sep='\t', header=F, stringsAsFactors=F), error=function(e) NA)
  if (is.na(subs.data)) {
    # Empty input
    return(NULL)
  }
  subs.data = subs.data[order(subs.data[,1], subs.data[,2]),]

  if (nrow(subs.data)!=nrow(unique(subs.data[, 1:2]))) {
	  stop("The input contains at least one pair of mutations at the same position. Chromosome and position must be unique for this pipeline to run.")
  }
  
  # Replace dinucleotides and longer with just the first base. Here we assume the depth of the second base is the same and the number of dinucleotides is so low that removing the second base is negligable
  subs.data[,3] = apply(as.data.frame(subs.data[,3]), 1, function(x) { substring(x, 1,1) })
  subs.data[,4] = apply(as.data.frame(subs.data[,4]), 1, function(x) { substring(x, 1,1) })
  
  subs.data.gr = GenomicRanges::GRanges(subs.data[,1], IRanges::IRanges(subs.data[,2], subs.data[,2]), rep('*', nrow(subs.data)))
  elementMetadata(subs.data.gr) = subs.data[,c(3,4)]
  
  alleleFrequencies = read.delim(allele_frequencies_file, sep='\t', header=T, quote=NULL, stringsAsFactors=F)
  alleleFrequencies = alleleFrequencies[order(alleleFrequencies[,1],alleleFrequencies[,2]),]
  print(head(alleleFrequencies))
  alleleFrequencies.gr = GenomicRanges::GRanges(alleleFrequencies[,1], IRanges::IRanges(alleleFrequencies[,2], alleleFrequencies[,2]), rep('*', nrow(alleleFrequencies)))
  elementMetadata(alleleFrequencies.gr) = alleleFrequencies[,3:7]
  
  # Subset the allele frequencies by the loci we would like to include
  overlap = findOverlaps(subs.data.gr, alleleFrequencies.gr)
  alleleFrequencies = alleleFrequencies[subjectHits(overlap),]
  
  nucleotides = c("A","C","G","T")
  ref.indices = match(subs.data[,3],nucleotides)
  alt.indices = match(subs.data[,4],nucleotides)
  WT.count = as.numeric(unlist(sapply(1:nrow(alleleFrequencies),function(v,a,i){v[i,a[i]+2]},v=alleleFrequencies,a=ref.indices)))
  mut.count = as.numeric(unlist(sapply(1:nrow(alleleFrequencies),function(v,a,i){v[i,a[i]+2]},v=alleleFrequencies,a=alt.indices)))
  
  combined = data.frame(chr=subs.data[,1],pos=subs.data[,2],WTCount=WT.count, mutCount=mut.count)
  colnames(combined) = c("chr","pos","WT.count","mut.count")
  
  combined.gr = GenomicRanges::GRanges(seqnames(subs.data.gr), ranges(subs.data.gr), rep('*', nrow(subs.data)))
  elementMetadata(combined.gr) = data.frame(WT.count=WT.count, mut.count=mut.count)
  
  combined.gr = sortSeqlevels(combined.gr)
  return(combined.gr)
}

##############################################
# GetDirichletProcessInfo
##############################################
#' Create the DPClust input file
#' 
#' Function that takes allele counts and a copy number profile to estimate mutation copy number,
#' cancer cell fraction and multiplicity for each point mutation.
#' @param loci_file Simple four column file with chromosome, position, reference allele and alternative allele
#' @param allele_frequencies_file Output file from alleleCounter on the specified loci
#' @param cellularity_file Full path to a Battenberg rho_and_psi output file
#' @param subclone_file Full path to a Battenberg subclones.txt output file
#' @param gender Specify male or female
#' @param SNP.phase.file Output file from mut_mut_phasing, supply NA (as char) when not available
#' @param mut.phase.file Output file from mut_cn_phasing, supply NA (as char) when not available
#' @param output_file Name of the output file
#' @author sd11
#' @export
runGetDirichletProcessInfo = function(loci_file, allele_frequencies_file, cellularity_file, subclone_file, gender, SNP.phase.file, mut.phase.file, output_file) {
  if(gender == 'male' | gender == 'Male') {
    isMale = T
  } else if(gender == 'female' | gender == 'Female') {
    isMale = F
  } else {
    stop("Unknown gender supplied, exit.")
  }
  info_counts = GetWTandMutCount(loci_file, allele_frequencies_file)
  cellularity = GetCellularity(cellularity_file)
  GetDirichletProcessInfo(output_file, cellularity, info_counts, subclone_file, is.male=isMale, SNP.phase.file=SNP.phase.file, mut.phase.file=mut.phase.file)
}

##############################################
# dpIn to VCF
##############################################
#' DPClust input file to vcf
#' 
#' Transform a dirichlet input file into a VCF with the same info. It filters out mutations in areas that are not contained in the supplied genome index (fai file) or are contained in the ignore file (ign file)
#' It takes the DP input file created by runGetDirichletProcessInfo and combines the columns with the vcf file supplied. Finally it gzips and indexes the file
#' @param vcf_infile Filename of the VCF file to use as a base
#' @param dpIn_file Filename of a DP input file
#' @param vcf_outfile Filename of the output file
#' @param fai_file Path to a reference genome index containing chromosome names
#' @param ign_file Path to a file containing contigs to ignore
#' @param genome Specify the reference genome for reading in the VCF
#' @author sd11
#' @export
dpIn2vcf = function(vcf_infile, dpIn_file, vcf_outfile, fai_file, ign_file, genome="hg19") {
  vcf = VariantAnnotation::readVcf(vcf_infile, genome=genome)
  
  # Remove muts on chroms not to look at
  fai = parseFai(fai_file)
  ign = parseIgnore(ign_file)
  allowed_chroms = which(!(fai$chromosome %in% ign$chromosome))
  vcf = vcf[as.vector(seqnames(vcf)) %in% allowed_chroms,]
  
  # Read in the to be annotated data
  dat = read.table(dpIn_file, header=T, stringsAsFactors=F)
  
  # Annotate the columns into the VCF object  
  vcf = addVcfInfoCol(vcf, dat$WT.count, 1, "Integer", "Number of reads carrying the wild type allele", "WC")
  vcf = addVcfInfoCol(vcf, dat$mut.count, 1, "Integer", "Number of reads carrying the mutant allele", "MC")
  vcf = addVcfInfoCol(vcf, dat$subclonal.CN, 1, "Float", "Total subclonal copynumber", "TSC")
  vcf = addVcfInfoCol(vcf, dat$nMaj1, 1, "Float", "Copynumber of major allele 1", "NMA1")
  vcf = addVcfInfoCol(vcf, dat$nMin1, 1, "Float", "Copynumber of minor allele 1", "NMI1")
  vcf = addVcfInfoCol(vcf, dat$frac1, 1, "Float", "Fraction of tumour cells containing copy number state 1", "FR1")
  vcf = addVcfInfoCol(vcf, dat$nMaj2, 1, "Float", "Copynumber of major allele 2", "NMA2")
  vcf = addVcfInfoCol(vcf, dat$nMin2, 1, "Float", "Copynumber of minor allele 2", "NMI2")
  vcf = addVcfInfoCol(vcf, dat$frac2, 1, "Float", "Fraction of tumour cells containing copy number state 2", "FR2")
  vcf = addVcfInfoCol(vcf, dat$phase, 1, "String", "Phase relation mutation and copynumber", "PHS")
  vcf = addVcfInfoCol(vcf, dat$mutation.copy.number, 1, "Float", "Mutation copy number", "MCN")
  vcf = addVcfInfoCol(vcf, dat$subclonal.fraction, 1, "Float", "Fraction of tumour cells carying this mutation", "CCF")
  vcf = addVcfInfoCol(vcf, dat$no.chrs.bearing.mut, 1, "Float", "Number of chromosomes bearing the mutation", "NCBM")
  
  # Write the output, gzip and index it
  VariantAnnotation::writeVcf(vcf, file=vcf_outfile, index=T)
}

#' Convenience function that annotates a column into the supplied VCF object
#' @noRd
addVcfInfoCol = function(vcf, data, number, type, description, abbreviation) {
  i = header(vcf)@header$INFO
  exptData(vcf)$header@header$INFO <- rbind(i, S4Vectors::DataFrame(Number=number, Type=type, Description=description, row.names=abbreviation))
  info(vcf)[,abbreviation] <- as(data, "CharacterList")
  return(vcf)
}

##############################################
# Produce project master file
##############################################
#' Creates a DPClust project master file
#' 
#' @param outputfile Full path with filename where the output will be written
#' @param donornames A vector with donor identifiers, use the same donor identifier to match multiple samplenames for a multi-sample DPClust run
#' @param samplenames A vector with sample identifiers
#' @param purities A vector with a purity value per sample
#' @param sex A vector with the sex of each donor
#' @param datafiles Vector with filenames in which the DPClust input is contained (Default: [samplename]_allDirichletProcessInfo.txt)
#' @param cndatafiles A vector with CNA DPClust input files (Default: NULL)
#' @param indeldatafiles A vector with indel DPClust input files (Default: NULL)
#' @author sd11
#' @export
createProjectFile = function(outputfile, donornames, samplenames, sex, purities=NULL, rho_and_psi_files=NULL, datafiles=paste(samplenames, "_allDirichletProcessInfo.txt", sep=""), cndatafiles=NULL, indeldatafiles=NULL) {
  if (is.null(purities) & is.null(rho_and_psi_files)) {
    stop("Please provide either a vector of purities or a vector of Battenberg rho_and_psi files")
  }
  
  .checklength = function(a) {
    if (!is.null(a)) {
      if (length(a)!=length(donornames)) {
        stop("All provided vectors must be of the same length")
      }
    }
  }
  .checklength(donornames)
  .checklength(samplenames)
  .checklength(sex)
  .checklength(purities)
  .checklength(rho_and_psi_files)
  .checklength(datafiles)
  .checklength(indeldatafiles)
  .checklength(cndatafiles)
  
  if (!is.null(rho_and_psi_files)) {
    purities = unlist(lapply(rho_and_psi_files, GetCellularity))
  }
  
  output = data.frame(sample=donornames, 
                      subsample=samplenames, 
                      datafile=datafiles, 
                      cellularity=purities, 
                      sex=sex, 
                      cnadatafile=ifelse(!is.null(cndatafiles), cndatafiles, NA),
                      indeldatafiles=ifelse(!is.null(indeldatafiles), indeldatafiles, NA), stringsAsFactors=F)
  write.table(output, file=outputfile, quote=F, row.names=F, sep="\t")
}


