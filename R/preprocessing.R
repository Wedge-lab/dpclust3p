ALLELECOUNTER = "alleleCounter"
LINKAGEPULL = "Linkage_pull.pl"

#'
#' Concatenate a series of files specified in a file of file names.
#'
concat_files = function(fofn, inputdir, outfile, haveHeader) {
  
  inputdir = paste(inputdir,"/", sep="")
  list_of_files = read.table(fofn, stringsAsFactors=F, header=F)[,1]
  
  output = data.frame()
  for (infile in list_of_files) {
    infile = paste(inputdir, infile, sep="")
    # Check if file is there and it contains data
    if (file.exists(infile) & file.info(infile)$size != 0) {
      dat = read.delim(infile, header=haveHeader, quote=NULL, stringsAsFactors=F)
      output = rbind(output, dat)
    }
  }
  
  write.table(output, file=outfile, col.names=haveHeader, row.names=F, sep="\t", quote=F)
}

#'
#' Split an input file per chromosome
#'
split_by_chrom = function(infile, prefix, postfix, outdir, chrom_file) {
  outdir = paste(outdir, "/", sep="")
  
  loci = read.delim(infile, stringsAsFactors=F, header=F)
  
  chroms = read.delim(chrom_file, stringsAsFactors=F, header=F)
  
  for (i in 1:nrow(chroms)) {
    selection = loci[loci[,1]==chroms[i,1],]
    chrom_id = chroms[i,2]
    write.table(selection, file=paste(outdir, prefix, chrom_id, postfix, sep=""), quote=F, row.names=F, col.names=F, sep="\t")
  }
}

############################################
# VCF 2 LOCI
############################################
parseFai = function(fai_file) {
  fai = read.table(fai_file, header=F, stringsAsFactor=F)
  colnames(fai) = c("chromosome", "length", "offset", "fasta_line_length", "line_blen")
  return(fai)
}

parseIgnore = function(ignore_file) {
  ign = read.table(ignore_file, header=F, stringsAsFactor=F)
  colnames(ign) = c("chromosome")
  return(ign)
}

vcf2loci = function(vcf_files, fai_file, ign_file, outfile) {
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
  # Remove duplicate entries
  combined.loci = unique(combined.loci)
  write.table(combined.loci, col.names=F, quote=F, row.names=F, file=outfile, sep="\t")
}

############################################
# Allele counting
############################################

#'
#' Count the alleles for specified locations in the loci file.
#'
alleleCount = function(locifile, bam, outfile, min_baq=20, min_maq=35) {

   cmd = paste(ALLELECOUNTER,
               "-b", bam,
               "-o", outfile,
               "-l", locifile,
               "-m", min_baq,
               "-q", min_maq, sep=" ")
   system(cmd,wait=T)
}

#' Dump allele counts stored in the sample columns of the VCF file. Output will go into a file
#' supplied as tumour_outfile and optionally normal_outfile. It will be a fully formatted
#' allele counts file as returned by alleleCounter.
dumpCounts.Sanger = function(vcf_infile, tumour_outfile, normal_outfile=NA, refence_genome="hg19") {
  # Helper function for writing the output  
  write.output = function(output, output_file) {
    write.table(output, file=output_file, col.names=T, quote=F, row.names=F, sep="\t")
  }
  
  # Read in the vcf and dump the tumour counts in the right format
  v = readVcf(vcf_infile, refence_genome)
  tumour = getCountsTumour(v)
  tumour = formatOutput(tumour, v)
  write.output(tumour, tumour_outfile)
  
  # Optionally dump the normal counts in the right format
  if (!is.na(normal_outfile)) {
    normal = getCountsNormal(v)
    normal = formatOutput(normal, v)
    write.output(normal, normal_outfile)
  }
}

#' Format a 4 column counts table into the alleleCounter format
formatOutput = function(counts_table, v) {
  output = data.frame(as.character(seqnames(v)), start(ranges(v)), counts_table, rowSums(counts_table))
  colnames(output) = c("#CHR","POS","Count_A","Count_C","Count_G","Count_T","Good_depth")
  return(output)
}

#' Returns an allele counts table for the normal sample
getCountsNormal = function(v) {
  return(getAlleleCounts.Sanger(v, 1))
}

#' Returns an allele counts table for the tumour sample
getCountsTumour = function(v) {
  return(getAlleleCounts.Sanger(v, 2))
}

#' Helper function that dumps the allele counts from a Sanger pipeline VCF file
getAlleleCounts.Sanger = function(v, sample_col) {
  return(cbind(geno(v)$FAZ[,sample_col]+geno(v)$RAZ[,sample_col], geno(v)$FCZ[,sample_col]+geno(v)$RCZ[,sample_col], geno(v)$FGZ[,sample_col]+geno(v)$RGZ[,sample_col], geno(v)$FTZ[,sample_col]+geno(v)$RTZ[,sample_col]))
}

############################################
# MUT 2 MUT phasing
############################################
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


mut_mut_phasing = function(loci_file, phased_file, bam_file, bai_file, max_distance) {
  # Check if there are lines in the file, otherwise it will crash this script
  if (file.info(loci_file)$size == 0) {
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







        #if (linked.muts$AF[i]>0.5) {
        #  
        #  if (linked.muts[i,ACGT[linked.muts$Var2[i]]] > 0 & linked.muts[i,ACGT[linked.muts$Ref2[i]]] == 0) {
        #    linked.muts$Parental[i] = "MUT_ON_RETAINED"
        #  }  else {
        #    if (linked.muts[i,ACGT[linked.muts$Var2[i]]] == 0 & linked.muts[i,ACGT[linked.muts$Ref2[i]]] > 0) {
        #      linked.muts$Parental[i] = "MUT_ON_DELETED"
        #    }
        #  }
        #}
        #else {	
        #  if (linked.muts[i,ACGT[linked.muts$Var2[i]]] > 0 & linked.muts[i,ACGT[linked.muts$Ref2[i]]] == 0) {
        #    linked.muts$Parental[i] = "MUT_ON_DELETED"
        #  } else {
        #    if (linked.muts[i,ACGT[linked.muts$Var2[i]]] == 0 & linked.muts[i,ACGT[linked.muts$Ref2[i]]] > 0) {
        #      linked.muts$Parental[i] = "MUT_ON_RETAINED"
        #    }
        #  }
        #}
        # 
        #220212
        #if (linked.muts$AFphased[i] < 0.5) {
        #  print(paste(linked.muts$Pos2[i],phased$AFphased[phased$Pos == linked.muts$Pos2[i]],sep=": "))
        #  if (!is.na(linked.muts$Parental[i]) & linked.muts$Parental[i] == "MUT_ON_DELETED") {
        #    linked.muts$Parental[i] = "MUT_ON_RETAINED"
        #  } else {
        #    if (!is.na(linked.muts$Parental[i]) & linked.muts$Parental[i] == "MUT_ON_RETAINED") {
        #      linked.muts$Parental[i] = "MUT_ON_DELETED"
        #    }
        #  }
        #}
      }
    }
  }
  
  write.table(linked.muts,outfile, sep="\t", quote=F, row.names=F)
}

############################################
# Combine all the steps into a DP input file
############################################
GetDirichletProcessInfo<-function(outputfile, cellularity, info, subclone.file, is.male = F, out.dir = NULL, SNP.phase.file = NULL, mut.phase.file = NULL){
 
  subclone.data = read.table(subclone.file,sep="\t",header=T,stringsAsFactors=F)
  #   subclone.data$subclonal.CN = (subclone.data$nMaj1_A + subclone.data$nMin1_A) * subclone.data$frac1_A
  subclone.data.gr = GenomicRanges::GRanges(subclone.data$chr, IRanges::IRanges(subclone.data$startpos, subclone.data$endpos), rep('*', nrow(subclone.data)))
  elementMetadata(subclone.data.gr) = subclone.data[,3:ncol(subclone.data)]
  # 	info2 = as.data.frame(cbind(as.data.frame(info), array(NA, c(length(info), 7))))
  #   colnames(info2) = c('chr','pos','WT.count','mut.count','subclonal.CN','nMaj1','nMin1','frac1','nMaj2','nMin2','frac2')
  
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
  expected.burden.for.MCN = mutationCopyNumberToMutationBurden(rep(1,length(info)),info$subclonal.CN,cellularity)
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
    prop.test(v$mut.count[i],v$mut.count[i] + v$WT.count[i],alternative="less")$p.value
  },v=info[non.zero.indices,], e= expected.burden.for.MCN[non.zero.indices])
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
  
  possible.zero.muts = intersect((1:length(info))[-non.zero.indices],which(!is.na(info$nMin1)))
  possible.zero.muts = c(possible.zero.muts,non.zero.indices[p.vals2>0.05])
  if(length(possible.zero.muts)>0){
    del.indices = which(info$nMin1[possible.zero.muts]==0 & !info$phase[possible.zero.muts]=="MUT_ON_RETAINED")
    info$subclonal.fraction[possible.zero.muts[del.indices]] = NA
    info$no.chrs.bearing.mut[possible.zero.muts[del.indices]] = 0
  }
  
  # convert GenomicRanges object to df
  df = data.frame(chr=as.data.frame(seqnames(info)),
                  start=start(info)-1,
                  end=end(info))
  df = cbind(df, as.data.frame(elementMetadata(info)))
  colnames(df)[1] = "chr"
  print(head(df))
  write.table(df, outputfile, sep="\t", row.names=F, quote=F)
}

GetCellularity <- function(rho_and_psi_file) {
  d = read.table(rho_and_psi_file, header=T, stringsAsFactors=F)
  return(d['FRAC_GENOME','rho'])
}

GetWTandMutCount <- function(loci_file, allele_frequencies_file) {
  #not.allowed.chroms = c("GL000191.1","GL000192.1","GL000193.1","GL000194.1","GL000195.1","GL000196.1","GL000197.1","GL000198.1","GL000199.1","GL000200.1","GL000201.1","GL000202.1","GL000203.1","GL000204.1","GL000205.1","GL000206.1","GL000207.1","GL000208.1","GL000209.1","GL000210.1","GL000211.1","GL000212.1","GL000213.1","GL000214.1","GL000215.1","GL000216.1","GL000217.1","GL000218.1","GL000219.1","GL000220.1","GL000221.1","GL000222.1","GL000223.1","GL000224.1","GL000225.1","GL000226.1","GL000227.1","GL000228.1","GL000229.1","GL000230.1","GL000231.1","GL000232.1","GL000233.1","GL000234.1","GL000235.1","GL000236.1","GL000237.1","GL000238.1","GL000239.1","GL000240.1","GL000241.1","GL000242.1","GL000243.1","GL000244.1","GL000245.1","GL000246.1","GL000247.1","GL000248.1","GL000249.1","hs37d5","MT","NC_007605")
  subs.data = read.table(loci_file, sep='\t', header=F, stringsAsFactors=F)
  subs.data = subs.data[order(subs.data[,1], subs.data[,2]),]
  #subs.data = subs.data[!(subs.data[,1] %in% not.allowed.chroms),]
  
  # Replace dinucleotides and longer with just the first base. Here we assume the depth of the second base is the same and the number of dinucleotides is so low that removing the second base is negligable
  subs.data[,3] = apply(as.data.frame(subs.data[,3]), 1, function(x) { substring(x, 1,1) })
  subs.data[,4] = apply(as.data.frame(subs.data[,4]), 1, function(x) { substring(x, 1,1) })
  
  subs.data.gr = GenomicRanges::GRanges(subs.data[,1], IRanges::IRanges(subs.data[,2], subs.data[,2]), rep('*', nrow(subs.data)))
  elementMetadata(subs.data.gr) = subs.data[,c(3,4)]
  
  #alleleFrequencies = read.table(allele_frequencies_file, sep='\t',header=F, stringsAsFactors=F)
  alleleFrequencies = read.delim(allele_frequencies_file, sep='\t', header=T, quote=NULL, stringsAsFactors=F)
  alleleFrequencies = alleleFrequencies[order(alleleFrequencies[,1],alleleFrequencies[,2]),]
  print(head(alleleFrequencies))
  alleleFrequencies.gr = GenomicRanges::GRanges(alleleFrequencies[,1], IRanges::IRanges(alleleFrequencies[,2], alleleFrequencies[,2]), rep('*', nrow(alleleFrequencies)))
  elementMetadata(alleleFrequencies.gr) = alleleFrequencies[,3:7]

  nucleotides = c("A","C","G","T")
  ref.indices = match(subs.data[,3],nucleotides)
  alt.indices = match(subs.data[,4],nucleotides)
  WT.count = as.numeric(unlist(sapply(1:nrow(alleleFrequencies),function(v,a,i){v[i,a[i]+2]},v=alleleFrequencies,a=ref.indices)))
  mut.count = as.numeric(unlist(sapply(1:nrow(alleleFrequencies),function(v,a,i){v[i,a[i]+2]},v=alleleFrequencies,a=alt.indices)))

  combined = data.frame(chr=subs.data[,1],pos=subs.data[,2],WTCount=WT.count, mutCount=mut.count)
  colnames(combined) = c("chr","pos","WT.count","mut.count")

  #   combined.gr = GenomicRanges::GRanges(subs.data[,1], IRanges::IRanges(subs.data[,2], subs.data[,2]+1), rep('*', nrow(subs.data)))
  combined.gr = GenomicRanges::GRanges(seqnames(subs.data.gr), ranges(subs.data.gr), rep('*', nrow(subs.data)))
  elementMetadata(combined.gr) = data.frame(WT.count=WT.count, mut.count=mut.count)
  return(combined.gr)
}

##############################################
# GetDirichletProcessInfo
##############################################
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
#' Transform a dirichlet input file into a VCF with the same info. It filters out mutations in areas that are not contained in the supplied genome index (fai file) or are contained in the ignore file (ign file)
#' It takes the DP input file created by runGetDirichletProcessInfo and combines the columns with the vcf file supplied. Finally it gzips and indexes the file
#' @param vcf_infile Filename of the VCF file to use as a base
#' @param dpIn_file Filename of a DP input file
#' @param vcf_outfile Filename of the output file
#' @param fai_file Path to a reference genome index containing chromosome names
#' @param ign_file Path to a file containing contigs to ignore
#' @param genome Specify the reference genome for reading in the VCF
#' @author Stefan Dentro
dpIn2vcf = function(vcf_infile, dpIn_file, vcf_outfile, fai_file, ign_file, genome="hg19") {
  vcf = readVcf(vcf_infile, genome=genome)
  
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
  writeVcf(vcf, file=vcf_outfile, index=T)
}

#' Convenience function that annotates a column into the supplied VCF object
addVcfInfoCol = function(vcf, data, number, type, description, abbreviation) {
  i = header(vcf)@header$INFO
  exptData(vcf)$header@header$INFO <- rbind(i, DataFrame(Number=number, Type=type, Description=description, row.names=abbreviation))
  info(vcf)[,abbreviation] <- as(data, "CharacterList")
  return(vcf)
}

##############################################
# Copy number preparations
##############################################
#' This function takes a Battenberg subclones.txt file, pulls out the A copy number profile,
#' then it calculates the total CN and a weighted ploidy before it subsets the input data
#' by columns: "chr", "startpos", "endpos", "nMaj1_A", "nMin1_A", "frac1_A", 
#' "nMaj2_A", "nMin2_A", "frac2_A", "SDfrac_A". 
#' 
#' Then logic is applied for classifying copy number aberrations into
#' a series of classes. There are two main groups: Clonal (starting with c) and
#' Subclonal (starting with s).
#' 
#' A CNA can be:
#' * HD : Homozygous deletion
#' * LOH : Loss of heterozygosity
#' * Amp : Amplification
#' * Gain : Gain
#' * NoCNV : Normal copy number without aberrations
#' * Loss : Loss
#' 
#' The copy number segments, classification, weighted ploidy and samplename are saved
#' to the specified output file
#' 
#' @param samplename, A String containing the samplename
#' @param bb_subclones_file, A String containing the full path to a Battenberg subclones.txt file
#' @param is_male, Boolean that should be TRUE if this sample is male, FALSE otherwise
#' @param outfile, A String with the full path to where the output should be written
#' @author Tom Mitchell
collate_bb_subclones = function(samplename, bb_subclones_file, is_male, outfile) {
  allsegs = NULL
  cndata = read.table(bb_subclones_file, header=T, stringsAsFactors=F)
  cndata = cbind(samplename, cndata)
  # Remove the X chromosome when sample is male as Battenberg cannot infer copy number there
  if (is_male) {
    cndata = cndata[!cndata$chr=="X",]
  }
  
  # Obtain the total copy number for each segment
  totalcn = rep(0, dim(cndata)[1])
  for (k in 1:nrow(cndata)) {
    # Clonal copy number
    if (cndata$frac1_A[k] == 1) {
      totalcn[k] = cndata$nMaj1_A[k] + cndata$nMin1_A[k]
    # Subclonal copy number - which is represented as a mixture of two states
    } else if (cndata$frac1_A[k]!="Inf" & cndata$frac1_A[k]!="-Inf") {
      totalcn[k] = (cndata$nMaj1_A[k] + cndata$nMin1_A[k]) * cndata$frac1_A[k] + (cndata$nMaj2_A[k] + cndata$nMin2_A[k]) * cndata$frac2_A[k]
    }
    cndataweighted = totalcn * (cndata$endpos - cndata$startpos) / (sum(as.numeric(cndata$endpos-cndata$startpos)))
  }
  # Ploidy weighted by segment size
  ploidy = round(sum(cndataweighted)/2)*2
  
  # Collating logical arguments for LOH, Amp, HD etc
  allsegs = data.frame(cndata[,1:4], cndata[,9:15], ploidy)
  names(allsegs) = c("Tumour_Name","chr", "startpos",
                       "endpos", "nMaj1_A", "nMin1_A",
                       "frac1_A", "nMaj2_A", "nMin2_A", "frac2_A",
                       "SDfrac_A", "tumour_ploidy")
  
  # Now classify all segments into a category
  tot = dim(allsegs)[1]
  # if you have clonal LOH, subclonal LOH is not counted!!!  Losses not counted
  allsegsa <- NULL
  CNA <- NULL
  for (i in 1:dim(allsegs)[1]) {
    # Clonal copy number
    if (is.na(allsegs$nMaj2_A[i])) {
      if (allsegs$nMin1_A[i] == 0 & allsegs$nMaj1_A[i] == 0) {
        allsegsa <- rbind(allsegsa, allsegs[i,c(1:7,11:12)])
        CNA <- c(CNA, "cHD")
      }
      if (xor(allsegs$nMin1_A[i] == 0, allsegs$nMaj1_A[i] == 0)) {
        allsegsa <- rbind(allsegsa, allsegs[i,c(1:7,11:12)])
        CNA <- c(CNA, "cLOH")
      }
      if ((allsegs$nMaj1_A[i] > allsegs$tumour_ploidy[i]/2) |
            (allsegs$nMin1_A[i] > allsegs$tumour_ploidy[i]/2)) {
        if ((allsegs$nMaj1_A[i] > allsegs$tumour_ploidy[i]*2) |
              (allsegs$nMin1_A[i] > allsegs$tumour_ploidy[i]*2)) {
          allsegsa <- rbind(allsegsa, allsegs[i,c(1:7,11:12)])
          CNA <- c(CNA, "cAmp")
        }
        else {
          allsegsa <- rbind(allsegsa, allsegs[i,c(1:7,11:12)])
          CNA <- c(CNA, "cGain")
        }
      }
      if ((allsegs$nMaj1_A[i] == allsegs$tumour_ploidy[i]/2) &
            (allsegs$nMin1_A[i] == allsegs$tumour_ploidy[i]/2)) {
        allsegsa <- rbind(allsegsa, allsegs[i,c(1:7,11:12)])
        CNA <- c(CNA, "NoCNV")
      }
      if ((allsegs$nMaj1_A[i] < allsegs$tumour_ploidy[i]/2) |
            (allsegs$nMin1_A[i] < allsegs$tumour_ploidy[i]/2)) {
        allsegsa <- rbind(allsegsa, allsegs[i,c(1:7,11:12)])
        CNA <- c(CNA, "cLoss")
      }
    }
    
    # Subclonal copy number
    if (!is.na(allsegs$nMaj2_A[i])) {
      if (allsegs$nMin1_A[i] == 0 & allsegs$nMaj1_A[i] == 0) {
        CNA <- c(CNA, "sHD")
        allsegsa <- rbind(allsegsa, allsegs[i,c(1:7,11:12)])
      }
      if ((allsegs$nMin1_A[i] == 0 | allsegs$nMaj1_A[i] == 0) &
            xor(allsegs$nMin2_A[i] == 0, allsegs$nMaj2_A[i] == 0)) {
        CNA <- c(CNA, "cLOH")
        tmp <- allsegs[i,c(1:4,8:12)]
        names(tmp) <- names(allsegs[i,c(1:7,11:12)])
        allsegsa <- rbind(allsegsa, tmp[,names(allsegs[i,c(1:7,11:12)])])
      }
      else if (xor(allsegs$nMin1_A[i] == 0, allsegs$nMaj1_A[i] == 0)) {
        CNA <- c(CNA, "sLOH")
        allsegsa <- rbind(allsegsa, allsegs[i,c(1:7,11:12)])
      }
      if ((allsegs$nMaj1_A[i] > allsegs$tumour_ploidy[i]/2 |
             allsegs$nMin1_A[i] > allsegs$tumour_ploidy[i]/2) &
            (allsegs$nMaj2_A[i] > allsegs$tumour_ploidy[i]/2 |
               allsegs$nMin2_A[i] > allsegs$tumour_ploidy[i]/2)) {
        if ((allsegs$nMaj1_A[i] > allsegs$tumour_ploidy[i]*2 |
               allsegs$nMin1_A[i] > allsegs$tumour_ploidy[i]*2) &
              (allsegs$nMaj2_A[i] > allsegs$tumour_ploidy[i]*2 |
                 allsegs$nMin2_A[i] > allsegs$tumour_ploidy[i]*2)) {
          CNA <- c(CNA, "cAmp")
          allsegsa <- rbind(allsegsa, allsegs[i,c(1:7,11:12)])
        }
        else {
          CNA <- c(CNA, "cGain")
          allsegsa <- rbind(allsegsa, allsegs[i,c(1:7,11:12)])
        }
      }
      if ((allsegs$nMin2_A[i] > allsegs$tumour_ploidy[i]/2 & allsegs$nMin1_A[i] != allsegs$nMin2_A[i])|
            (allsegs$nMaj2_A[i] > allsegs$tumour_ploidy[i]/2 & allsegs$nMaj1_A[i] != allsegs$nMaj2_A[i])) {
        if ((allsegs$nMin2_A[i] > allsegs$tumour_ploidy[i]*2 & allsegs$nMin1_A[i] != allsegs$nMin2_A[i]) |
              (allsegs$nMaj2_A[i] > allsegs$tumour_ploidy[i]*2 & allsegs$nMaj1_A[i] != allsegs$nMaj2_A[i])) {
          CNA <- c(CNA, "sAmp")
          tmp <- allsegs[i,c(1:4,8:12)]
          names(tmp) <- names(allsegs[i,c(1:7,11:12)])
          allsegsa <- rbind(allsegsa, tmp[,names(allsegs[i,c(1:7,11:12)])])
        }
        else {
          CNA <- c(CNA, "sGain")
          tmp <- allsegs[i,c(1:4,8:12)]
          names(tmp) <- names(allsegs[i,c(1:7,11:12)])
          allsegsa <- rbind(allsegsa, tmp[,names(allsegs[i,c(1:7,11:12)])])
        }
      }
      if ((allsegs$nMaj1_A[i] < allsegs$tumour_ploidy[i]/2 |
             allsegs$nMin1_A[i] < allsegs$tumour_ploidy[i]/2) &
            (allsegs$nMaj2_A[i] < allsegs$tumour_ploidy[i]/2 |
               allsegs$nMin2_A[i] < allsegs$tumour_ploidy[i]/2)) {
        CNA <- c(CNA, "cLoss")
        tmp <- allsegs[i,c(1:4,8:12)]
        names(tmp) <- names(allsegs[i,c(1:7,11:12)])
        allsegsa <- rbind(allsegsa, tmp[,names(allsegs[i,c(1:7,11:12)])])
      }
      if (((allsegs$nMaj1_A[i] < allsegs$tumour_ploidy[i]/2 & allsegs$nMaj1_A[i] != allsegs$nMaj2_A[i]) |
             (allsegs$nMin1_A[i] < allsegs$tumour_ploidy[i]/2 & allsegs$nMin1_A[i] != allsegs$nMin2_A[i]))) {
        allsegsa <- rbind(allsegsa, allsegs[i,c(1:7,11:12)])
        CNA <- c(CNA, "sLoss")      
      }
    }
    if(i %% 100 ==0){
      print(paste(i,"/",tot))
    }
  }  
  
  allsegsa <- cbind(allsegsa, CNA)
  allsegsa$frac1_A[allsegsa$CNA == "cGain"|allsegsa$CNA == "cAmp"|allsegsa$CNA == "cLOH"|allsegsa$CNA == "cHD"|allsegsa$CNA == "cLoss"] = 1
  
  write.table(allsegsa, file=outfile, sep="\t", quote=F, row.names=F)
}

##############################################
# QC
##############################################
createPng <- function(p, filename, width, height) {
  png(filename=filename, width=width, height=height)
  print(p)
  dev.off()
}

meltFacetPlotData <- function(data, subsamplenames) {
  d = as.data.frame(data)
  colnames(d) = subsamplenames
  data.m = reshape2::melt(d)
  return(data.m)
}

createHistFacetPlot <- function(data, title, xlab, ylab, binwidth) {
  p = ggplot2::ggplot(data) + ggplot2::aes(x=value, y=..count..) + ggplot2::geom_histogram(binwidth=binwidth, colour="black", fill="gray") + ggplot2::facet_grid(variable ~ .)
  p = p + ggplot2::theme_bw(base_size=35) + ggplot2::ggtitle(title) + ggplot2::xlab(xlab) + ggplot2::ylab(ylab)
  return(p)
}

createBoxFacetPlot <- function(data, title, xlab, ylab) {
  p = ggplot2::ggplot(data) + ggplot2::aes(x=variable, y=value) + ggplot2::geom_boxplot() + ggplot2::facet_grid(subsample ~ .)
  p = p + ggplot2::theme_bw() + ggplot2::ggtitle(title) + ggplot2::xlab(xlab) + ggplot2::ylab(ylab)
  return(p)
}

createQCDocument <- function(res, samplename, subsamplenames, outpath, cellularity) {
  p = createHistFacetPlot(meltFacetPlotData(res$totalCopyNumber, subsamplenames), paste(samplename, "totalCopyNumber"), "Copynumber", "Count", binwidth=1)
  createPng(p, paste(outpath, samplename, "_totalCopyNumber.png", sep=""), width=1500, height=500*length(subsamplenames))
  
  p = createHistFacetPlot(meltFacetPlotData(res$mutation.copy.number, subsamplenames), paste(samplename, "mutation.copy.number"), "mutation.copy.number", "Count", binwidth=0.1)
  p = p + ggplot2::xlim(0,5)
  createPng(p, paste(outpath, samplename, "_mutation.copy.number.png", sep=""), width=1500, height=500*length(subsamplenames))
  
  p = createHistFacetPlot(meltFacetPlotData(res$copyNumberAdjustment, subsamplenames), paste(samplename, "copyNumberAdjustment"), "copyNumberAdjustment", "Count (log10)", binwidth=1)
  p = p + ggplot2::scale_y_log10()
  createPng(p, paste(outpath, samplename, "_copyNumberAdjustment.png", sep=""), width=1500, height=500*length(subsamplenames))
  
  p = createHistFacetPlot(meltFacetPlotData(res$mutCount/(res$mutCount+res$WTCount), subsamplenames), paste(samplename, "alleleFrequency"), "Allele Frequency", "Count", binwidth=0.01)
  createPng(p, paste(outpath, samplename, "_alleleFrequency.png", sep=""), width=1500, height=500*length(subsamplenames))
  
  p = createHistFacetPlot(meltFacetPlotData(res$kappa, subsamplenames), paste(samplename, "kappa"), "Kappa", "Count", binwidth=0.01)
  createPng(p, paste(outpath, samplename, "_kappa.png", sep=""), width=1500, height=500*length(subsamplenames))
  
  #   p = createHistFacetPlot(meltFacetPlotData(res$subclonal.fraction, subsamplenames), paste(samplename, "subclonal.fraction"), "Subclonal fraction", "Count", binwidth=0.01)
  #   createPng(p, paste(outpath, samplename, "_subclonal.fraction.png", sep=""), width=1500, height=500*length(subsamplenames))
  
  p = createHistFacetPlot(meltFacetPlotData((res$mutCount+res$WTCount), subsamplenames), paste(samplename, "depth"), "Depth", "Count", binwidth=5)
  createPng(p, paste(outpath, samplename, "_depth.png", sep=""), width=1500, height=500*length(subsamplenames))
  
  #   fractionOfCells = res$mutation.copy.number / res$copyNumberAdjustment
  #   meltFacetPlotData(fractionOfCells, subsamplenames)
  p = createHistFacetPlot(meltFacetPlotData(res$subclonal.fraction, subsamplenames), paste(samplename, "Fraction Of Cells"), "Fraction of Cells", "Count", binwidth=0.05)
  p = p + ggplot2::geom_vline(xintercept=0.5, colour="red", linetype="longdash", size=2) + ggplot2::geom_vline(xintercept=1.5, colour="red", linetype="longdash", size=2) + ggplot2::xlim(0,3) #scale_y_log10()
  createPng(p, paste(outpath, samplename, "_fractionOfCells.png", sep=""), width=1500, height=500*length(subsamplenames))
  
  #   manualMutCopyNum = mutationBurdenToMutationCopyNumber(res$mutCount/(res$mutCount+res$WTCount),res$totalCopyNumber ,cellularity)
  #   p = createHistFacetPlot(meltFacetPlotData(manualMutCopyNum, subsamplenames), paste(samplename, "manualMutCopyNum"), "Mutation copy number", "Count", binwidth=0.1)
  #   createPng(p, paste(outpath, samplename, "_manualMutCopyNum.png", sep=""), width=500, height=250*length(subsamplenames))
  #   
  #   manualFractionOfCells = mutationBurdenToMutationCopyNumber(res$mutCount/(res$mutCount+res$WTCount),res$totalCopyNumber ,cellularity) / res$copyNumberAdjustment
  #   p = createHistFacetPlot(meltFacetPlotData(manualFractionOfCells, subsamplenames), paste(samplename, "manualFractionOfCells"), "Fraction of Cells", "Count", binwidth=0.01)
  #   createPng(p, paste(outpath, samplename, "_manualFractionOfCells.png", sep=""), width=1500, height=500*length(subsamplenames))
  
  # Manually melt the data
  d.m = data.frame()
  for (i in 1:length(subsamplenames)) {
    d.loc = data.frame(subsample=subsamplenames[i],variable=factor(dataset$chromosome[,i]), value=dataset$subclonal.fraction[,i])
    d.m = rbind(d.m, d.loc)
  }
  p = createBoxFacetPlot(d.m, paste(samplename, "subclonal fraction per chrom"), "Chromosome", "Subclonal fraction")
  createPng(p, paste(outpath, samplename, "_subclonalFractionPerChromosome.png", sep=""), width=1500, height=500*length(subsamplenames))
  
  # Manually melt the data
  d = as.data.frame(dataset$chromosome)
  colnames(d) = subsamplenames
  d.m = data.frame()
  for (i in 1:ncol(d)) {
    d.loc = d[dataset$subclonal.fraction > 1.5,i]
    if(length(d.loc) > 0) { 
      d.m = rbind(d.m, data.frame(variable=colnames(d)[i], value=factor(d.loc)))
    }
  }
  
  if (nrow(d.m) > 0) {
    p = createHistFacetPlot(d.m, paste(samplename, "subclonal fraction > 1.5"), "Chromosome", "Count", binwidth=1)
    createPng(p, paste(outpath, samplename, "_large.subclonal.fraction.by.chrom.png", sep=""), width=1500, height=500*length(subsamplenames))
  }
}

#'
#' Run the QC on specified data
#'
runQc = function(infile, datpath, outpath) {
  sample2purity = read.table(infile, header=T, stringsAsFactors=F)
  for (samplename in unique(sample2purity$sample)) {
    print(samplename)
    datafiles = sample2purity[sample2purity$sample==samplename,]$datafile
    subsamples = sample2purity[sample2purity$sample==samplename,]$subsample
    cellularity = sample2purity[sample2purity$sample==samplename,]$cellularity
    if (file.exists(paste(datpath,datafiles, sep=""))) {
      dataset = load.data(datpath,"",datafiles, cellularity=cellularity, Chromosome="chr", position="end", WT.count="WT.count", mut.count="mut.count", subclonal.CN="subclonal.CN", no.chrs.bearing.mut="no.chrs.bearing.mut", mutation.copy.number="mutation.copy.number", subclonal.fraction="subclonal.fraction", data_file_suffix="")
      createQCDocument(dataset, samplename, subsamples, outpath, cellularity)
    }
  }
}
