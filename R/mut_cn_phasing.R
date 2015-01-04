# args=commandArgs(TRUE)
# 
# loci_file = toString(args[1])
# phased_file = toString(args[2])
# hap_file = toString(args[3])
# bam_file = toString(args[4])
# bai_file = toString(args[5])
# outfile = toString(args[6])
# max_distance = as.integer(args[7])

LINKAGEPULL = "Linkage_pull.pl"

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


mut_cn_phasing = function(loci_file, phased_file, hap_file, bam_file, bai_file, outfile, max.distance) {
  
  if (file.info(loci_file)$size == 0) {
    linked.muts = data.frame(matrix(rep(NA, 13), nrow=1))
    colnames(linked.muts) = c("Chr","Pos1","Ref1","Var1","Pos2","Ref2","Var2","AF","AFphased","Num_linked_to_A","Num_linked_to_C","Num_linked_to_G","Num_linked_to_T")
    linked.muts = na.omit(linked.muts)
  } else if(nrow(read.delim(loci_file, header=T, stringsAsFactors=F, fill=T))==0) {
    linked.muts = data.frame(matrix(rep(NA, 13), nrow=1))
    colnames(linked.muts) = c("Chr","Pos1","Ref1","Var1","Pos2","Ref2","Var2","AF","AFphased","Num_linked_to_A","Num_linked_to_C","Num_linked_to_G","Num_linked_to_T")
    linked.muts = na.omit(linked.muts)
  } else {
    # TODO: is this filtering required when just supplying loci files from a single chromosome?
    chr.muts = read.delim(loci_file, header=T, stringsAsFactors=F, fill=T)
    names(chr.muts) = c("CHR","POSITION","REF_BASE","MUT_BASE")
    
    # Match phased SNPs and their haplotypes together
    phased = read.delim(phased_file, header=T, stringsAsFactors=F, quote="\"")
    #220212
    colnames(phased) = c("SNP", "Chr","Pos", "AF", "AFphased", "AFsegmented")
    
    # TODO: check that chromosomes are using the same names between loci and phased files
    phased = phased[phased$Chr %in% chr.muts$CHR,]
    
    hap.info = read.delim(hap_file, sep=" ", header=F, row.names=NULL, stringsAsFactors=F)
    colnames(hap.info) = c("SNP","dbSNP","pos","ref","mut","ref_count","mut_count")
    # get haplotypes that match phased heterozygous SNPs
    hap.info = hap.info[match(phased$Pos,hap.info$pos),]
    print(sum(hap.info[,6]==1))
    #220212
    phased$AF[hap.info$ref_count==1] = 1-phased$AF[hap.info$ref_count==1]
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
        if (linked.muts$AF[i]>0.5) {
          
          if (linked.muts[i,ACGT[linked.muts$Var2[i]]] > 0 & linked.muts[i,ACGT[linked.muts$Ref2[i]]] == 0) {
            linked.muts$Parental[i] = "MUT_ON_RETAINED"
          }	else {
            if (linked.muts[i,ACGT[linked.muts$Var2[i]]] == 0 & linked.muts[i,ACGT[linked.muts$Ref2[i]]] > 0) {
              linked.muts$Parental[i] = "MUT_ON_DELETED"
            }
          }
        }
        else {	
          if (linked.muts[i,ACGT[linked.muts$Var2[i]]] > 0 & linked.muts[i,ACGT[linked.muts$Ref2[i]]] == 0) {
            linked.muts$Parental[i] = "MUT_ON_DELETED"
          } else {
            if (linked.muts[i,ACGT[linked.muts$Var2[i]]] == 0 & linked.muts[i,ACGT[linked.muts$Ref2[i]]] > 0) {
              linked.muts$Parental[i] = "MUT_ON_RETAINED"
            }
          }
        }
        
        #220212
        if (linked.muts$AFphased[i] < 0.5) {
          print(paste(linked.muts$Pos2[i],phased$AFphased[phased$Pos == linked.muts$Pos2[i]],sep=": "))
          if (!is.na(linked.muts$Parental[i]) & linked.muts$Parental[i] == "MUT_ON_DELETED") {
            linked.muts$Parental[i] = "MUT_ON_RETAINED"
          } else {
            if (!is.na(linked.muts$Parental[i]) & linked.muts$Parental[i] == "MUT_ON_RETAINED") {
              linked.muts$Parental[i] = "MUT_ON_DELETED"
            }
          }
        }
      }
    }
  }
  
  write.table(linked.muts,outfile, sep="\t", quote=F, row.names=F)
}

