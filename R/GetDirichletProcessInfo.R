# library(GenomicRanges)
# 
# args=commandArgs(TRUE)
# #                                                  Battenberg output                                GetAlleleFrequencies input  Battenberg output                               GetAlleleFrequencis output
# # R CMD BATCH '--no-restore-data --args PD7404a ../../battenberg/PD7404a/PD7404a_rho_and_psi.txt ../mutation_loci/PD7404a.loci ../../battenberg/PD7404a/PD7404a_subclones.txt PD7404a_alleleFrequencies.txt female' ~/repo/dirichlet/GetDirichletProcessInfo.R PD7404a.Rout
# 
# # First: /lustre/scratch110/sanger/dw9/life_history/PhaseSNPalleles.R (uses *BAFsegmented.txt and *_allHaplotypeInfo.txt to produce *chr17_amplifiedAlleles.txt ??)
# # *BAFsegmented.txt and *_allHaplotypeInfo.txt are BB output
# # Phase mut and CN: /lustre/scratch110/sanger/dw9/life_history/SNPs_analysis_PD9769.R (uses *BAFsegmented.txt and *_allHaplotypeInfo.txt *_muts_linkedMuts_segmented_*, which is read into this file)
# # Data here: /nfs/team78pc12/dw9/haplotyped_data/Lucy/    PD9769_17q_muts.txt
# 
# # samplename = toString(args[1])
# cellularity_file = toString(args[1])
# loci_file = toString(args[2])
# subclone_file = toString(args[3])
# allele_frequencies_file = toString(args[4])
# gender = toString(args[5])
# snp_phase_file = toString(args[6])
# mut_phase_file = toString(args[7])
# output_file = toString(args[8])

# wd = getwd()
# setwd("~/repo/dirichlet/dp_combined")
# source("interconvertMutationBurdens.R")
# setwd(wd)

# if(gender == 'male' | gender == 'Male') {
#   isMale = T
# } else if(gender == 'female' | gender == 'Female') {
#   isMale = F
# } else {
#   stop("Unknown gender supplied, exit.")
# }

##################################
# See script below these functions
##################################

.GetDirichletProcessInfo<-function(outputfile, cellularity, info, subclone.file, is.male = F, out.dir = NULL, phase.dir = NULL, SNP.phase.file = NULL, mut.phase.file = NULL){
  
  subclone.data = read.table(subclone.file,sep="\t",header=T,row.names=1, stringsAsFactors=F)
  # 	subclone.data$subclonal.CN = (subclone.data$nMaj1_A + subclone.data$nMin1_A) * subclone.data$frac1_A
  subclone.data.gr = GenomicRanges::GRanges(subclone.data$chr, IRanges(subclone.data$startpos, subclone.data$endpos), rep('*', nrow(subclone.data)))
  elementMetadata(subclone.data.gr) = subclone.data[,3:ncol(subclone.data)]
  # 	info2 = as.data.frame(cbind(as.data.frame(info), array(NA, c(length(info), 7))))
  #   colnames(info2) = c('chr','pos','WT.count','mut.count','subclonal.CN','nMaj1','nMin1','frac1','nMaj2','nMin2','frac2')
  
  info_anno = as.data.frame(cbind(array(NA, c(length(info), 7)))) 
  colnames(info_anno) = c('subclonal.CN','nMaj1','nMin1','frac1','nMaj2','nMin2','frac2')
  inds = findOverlaps(info, subclone.data.gr)  
  info_anno[queryHits(inds),2:7] = subclone.data[subjectHits(inds),][,c("nMaj1_A", "nMin1_A", "frac1_A", "nMaj2_A", "nMin2_A", "frac2_A")]
  print("Before CN")
  CN1 = (info_anno[queryHits(inds),]$nMaj1 + info_anno[queryHits(inds),]$nMin1) * info_anno[queryHits(inds),]$frac1
  # If frac is not one for allele 1 (i.e. not only CN data for allele 1), add the CN contribution of allele 2 as well
  CN2 = (info_anno[queryHits(inds),]$nMaj2 + info_anno[queryHits(inds),]$nMin2) * info_anno[queryHits(inds),]$frac2 * ifelse(info_anno[queryHits(inds),]$frac1 != 1, 1, 0)
  CN2[is.na(CN2)] = 0
  info_anno[queryHits(inds),]$subclonal.CN = CN1 + CN2
  elementMetadata(info) = cbind(elementMetadata(info), info_anno)
  
  info$phase="unphased"
  # Code to turn on when haplotype pipeline is fully operational
  phasing = read.table(SNP.phase.file, header=T, stringsAsFactors=F) #header=T, skip=1, 
  print(head(phasing))
  print(phasing$Pos1[!is.numeric(phasing$Pos1)])
  phasing.gr = GenomicRanges::GRanges(phasing$Chr, IRanges(phasing$Pos1, phasing$Pos1))
  phasing.gr$phasing = phasing$Parental
  inds = findOverlaps(info, phasing.gr)  
  info$phase[queryHits(inds)] = phasing.gr$phasing[subjectHits(inds)]
  
  info$phase[is.na(info$phase)]="unphased"
  print("M.C.N.")
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
  df = data.frame(chr=seqnames(info),
                  start=start(info)-1,
                  end=end(info))
  #   df = data.frame(chr=seqnames(info),
  #                 pos=start(info))
  df = cbind(df, elementMetadata(info))
  
  #   out.file = paste(samplename,"_allDirichletProcessInfo.txt",sep="")
  # 	if(!is.null(out.dir)){
  # 		out.file = paste(out.dir,'/',out.file,sep="")
  # 	}
  # 	write.table(df, out.file, sep="\t", row.names=F, quote=F)
  write.table(df, outputfile, sep="\t", row.names=F, quote=F)
}

.GetCellularity <- function(rho_and_psi_file) {
  d = read.table(rho_and_psi_file, header=T, stringsAsFactors=F)
  return(d['FRAC_GENOME','rho'])
}

.GetWTandMutCount <- function(loci_file, allele_frequencies_file) {
  not.allowed.chroms = c("GL000191.1","GL000192.1","GL000193.1","GL000194.1","GL000195.1","GL000196.1","GL000197.1","GL000198.1","GL000199.1","GL000200.1","GL000201.1","GL000202.1","GL000203.1","GL000204.1","GL000205.1","GL000206.1","GL000207.1","GL000208.1","GL000209.1","GL000210.1","GL000211.1","GL000212.1","GL000213.1","GL000214.1","GL000215.1","GL000216.1","GL000217.1","GL000218.1","GL000219.1","GL000220.1","GL000221.1","GL000222.1","GL000223.1","GL000224.1","GL000225.1","GL000226.1","GL000227.1","GL000228.1","GL000229.1","GL000230.1","GL000231.1","GL000232.1","GL000233.1","GL000234.1","GL000235.1","GL000236.1","GL000237.1","GL000238.1","GL000239.1","GL000240.1","GL000241.1","GL000242.1","GL000243.1","GL000244.1","GL000245.1","GL000246.1","GL000247.1","GL000248.1","GL000249.1","hs37d5","MT","NC_007605")
  subs.data = read.table(loci_file, sep='\t', header=F, stringsAsFactors=F)
  subs.data = subs.data[!(subs.data[,1] %in% not.allowed.chroms),]
  
  # Replace dinucleotides and longer with just the first base. Here we assume the depth of the second base is the same and the number of dinucleotides is so low that removing the second base is negligable
  subs.data[,3] = apply(as.data.frame(subs.data[,3]), 1, function(x) { substring(x, 1,1) })
  subs.data[,4] = apply(as.data.frame(subs.data[,4]), 1, function(x) { substring(x, 1,1) })
  
  subs.data.gr = GenomicRanges::GRanges(subs.data[,1], IRanges(subs.data[,2], subs.data[,2]), rep('*', nrow(subs.data)))
  elementMetadata(subs.data.gr) = subs.data[,c(3,4)]
  
  alleleFrequencies = read.table(allele_frequencies_file, sep='\t',header=T, stringsAsFactors=F)
  print(head(alleleFrequencies))
  alleleFrequencies.gr = GenomicRanges::GRanges(alleleFrequencies[,1], IRanges(alleleFrequencies[,2], alleleFrequencies[,2]), rep('*', nrow(alleleFrequencies)))
  elementMetadata(alleleFrequencies.gr) = alleleFrequencies[,3:7]
  print("1")
  ref.indices = match(subs.data[,3],nucleotides)
  alt.indices = match(subs.data[,4],nucleotides)
  print(head(sapply(1:nrow(alleleFrequencies),function(v,a,i){v[i,a[i]+2]},v=alleleFrequencies,a=ref.indices)))
  print(head(sapply(1:nrow(alleleFrequencies),function(v,a,i){v[i,a[i]+2]},v=alleleFrequencies,a=alt.indices)))
  WT.count = as.numeric(sapply(1:nrow(alleleFrequencies),function(v,a,i){v[i,a[i]+2]},v=alleleFrequencies,a=ref.indices))
  mut.count = as.numeric(sapply(1:nrow(alleleFrequencies),function(v,a,i){v[i,a[i]+2]},v=alleleFrequencies,a=alt.indices))
  print("2")  
  combined = data.frame(chr=subs.data[,1],pos=subs.data[,2],WTCount=WT.count, mutCount=mut.count)
  colnames(combined) = c("chr","pos","WT.count","mut.count")
  print("3")
  #   combined.gr = GenomicRanges::GRanges(subs.data[,1], IRanges(subs.data[,2], subs.data[,2]+1), rep('*', nrow(subs.data)))
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
  nucleotides = c("A","C","G","T")
  info_counts = GetWTandMutCount(loci_file, allele_frequencies_file)
  cellularity = GetCellularity(cellularity_file)
  GetDirichletProcessInfo(output_file, cellularity, info_counts, subclone_file, is.male=isMale, SNP.phase.file=snp_phase_file, mut.phase.file=mut_phase_file)
}
