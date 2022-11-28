


args = commandArgs(T)
samplename = args[1]
bam_file = args[2]
loci_file = args[3]
chr = args[4]
outdir = args[5]

minMapQ = 20

setwd(outdir)
library(Rsamtools)
library(VariantAnnotation)

get_allele_combination_counts = function(bam, dat_pair) {
  alleles = data.frame()
  
  # alleles = lapply(unique(bam$groupid), function(i) {
  #for (j in seq(1, (length(bam$seq)), 2)) {
  for (i in unique(bam$groupid)) {
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
       
    if (!is.na(rel_pos_snv1)) {
      base_snv1 = substr(as.character(bam$seq[read_snv1]), rel_pos_snv1, rel_pos_snv1)
    } else {  
      base_snv1 = NA
    }
    
    alleles = rbind(alleles, data.frame(snv1=base_snv1, read_1=read_1))
  }#)
  
  # alleles = do.call(rbind, alleles)
  alleles_m = na.omit(alleles)
  return(alleles_m)
}

print("START")
muts <- read.delim(loci_file, header=F, row.names=NULL, stringsAsFactors=F)
names(muts) = c("CHR","POSITION","WT","MT")
chr.names = unique(muts$CHR)
chr.muts = muts[muts$CHR==chr,]
output <- data.frame(Chr = vector(mode="character",length=0), Pos1 = vector(mode="numeric",length=0), Ref1 = vector(mode="character",length=0), Var1 = vector(mode="character",length=0), Pos2 = vector(mode="numeric",length=0), Ref2 = vector(mode="character",length=0), Var2 = vector(mode="character",length=0))
count.data = matrix(NA,nrow=dim(chr.muts)[1],ncol=2)


for (i in 1:dim(chr.muts)[1]) {
  
  skip_to_next = FALSE
  
  tryCatch({
    dat_pair = data.frame(chrom=chr.muts$CHR[i], start=chr.muts$POSITION[i], end=chr.muts$POSITION[i])
    dat_pair = makeGRangesFromDataFrame(dat_pair)    
    wt = chr.muts$WT[i]
    mut = chr.muts$MT[i]
    
    flag = scanBamFlag(isPaired=T, hasUnmappedMate=F, isDuplicate=F, isUnmappedQuery=F) #, isProperPair=T
    param <- ScanBamParam(which=dat_pair[1,], what=scanBamWhat(), flag=flag, mapqFilter=minMapQ) #, mapqFilter=minMapQ
    bamfile <- BamFile(bam_file, asMates=TRUE)
    bam <- scanBam(bamfile, param=param)[[1]]
    alleles = get_allele_combination_counts(bam, dat_pair)
  }, error=function(e){
    cat("ERROR :",conditionMessage(e), "\n")
    skip_to_next <<-  TRUE})

  if(skip_to_next){
    next
  }
  
  allele_pairs = table(factor(paste(alleles$snv1), levels=c(wt, mut)))
  count.data[i,1:2] = c(allele_pairs[[wt]], allele_pairs[[mut]])

}

chr.muts$wt = count.data[,1]
chr.muts$mut = count.data[,2]

chr.muts = chr.muts[which(!is.na(chr.muts[,1])),]
write.table(chr.muts,file=paste0(outdir,"Count_Rsamtools_",chr,"_Reads.txt"),sep="\t",quote=F,row.names=F)

print("DONE")