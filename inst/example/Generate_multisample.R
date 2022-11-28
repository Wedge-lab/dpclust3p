
args = commandArgs(T)

samplenames = args[1]
count_path = args[2]
Battenberg_Path = args[3]
outdir = args[4]

library(tidyverse)
library(tidyr)
source("/gpfs3/well/woodcock/users/dus183/sim_data_Aug/script/interconvertMutationBurdens.R")

setwd(outdir)
#setwd("/gpfs3/well/woodcock/users/dus183/Real_sample/DPClust_input")


ac_dir = "./01_Allele_Count/"

if(!file.exists(ac_dir)){
  dir.create(ac_dir)
}

for(samplename in samplenames){
  if(file.exists(paste0(count_path,"/",samplename,"/Count_Rsamtools_1_Reads.txt"))){
    data = read.delim(paste0(count_path,"/",samplename,"/Count_Rsamtools_1_Reads.txt"),header=TRUE, stringsAsFactor=F)
    for(i in c(2:22)){
      data_new = read.delim(paste0(count_path,"/",samplename,"/Count_Rsamtools_",i,"_Reads.txt"),header=TRUE, stringsAsFactor=F)
      data = rbind(data,data_new)
    }
  }else{
    data = read.delim(paste0(count_path,"/",samplename,"/Count_Rsamtools_chr1_Reads.txt"),header=TRUE, stringsAsFactor=F)
    for(i in c(2:22)){
      data_new = read.delim(paste0(count_path,"/",samplename,"/Count_Rsamtools_chr",i,"_Reads.txt"),header=TRUE, stringsAsFactor=F)
      data = rbind(data,data_new)
    }
  }
  alt_reads = data[,6]
  ref_reads = data[,5]
  
  A = C = G = T = rep(0,dim(data)[1])
  A[which(data[,3]=="A")] = ref_reads[which(data[,3]=="A")]
  A[which(data[,4]=="A")] = alt_reads[which(data[,4]=="A")]
  
  C[which(data[,3]=="C")] = ref_reads[which(data[,3]=="C")]
  C[which(data[,4]=="C")] = alt_reads[which(data[,4]=="C")]
  
  G[which(data[,3]=="G")] = ref_reads[which(data[,3]=="G")]
  G[which(data[,4]=="G")] = alt_reads[which(data[,4]=="G")]
  
  T[which(data[,3]=="T")] = ref_reads[which(data[,3]=="T")]
  T[which(data[,4]=="T")] = alt_reads[which(data[,4]=="T")]
  
  D = A + C + G + T
  
  ID = rep(t, times= dim(data)[1])
  
  ac = cbind(as.vector(data[,1]), as.vector(data[,2]), ID, A, C, G, T, D, as.character(data[,3]), as.character(data[,4]))
  
  colnames(ac) = c("chr", "pos", "ID", "A", "C", "G", "T","depth","ref","alt")
  
  write.table(ac,
              file = paste0(ac_dir, t, "_ac.txt"),
              sep ='\t',
              col.names = c("chr", "pos", "ID", "A", "C", "G", "T","depth","ref","alt"),
              row.names = F,
              quote = F)
  
  m = as.data.frame(read.table(paste0(ac_dir, t, "_ac.txt"),header=TRUE,sep='\t',stringsAsFactors=F))
  
  colnames(m) = c("chr", "pos", "ID", "A", "C", "G", "T","depth","ref","alt")
  
  for (i in 1:nrow(m)){
    for (j in 1:ncol(m)){
      for (k in 1:ncol(m)){
        if (!is.na(match(m$ref[i],names(m)[j])) & !is.na(match(m$alt[i],names(m)[k]))){
          m$ref_count[i]=m[i,j]
          m$alt_count[i]=m[i,k]
          #print(paste(i))
        }
      }}}
  
  write.table(m,
              file = paste0(ac_dir, t, "_ac_match.txt"),
              sep ='\t',
              col.names = c("chr", "pos", "ID", "A", "C", "G", "T","depth","ref","alt","ref_count","alt_count"),
              row.names = F,
              quote = F)
  
  id = paste(m[,1], m[,2], sep="_")
  
  write.table(cbind(id, as.numeric(m[,"ref_count"])), file = paste0(ac_dir, t, "_WTCount.txt"), sep ='\t', col.names = c("ID", t), row.names = FALSE, quote= FALSE)
  
  write.table(cbind(id, as.numeric(m[,"alt_count"])), file = paste0(ac_dir, t, "_mutCount.txt"), sep ='\t', col.names = c("ID", t), row.names = FALSE, quote= FALSE)
  
}


#step2

nucleotides = c("A","C","G","T")

no.subsamples = length(samplenames)
  
for(i.sample in no.subsamples){
  
  print(samplenames[i.sample])
  
  ce= as.data.frame(read.table(paste0(Battenberg_Path, samplenames[i,samplename], "_rho_and_psi.txt"),header=TRUE,sep='\t'))
  cellularity = as.vector(as.numeric(ce[1,"rho"]))
  
  subclone.files = as.vector(paste0(Battenberg_Path, samplenames[i,samplename], "_subclones.txt"))
  
  INFO= as.data.frame(read.table(paste0("./01_allele_count/", samplenames[i,samplename], "_ac_match.txt"),header=TRUE,sep='\t'))[,1:2]
  
  info= INFO
  
  #colnames(INFO) = c("chr", "pos")
  
  mu=read.table(paste0("./01_allele_count/", samplenames[i,samplename], "_mutCount.txt"), sep='\t',header=TRUE,stringsAsFactors=F)
  
  mutCount= data.frame(mu[,-1])
  
  colnames(mutCount) = samplenames[i,samplename]
  
  w=read.table(paste0("./01_allele_count/", samplenames[i,samplename], "_WTCount.txt"), sep='\t',header=TRUE,stringsAsFactors=F)
  
  WTCount= data.frame(w[,-1])
  
  colnames(WTCount) = samplenames[i,samplename]
  
}

is.male = TRUE
out.files = NULL
phase.dir = NULL
SNP.phase.file = NULL
mut.phase.file = NULL
keep.muts.not.explained.by.CN=T

non.zero.indices=list()
p.vals2 = list()

if(is.null(out.files)){
  out.files = cbind(paste(samplename,"_ssDPI.txt",sep=""),paste(samplenames,"_ssDPI_WCL.txt",sep=""))
}

if(!file.exists("./02_Input_all_SNVs")){
  dir.create("./02_Input_all_SNVs")
}


data = list()

for(s in 1:no.subsamples){
  data[[s]] = data.frame(cbind(info,WTCount[,s],mutCount[,s]))
  names(data[[s]])[(ncol(data[[s]])-1):ncol(data[[s]])] =c("WT.count","mut.count")
  full.samplename = samplename[s]
  data[[s]]$pos=as.numeric(data[[s]]$pos) #naser: read as character and not used in mathematical equations below -> thus NA output		
  subclone.data = read.table(subclone.files[s],header=TRUE,stringsAsFactors=F)	
  subclone.data[,1] = substring(subclone.data[,1],4)

  data[[s]]$subclonal.CN = NA
  data[[s]]$nMaj1 = NA
  data[[s]]$nMin1 = NA
  data[[s]]$frac1 = NA
  data[[s]]$nMaj2 = NA
  data[[s]]$nMin2 = NA
  data[[s]]$frac2 = NA
  
  print("subclones file OK")
  
  for(r in 1:nrow(subclone.data)){
    CN = (subclone.data$nMaj1_A[r] + subclone.data$nMin1_A[r]) * subclone.data$frac1_A[r]
    if(subclone.data$frac1_A[r] != 1){
      CN = CN + (subclone.data$nMaj2_A[r] + subclone.data$nMin2_A[r]) * subclone.data$frac2_A[r]
    }
    data[[s]]$subclonal.CN[data[[s]][,1]==subclone.data$chr[r] & data[[s]][,2]>=subclone.data$startpos[r] & data[[s]][,2]<=subclone.data$endpos[r]] = CN
    data[[s]]$nMaj1[data[[s]][,1]==subclone.data$chr[r] & data[[s]][,2]>=subclone.data$startpos[r] & data[[s]][,2]<=subclone.data$endpos[r]] = subclone.data$nMaj1_A[r]
    data[[s]]$nMin1[data[[s]][,1]==subclone.data$chr[r] & data[[s]][,2]>=subclone.data$startpos[r] & data[[s]][,2]<=subclone.data$endpos[r]] = subclone.data$nMin1_A[r]
    data[[s]]$frac1[data[[s]][,1]==subclone.data$chr[r] & data[[s]][,2]>=subclone.data$startpos[r] & data[[s]][,2]<=subclone.data$endpos[r]] = subclone.data$frac1_A[r]
    data[[s]]$nMaj2[data[[s]][,1]==subclone.data$chr[r] & data[[s]][,2]>=subclone.data$startpos[r] & data[[s]][,2]<=subclone.data$endpos[r]] = subclone.data$nMaj2_A[r]
    data[[s]]$nMin2[data[[s]][,1]==subclone.data$chr[r] & data[[s]][,2]>=subclone.data$startpos[r] & data[[s]][,2]<=subclone.data$endpos[r]] = subclone.data$nMin2_A[r]
    data[[s]]$frac2[data[[s]][,1]==subclone.data$chr[r] & data[[s]][,2]>=subclone.data$startpos[r] & data[[s]][,2]<=subclone.data$endpos[r]] = subclone.data$frac2_A[r]	
  }
  
  data[[s]]$phase="unphased"
  #for(chr in unique(data[[s]][,1])){
  #	chr.inds = which(data[[s]][,1]==chr)
  #	ph.file = paste("phasing_data/",samplename,"_muts_linkedMuts_segmented_chr",chr,".txt",sep="")
  #	if(file.exists(ph.file)){
  #		phase.data[[s]] = read.table(ph.file,sep="\t",header=F,row.names=NULL,stringsAsFactors=F)
  #		indices = match(data[[s]][chr.inds,2],phase.data[[s]][,2])
  #		data[[s]]$phase[chr.inds[!is.na(indices)]] = phase.data[[s]][indices[!is.na(indices)],14]
  #	}
  #}
  #data[[s]]$phase[is.na(data[[s]]$phase)]="unphased"
  
  data[[s]]$normalCopyNumber = 2
  #assume that col1 contains chromosome
  data[[s]]$normalCopyNumber[data[[s]][,1]=="X"] = 1
  data[[s]]$mutation.copy.number = mutationBurdenToMutationCopyNumber(data[[s]]$mut.count/ (data[[s]]$mut.count + data[[s]]$WT.count) , data[[s]]$subclonal.CN, cellularity[s],data[[s]]$normalCopyNumber)
  
  data[[s]]$mutation.copy.number[is.nan(data[[s]]$subclonal.CN)]=NA
  data[[s]]$subclonal.CN[is.nan(data[[s]]$subclonal.CN)] = NA
  
  #convert MCN to subclonal fraction - tricky for amplified mutations
  data[[s]]$subclonal.fraction = data[[s]]$mutation.copy.number
  expected.burden.for.MCN = mutationCopyNumberToMutationBurden(rep(1,nrow(data[[s]])),data[[s]]$subclonal.CN,cellularity[s],data[[s]]$normalCopyNumber)
  non.zero.indices[[s]] = which(data[[s]]$mut.count>0 & !is.na(expected.burden.for.MCN))
  #test for mutations in more than 1 copy
  if(length(non.zero.indices[[s]])>0){
    p.vals = sapply(1:length(non.zero.indices[[s]]),function(v,e,i){prop.test(v$mut.count[i],v$mut.count[i] + v$WT.count[i],e[i],alternative="greater")$p.value},v=data[[s]][non.zero.indices[[s]],], e= expected.burden.for.MCN[non.zero.indices[[s]]])
    amplified.muts = non.zero.indices[[s]][p.vals<=0.05]
  }else{
    amplified.muts = integer(0)
  }
  
  data[[s]]$no.chrs.bearing.mut = 1	
  
  #copy numbers of subclones can only differ by 1 or 0 (as assumed when calling subclones)
  if(length(amplified.muts)>0){		
    for(a in 1:length(amplified.muts)){
      max.CN2=0
      #use phasing info - if on 'deleted' (lower CN chromosome), use the minor copy number
      if(data[[s]]$phase[amplified.muts[a]]=="MUT_ON_DELETED"){
        print("mut on minor chromosome")
        max.CN1 = data[[s]]$nMin1[amplified.muts[a]]
        frac1 = data[[s]]$frac1[amplified.muts[a]]
        frac2=0
        if(!is.na(data[[s]]$nMin2[amplified.muts[a]])){
          #swap subclones, so that the one with the higher CN is first
          if(data[[s]]$nMin2[amplified.muts[a]]>max.CN1){
            max.CN2 = max.CN1
            max.CN1 = data[[s]]$nMin2[amplified.muts[a]]
            frac2 = frac1
            frac1 = data[[s]]$frac2[amplified.muts[a]]
          }else{
            max.CN2 = data[[s]]$nMin2[amplified.muts[a]]
            frac2 = data[[s]]$frac2[amplified.muts[a]]
          }
        }					
      }else{
        max.CN1 = data[[s]]$nMaj1[amplified.muts[a]]
        frac1 = data[[s]]$frac1[amplified.muts[a]]
        frac2=0
        if(!is.na(data[[s]]$nMaj2[amplified.muts[a]])){
          #swap subclones, so that the one with the higher CN is first
          if(data[[s]]$nMaj2[amplified.muts[a]]>max.CN1){
            max.CN2 = max.CN1
            max.CN1 = data[[s]]$nMaj2[amplified.muts[a]]
            frac2 = frac1
            frac1 = data[[s]]$frac2[amplified.muts[a]]						
          }else{
            max.CN2 = data[[s]]$nMaj2[amplified.muts[a]]
            frac2 = data[[s]]$frac2[amplified.muts[a]]
          }
        }	
      }
      best.err = data[[s]]$mutation.copy.number[amplified.muts[a]] - 1
      best.CN=1
      for(j in 1:max.CN1){
        for(k in (j-1):min(j,max.CN2)){
          potential.CN = j * frac1 + k * frac2
          err = abs(data[[s]]$mutation.copy.number[amplified.muts[a]]/potential.CN-1)
          if(err<best.err){
            data[[s]]$no.chrs.bearing.mut[amplified.muts[a]] = potential.CN
            best.err=err
            best.CN = potential.CN
          }
        }
      }
      data[[s]]$subclonal.fraction[amplified.muts[a]] = data[[s]]$mutation.copy.number[amplified.muts[a]] / best.CN
    }
  }
  
  ##########################################################################
  
  #test for subclonal mutations
  #test whether mut burden is less than expected value for MCN = 1
  if(length(non.zero.indices[[s]])>0){
    p.vals1 = sapply(1:length(non.zero.indices[[s]]),function(v,e,i){prop.test(v$mut.count[i],v$mut.count[i] + v$WT.count[i],alternative="less")$p.value},v=data[[s]][non.zero.indices[[s]],], e= expected.burden.for.MCN[non.zero.indices[[s]]])
    #test whether mut burden is above error rate (assumed to be 1 in 200)
    p.vals2[[s]] = sapply(1:length(non.zero.indices[[s]]),function(v,i){prop.test(v$mut.count[i],v$mut.count[i] + v$WT.count[i],0.005,alternative="greater")$p.value},v=data[[s]][non.zero.indices[[s]],])
    
    #test for non-zero muts causes problems for muts that have been lost in most cells
    #subclonal.muts = non.zero.indices[[s]][p.vals1<=0.05 & p.vals2[[s]]<=0.05]
    subclonal.muts = non.zero.indices[[s]][p.vals1<=0.05]
  }else{
    subclonal.muts = integer(0)
  }
  
  # use subclonal CN that minimises the difference in subclonal fraction from 1
  if(length(subclonal.muts)>0){
    for(a in 1:length(subclonal.muts)){
      #if there are no subclonal CNVs, don't adjust subclonal fraction
      if(is.na(data[[s]]$frac2[subclonal.muts[a]])){next}
      #assume subclonal muts are on one chromosome copy, therefore mutation copy number must be subclonal fraction of the higher CN subclone (i.e. lost in lower CN subclone) or 1 (i.e. present in both subclones)
      if(data[[s]]$nMaj1[subclonal.muts[a]]+data[[s]]$nMin1[subclonal.muts[a]] > data[[s]]$nMaj2[subclonal.muts[a]]+data[[s]]$nMin2[subclonal.muts[a]]){	
        possible.subclonal.fractions = c(data[[s]]$frac1[subclonal.muts[a]],1)
      }else{
        possible.subclonal.fractions = c(data[[s]]$frac2[subclonal.muts[a]],1)
      }
      best.CN = possible.subclonal.fractions[which.min(abs(data[[s]]$mutation.copy.number[subclonal.muts[a]]/possible.subclonal.fractions - 1))]
      #extra test 200313 - check whether subclonal CN results in clonal mutation, otherwise subclonal CN doesn't explain subclonal MCN
      if(best.CN != 1 & prop.test(data[[s]]$mut.count[subclonal.muts[a]],data[[s]]$mut.count[subclonal.muts[a]]+data[[s]]$WT.count[subclonal.muts[a]],expected.burden.for.MCN[subclonal.muts[a]] * best.CN)$p.value > 0.05){
        data[[s]]$subclonal.fraction[subclonal.muts[a]] = data[[s]]$mutation.copy.number[subclonal.muts[a]] / best.CN
        data[[s]]$no.chrs.bearing.mut[subclonal.muts[a]] = best.CN
      }
    }
  }
  
  colnames(data[[s]])[2] = "end"
  write.table(data[[s]], paste0("./02_Input_all_SNVs/",out.files[s,1]),sep="\t",row.names=F,quote=F)
}
