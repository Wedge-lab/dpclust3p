args = commandArgs(T)

samplename = args[1]
phasing_path = args[2]
batt_path = args[3]
dp_input_path = args[4]
outdir = args[5]

copy.major = args[6]
copy.frac= args[7]

if (length(args) > 7) {
  copy.minor = args[8]
} else {
  copy.minor = NA
}


check_copy.statue = function(phasing.input,subclone.used){
  chr.temp = phasing.input[,1]
  loc.1 = phasing.input[,2]
  loc.2= phasing.input[,5]
  
  subclone.chr.index = which(subclone.used[,1]==chr.temp)
  
  if(length(which(subclone.used[,1]==chr.temp))==0){
    return("FALSE")
    break
  }
  
  subclone.seg = subclone.used[subclone.chr.index,]
  
  subclone.slice.index = which(subclone.seg[,2]<=loc.1 & subclone.seg[,3]>=loc.1)
  
  if(length(subclone.slice.index)==0){
    return("FALSE")
    break
  }
  
  seg.beg = subclone.seg[subclone.slice.index,2]
  seg.end = subclone.seg[subclone.slice.index,3]
  
  if(seg.beg<=loc.2 & seg.end>=loc.2){
    return("TRUE")
  }else{
    retrun("FALSE")
  }
}



subclone = read.delim(file.path(batt_path,paste0(samplename,"_subclones.txt")),header=TRUE,stringsAsFactors=F)

if(!is.na(copy.minor)){
  subclone.used.index = which(subclone[,8]==copy.major & subclone[,9]==copy.minor & subclone[,10]==copy.frac)
}else{
  subclone.used.index = which(subclone[,8]==copy.major & subclone[,10]==copy.frac)
}


subclone.used = subclone[subclone.used.index,]

phasing = read.delim(file.path(phasing_path,paste0(samplename,"_mutmut_phasing.txt")),header=TRUE,stringsAsFactors = F)


chr.loc = sapply(1:dim(phasing)[1],function(t){check_copy.statue(phasing[t,],subclone.used)})
cl.phasing.pass.index = which(chr.loc=="TRUE")

dp_input = read.delim(file.path(dp_input_path, paste0(samplename,"_ssDPI.txt")) ,header=TRUE,stringsAsFactors = F)


cl.phasing.pass = matrix(NA,nrow=length(cl.phasing.pass.index),ncol=17)
colnames(cl.phasing.pass)=c("chr","end_1","WT.count1","mut.count1","subclonal.fraction1",
                            "end_2","WT.count2","mut.count2","subclonal.fraction2",
                            "Num_WT_WT","Num_MUT_MUT","Num_WT_MUT","Num_MUT_WT","Copy_Major","Copy_Minor","Mut1_index","Mut2_index" )

print(cl.phasing.pass.index)
for(i.index in 1:length(cl.phasing.pass.index)){
  cl.phasing.pass[i.index,1:2] = unlist(phasing[cl.phasing.pass.index[i.index],1:2])
  cl.phasing.pass[i.index,6] = phasing[cl.phasing.pass.index[i.index],5]
  cl.phasing.pass[i.index,10:13] = unlist(phasing[cl.phasing.pass.index[i.index],8:11])
  mut1.index = which(dp_input[,1]==as.numeric(phasing[cl.phasing.pass.index[i.index],1]) & 
                    dp_input[,2]==as.numeric(phasing[cl.phasing.pass.index[i.index],2]))
  mut2.index = which(dp_input[,1]==as.numeric(phasing[cl.phasing.pass.index[i.index],1]) & 
                    dp_input[,2]==as.numeric(phasing[cl.phasing.pass.index[i.index],5]))

  if(length(mut1.index)==0 | length(mut2.index)==0){
    next
  }else{
    cl.phasing.pass[i.index,3:5] = unlist(dp_input[mut1.index,c(3,4,15)])
    cl.phasing.pass[i.index,7:9] = unlist(dp_input[mut2.index,c(3,4,15)])
    cl.phasing.pass[i.index,10:13] = unlist(phasing[cl.phasing.pass.index[i.index],8:11])
    cl.phasing.pass[i.index,14:15] = unlist(dp_input[mut1.index,6:7])
    cl.phasing.pass[i.index,16] = mut1.index
    cl.phasing.pass[i.index,17] = mut2.index

  }
}

cl.phasing.pass = cl.phasing.pass[which(!is.na(cl.phasing.pass[,16])),]
cl.phasing.pass = cl.phasing.pass[which(!is.na(cl.phasing.pass[,17])),]
cl.phasing.pass = as.data.frame(cl.phasing.pass)
save(cl.phasing.pass,file=file.path(outdir,paste0(samplename,"_cl_phasing_pass.RData")))

write.table(cl.phasing.pass,file=file.path(outdir,paste0(samplename,"_cl_phasing.txt")),sep="\t",row.names=F,quote=F)


