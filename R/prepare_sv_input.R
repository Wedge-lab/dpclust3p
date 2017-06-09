# suppressMessages(library(GenomicRanges))


prepare_sv_dpinput = function(svclone_file) {
  dat = read.table(svclone_file, header=T, stringsAsFactors=F, sep="\t")
  # bb = read.table("../../../copynumberConsensus/preliminary_consensus/preliminary_consensus_battenberg/GBM-US/1e27cc8a-5394-4958-9af6-5ece1fe24516/1e27cc8a-5394-4958-9af6-5ece1fe24516_subclones.txt", header=T, stringsAsFactors=F)
  
  mutCount = dat$adjusted_depth*dat$adjusted_vaf
  WTCount = dat$adjusted_depth-mutCount
  
  sv_chrom_pos = data.frame()
  for (i in 1:nrow(dat)) {
    # Preferred copy number
    if (dat$preferred_side[i]==0) {
      sv_chrom_pos = rbind(sv_chrom_pos, data.frame(chrom=dat$chr1[i], pos=dat$pos1[i]))
    } else {
      sv_chrom_pos = rbind(sv_chrom_pos, data.frame(chrom=dat$chr2[i], pos=dat$pos2[i]))
    }
  }
  
  all_dat = GRanges(sv_chrom_pos$chrom, IRanges(sv_chrom_pos$pos, sv_chrom_pos$pos), mut.count=mutCount, WT.count=WTCount)
  return(all_dat)
}

runGetDirichletProcessInfo_sv = function(svclone_file, cellularity_file, subclone_file, gender, output_file) {
  if(gender == 'male' | gender == 'Male') {
    isMale = T
  } else if(gender == 'female' | gender == 'Female') {
    isMale = F
  } else {
    stop("Unknown gender supplied, exit.")
  }
  info_counts = prepare_sv_dpinput(svclone_file)
  cellularity = dpclust3p:::GetCellularity(cellularity_file)
  
  print(info_counts)
  print(cellularity)
  
  dpclust3p:::GetDirichletProcessInfo(output_file, cellularity, info_counts, subclone_file, is.male=isMale, SNP.phase.file="NA", mut.phase.file="NA")
}
  

# args = commandArgs(T)
# svclone_file = args[1]
# cellularity_file = args[2]
# subclone_file = args[3]
# gender = args[4]
# output_file = args[5]
# 
#   
# # cellularity_file = "../../../copynumberConsensus/preliminary_consensus/preliminary_consensus_battenberg/GBM-US/1e27cc8a-5394-4958-9af6-5ece1fe24516/1e27cc8a-5394-4958-9af6-5ece1fe24516_rho_and_psi.txt"
# # subclone_file = "../../../copynumberConsensus/preliminary_consensus/preliminary_consensus_battenberg/GBM-US/1e27cc8a-5394-4958-9af6-5ece1fe24516/1e27cc8a-5394-4958-9af6-5ece1fe24516_subclones.txt"
# # gender = "female"
# # output_file = "1e27cc8a-5394-4958-9af6-5ece1fe24516_test.txt"
# runGetDirichletProcessInfo_sv(svclone_file, cellularity_file, subclone_file, gender, output_file)