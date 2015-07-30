#' THIS FUNCTION IS DEPRECATED - SHOULD BE MERGED WITH THE ONE IN DPCLUST
#' 
#' load.data
#'
#' Load data files for a specific set of samples belonging to the same individual/experiment.
#'
#' @param datpath Path to where the input datafiles can be found, type string.
#' @param samplename A unique identifier for this particular individual/experiment, type string.
#' @param list_of_data_files A list of file names that are available in datpath, list of strings.
#' @param cellularity Vector containing the cellularity estimate for each sample in this set, vector of floats.
#' @param Chromosome The name of the column containing the chromosome information in the input, type string.
#' @param position The name of the column containing the position information in the input, type string.
#' @param WT.count The name of the column containing the wild type allele count information in the input, type string.
#' @param mut.count The name of the column containing the mutation allele count information in the input, type string.
#' @param subclonal.CN The name of the column containing the total subclonal copynumber information in the input, type string.
#' @param no.chrs.bearing.mut The name of the column containing the number of chromsomes bearing the mutation information in the input, type string.
#' @param data_file_suffix A suffix to be appended to each item in list_of_data_files to complete the filenames (for example dirichletInput.R), type string.
#' @return A named list of matrices and vectors with the following items: chromosome, position, WTCount, mutCount, totalCopyNumber, copyNumberAdjustment, non.deleted.muts, kappa, mutation.copy.number, subclonal.fraction, removed_indices, chromosome.not.filtered, mut.position.not.filtered
#' @author Stefan Dentro
load.data <- function(datpath, samplename, list_of_data_files, cellularity, Chromosome, position, WT.count, mut.count, subclonal.CN, no.chrs.bearing.mut, mutation.copy.number, subclonal.fraction, data_file_suffix) {
  warning("The load.data function is deprecated and should be merged with the one in DPCLUST")
  data=list()
  for(s in 1:length(list_of_data_files)){
    data[[s]] = read.table(paste(datpath,samplename,list_of_data_files[s],data_file_suffix,sep=""),header=T,sep="\t",stringsAsFactors=F)
  }

  no.subsamples = length(list_of_data_files)
  no.muts = nrow(data[[1]])
  
  chromosome = matrix(0,no.muts,no.subsamples)
  mut.position = matrix(0,no.muts,no.subsamples)
  WTCount = matrix(0,no.muts,no.subsamples)
  mutCount = matrix(0,no.muts,no.subsamples)
  totalCopyNumber = matrix(0,no.muts,no.subsamples)
  copyNumberAdjustment = matrix(0,no.muts,no.subsamples)
  non.deleted.muts = vector(mode="logical",length=nrow(data[[1]]))
  mutationCopyNumber = matrix(NA,no.muts,no.subsamples)
  subclonalFraction = matrix(NA,no.muts,no.subsamples)
  for(s in 1:length(list_of_data_files)){
    chromosome[,s] = data[[s]][,Chromosome]
    mut.position[,s] = as.numeric(data[[s]][,position])
    WTCount[,s] = as.numeric(data[[s]][,WT.count])
    mutCount[,s] = as.numeric(data[[s]][,mut.count])
    totalCopyNumber[,s] = as.numeric(data[[s]][,subclonal.CN])
    copyNumberAdjustment[,s] = as.numeric(data[[s]][,no.chrs.bearing.mut])
    non.deleted.muts[data[[s]][,no.chrs.bearing.mut]>0]=T
    mutationCopyNumber[,s] = as.numeric(data[[s]][,mutation.copy.number])
    subclonalFraction[,s] = as.numeric(data[[s]][,subclonal.fraction])
  }
  
  kappa = matrix(1,no.muts,no.subsamples)
  for(i in 1:length(list_of_data_files)){
    # multiply by no.chrs.bearing.mut, so that kappa is the fraction of reads required for fully clonal mutations, rather than mutation at MCN = 1
    kappa[,i] = mutationCopyNumberToMutationBurden(1,data[[i]][,subclonal.CN],cellularity[i]) * data[[i]][,no.chrs.bearing.mut]
  }

  # Remove those mutations that have missing values
  not.there.wt = apply(WTCount, 1, function(x) { sum(is.na(x))>0 })
  not.there.mut = apply(mutCount, 1, function(x) { sum(is.na(x))>0 })
  not.there.cn = apply(totalCopyNumber, 1, function(x) { sum(is.na(x))>0 })
  not.there.cna = apply(copyNumberAdjustment, 1, function(x) { sum(is.na(x))>0 })
  not.there.kappa = apply(kappa, 1, function(x) { sum(is.na(x))>0 })
  # Remove those mutations that have no coverage. These cause for trouble lateron.
  not.coverage = apply(WTCount+mutCount, 1, function(x) { any(x==0) })
  not.cna = apply(copyNumberAdjustment, 1, function(x) { any(x==0) })
  not.on.supported.chrom = apply(chromosome, 1, function(x) { ! any(x %in% as.character(1:22)) })

  cov.mean = mean(colMeans(WTCount+mutCount))
  cov.std = mean(apply((WTCount+mutCount), 2, sd))
  too.high.coverage = apply(WTCount+mutCount, 1, function(x) { any(x > cov.mean+6*cov.std)})

  print(paste("Removed", sum(not.there.wt),"with missing WTCount", sep=" "))
  print(paste("Removed", sum(not.there.mut),"with missing mutCount", sep=" "))
  print(paste("Removed", sum(not.there.cn),"with missing totalCopyNumber", sep=" "))
  print(paste("Removed", sum(not.there.cna),"with missing copyNumberAdjustment", sep=" "))
  print(paste("Removed", sum(not.there.kappa),"with missing kappa", sep=" "))
  print(paste("Removed", sum(not.coverage),"with no coverage", sep=" "))
  print(paste("Removed", sum(not.cna),"with zero copyNumberAdjustment", sep=" "))
  print(paste("Removed", sum(not.on.supported.chrom), "on not supported genomic regions", sep=" "))
  print(paste("Removed", sum(too.high.coverage), "mutations with coverage over",cov.mean+6*cov.std, sep=" "))

  select = !(not.there.wt | not.there.mut | not.there.cn | not.there.cna | not.there.kappa | not.coverage | not.cna | not.on.supported.chrom | too.high.coverage)
  
  # Keep indices of removed mutations to 'spike in' lateron when constructing the output
  removed_indices = which(!select)
  chromosome.not.filtered = chromosome
  mut.position.not.filtered = mut.position
  
  # Remove mutations that have been flagged for various reasons
  chromosome = as.matrix(chromosome[select,])
  mut.position = as.matrix(mut.position[select,])
  WTCount = as.matrix(WTCount[select,])
  mutCount = as.matrix(mutCount[select,])
  totalCopyNumber = as.matrix(totalCopyNumber[select,])
  copyNumberAdjustment = as.matrix(copyNumberAdjustment[select,])
  non.deleted.muts = as.matrix(non.deleted.muts[select])
  kappa = as.matrix(kappa[select,])
  mutationCopyNumber = as.matrix(mutationCopyNumber[select,])
  subclonalFraction = as.matrix(subclonalFraction[select,])
  print(paste("Removed",no.muts-nrow(WTCount), "mutations with missing data"))
  
  return(list(chromosome=chromosome, position=mut.position, WTCount=WTCount, mutCount=mutCount, 
              totalCopyNumber=totalCopyNumber, copyNumberAdjustment=copyNumberAdjustment, 
              non.deleted.muts=non.deleted.muts, kappa=kappa, mutation.copy.number=mutationCopyNumber, 
              subclonal.fraction=subclonalFraction, removed_indices=removed_indices,
              chromosome.not.filtered=chromosome.not.filtered, mut.position.not.filtered=mut.position.not.filtered))
}
