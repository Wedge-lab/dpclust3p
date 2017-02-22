#' Concatenate split files
#' 
#' Convenience function to concatenate a series of files specified in a file of file names.
#' This function assumes all files have the same layout.
#' @param fofn A file of file names to be concatenated
#' @param inputdir Full path to where the input files are stored
#' @param outfile Full path to where the output should be written
#' @param haveHeader Boolean that specifies whether the input files have a header
#' @author sd11
#' @export
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

#' Split a file per chromosome
#'
#' Convenience function to split an input file per chromosome. All it requires is that
#' the infile has as first column chromosome specification. The output files will be named
#' outdir/prefixCHROMNUMBERpostfix
#' @param infile The file to be split
#' @param prefix Prefix of the output file
#' @param postfix Postfix of the output file
#' @param outdir Directory where the output files are to be written
#' @param chrom_file A simple list of chromosomes to be considered
#' @author sd11
#' @export
split_by_chrom = function(infile, prefix, postfix, outdir, chrom_file) {
  outdir = paste(outdir, "/", sep="")
  
  # Check if there are lines in the file, otherwise it will crash this script
  if (file.info(infile)$size == 0) {
    print("No lines in loci file")
    q(save="no")
  }
  
  loci = read.delim(infile, stringsAsFactors=F, header=F)
  
  chroms = read.delim(chrom_file, stringsAsFactors=F, header=F)
  
  for (i in 1:nrow(chroms)) {
    selection = loci[loci[,1]==chroms[i,1],]
    chrom_id = chroms[i,2]
    write.table(selection, file=paste(outdir, prefix, chrom_id, postfix, sep=""), quote=F, row.names=F, col.names=F, sep="\t")
  }
}

#' Calculate power to call subclones. 
#' 
#' This function calculates the 
#' average number of reads per clonal chromosome copy.
#' @param purity The sample purity
#' @param ploidy The tumour ploidy
#' @param coverage The tumour sample coverage
#' @return The average reads per chromosome copy, a.k.a. power
#' @author sd11
calc_power = function(purity, ploidy, coverage) {
  return(round((purity) / (purity*ploidy + (1-purity)*2) * coverage, 3))
}