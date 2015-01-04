# args=commandArgs(TRUE)
# infile = toString(args[1]) # Where data is stored
# prefix = toString(args[2]) # _muts_linkedMuts_segmented_
# postfix = toString(args[3]) # .txt
# outdir = toString(args[4]) # Where to write the output
# chrom_file = toString(args[5]) # A file containing the chromosome name (first column) and it's integer number (second column)

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