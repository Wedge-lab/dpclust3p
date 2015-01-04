# args=commandArgs(TRUE)
# fofn = toString(args[1]) # List of input files
# inputdir = toString(args[2]) # Dir where input files reside
# outfile = toString(args[3]) # Full path to where output should be written
# haveHeader = as.logical(args[4]) # TRUE or FALSE depending on whether the fofn files have a header

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