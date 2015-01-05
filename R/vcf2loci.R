#library(VariantAnnotation)
#args=commandArgs(TRUE)
#vcf_file = toString(args[1]) # Full path to vcf file
#fai_file= toString(args[2]) # Full path to fai file
#ign_file = toString(args[3]) # Full path to ignore file
#outfile = toString(args[4]) # Full path to output file

parseFai = function(fai_file) {
                fai = read.table(fai_file, header=F, stringsAsFactor=F)
        colnames(fai) = c("chromosome", "length", "offset", "fasta_line_length", "line_blen")
                return(fai)
}

parseIgnore = function(ignore_file) {
                ign = read.table(ignore_file, header=F, stringsAsFactor=F)
        colnames(ign) = c("chromosome")
                return(ign)
}

vcf2loci = function(vcf_file, fai_file, ign_file, outfile) {
#fai = parseFai("/lustre/scratch110/sanger/sd11/Documents/GenomeFiles/refs_icgc_pancan/genome.fa.fai")
#ign = parseIgnore("/lustre/scratch110/sanger/sd11/Documents/GenomeFiles/battenberg_ignore/ignore.txt")
fai = parseFai(fai_file)
ign = parseIgnore(ign_file)
allowed_chroms = which(!(fai$chromosome %in% ign$chromosome))

# Setup that uses VariantAnnotation but takes too long
# Make sure only the allowed chromosomes are read in
#param = ScanVcfParam(which=GRanges(fai$chromosome[allowed_chroms], IRanges(rep(1, length(allowed_chroms)), fai$length[allowed_chroms])))
#vcf = readVcf("PD7404a.cave.annot.vcf.gz", "hg19", param=param)
# Parse the vcf into the loci format and write to disk
#vcf.loci = as.data.frame(head(rowData(vcf))[,c("REF","ALT")])
#vcf.loci = data.frame(seqnames(vcf), start(vcf), fixed(vcf)[, c("REF","ALT")])
#vcf.loci = cbind(as.data.frame(head(ranges(vcf))), head(fixed(vcf))[, c("REF","ALT")])
#vcf.loci = cbind(as.data.frame(ranges(vcf)), fixed(vcf)[, c("REF","ALT")])

vcf.cols = ncol(read.delim(vcf_file, comment.char="#", header=F, stringsAsFactor=F, nrows=1))
vcf.cols.default = 10 # vcf file standard contains 10 columns
vcf.colClasses = c(NA, NA, "NULL", NA, NA, rep("NULL", 5+(vcf.cols-vcf.cols.default)))
vcf.loci = read.delim(vcf_file, comment.char="#", header=F, stringsAsFactor=F, colClasses=vcf.colClasses)
colnames(vcf.loci) = c("chromosome", "pos", "ref","alt")
vcf.loci.sel = subset(vcf.loci, chromosome %in% fai$chromosome[allowed_chroms])
write.table(vcf.loci.sel, col.names=F, quote=F, row.names=F, file=outfile, sep="\t")
}
