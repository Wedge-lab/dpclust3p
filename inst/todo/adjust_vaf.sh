#!/bin/bash
infile=$1
samplename=`basename ${infile} | cut -f 1 -d _`
cell=`grep ${samplename} ../dirichlet/conscna.txt  | cut -f 4`
Rscript adjust_vaf.R ${infile} ${cell}
