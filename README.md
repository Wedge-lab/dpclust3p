
# DPClust pre-processing

This R package contains various functions to produce input data for [DPClust](https://github.com/Wedge-Oxford/dpclust) using SNV variant calls and copy number data from [Battenberg](https://github.com/Wedge-Oxford/battenberg). Most importantly, it contains the `runGetDirichletProcessInfo` function that produces the input data for SNV based clustering.

## Installation instructions
dpclust3p is an R package and can be installed with the commands right below. It also requires the [alleleCounter](https://github.com/cancerit/allelecount) tool to be in `$PATH`.
```
source("http://bioconductor.org/biocLite.R"); biocLite(c("optparse","VariantAnnotation","GenomicRanges","Rsamtools","ggplot2","IRanges","S4Vectors","reshape2"))'
devtools::install_github("Wedge-Oxford/dpclust3p")
```

## Running pre-processing

The typical usage is to create the DPClust input data. See `inst/example` for a few example pipelines. A pipeline typically consists of three steps:
 * Transform loci from a VCF file into a loci file
 * Obtain allele counts for all mutations, either by invoking alleleCount or by dumping counts from the VCF file
 * Convert allele counts and copy number information into DPClust input

The R package contains many functions from which one can build their own pipeline

| File | Description |
|---|---|
| preprocessing.R | Main preprocessing functions to create DPClust input, perform mutation phasing, filter by mutational signature |
| allelecount.R | Functions to count alleles in a BAM file, or dump counts from a range of VCF formats |
| kataegis.R | Functions to identify kataegis events (requires fastPCF.R) |
| copynumber.R | Various functions related to copy number |
| qualitycontrol.R | Create plots that can be used for QCing |
| interconvertMutationBurdens.R | Basic functions for data transformations |
| util.R | Various utility functions |

## Docker

This package has been Dockerised, build as follows:
```
docker build -t dpclust3p:1.0.8 .
```

