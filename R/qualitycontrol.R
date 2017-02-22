##############################################
# QC
##############################################
#' Create a PNG file
#' @param p The figure
#' @param filename Where to save the image
#' @param width Width
#' @param height Height
#' @export
createPng <- function(p, filename, width, height) {
  png(filename=filename, width=width, height=height)
  print(p)
  dev.off()
}

#' Convenience function used for creating figures
#' @noRd
meltFacetPlotData <- function(data, subsamplenames) {
  d = as.data.frame(data)
  colnames(d) = subsamplenames
  data.m = reshape2::melt(d)
  return(data.m)
}

#' Template to create histograms
#' @noRd
createHistFacetPlot <- function(data, title, xlab, ylab, binwidth) {
  p = ggplot2::ggplot(data) + ggplot2::aes(x=value, y=..count..) + ggplot2::geom_histogram(binwidth=binwidth, colour="black", fill="gray") + ggplot2::facet_grid(variable ~ .)
  p = p + ggplot2::theme_bw(base_size=35) + ggplot2::ggtitle(title) + ggplot2::xlab(xlab) + ggplot2::ylab(ylab)
  return(p)
}

#' Template to create box plots
#' @noRd
createBoxFacetPlot <- function(data, title, xlab, ylab) {
  p = ggplot2::ggplot(data) + ggplot2::aes(x=variable, y=value) + ggplot2::geom_boxplot() + ggplot2::facet_grid(subsample ~ .)
  p = p + ggplot2::theme_bw() + ggplot2::ggtitle(title) + ggplot2::xlab(xlab) + ggplot2::ylab(ylab)
  return(p)
}

#' Plot mutation copy number histogram
#' 
#' @param data A DPClust input table
#' @param samplename Name of the sample
#' @param outdir Directory where the figure is to be stored
#' @author sd11
#' @export
plot_mcn_hist = function(data, samplename, outdir) {
  p = createHistFacetPlot(meltFacetPlotData(data$mutation.copy.number, samplename), paste(samplename, "mutation.copy.number"), "mutation.copy.number", "Count", binwidth=0.1)
  p = p + ggplot2::xlim(0,5)
  createPng(p, file.path(outpath, paste(samplename, "_mutation.copy.number.png", sep="")), width=1500, height=500)
}

#' Plot cancer cell fraction histogram
#' 
#' @param data A DPClust input table
#' @param samplename Name of the sample
#' @param outdir Directory where the figure is to be stored
#' @author sd11
#' @export
plot_ccf_hist = function(data, samplename, outdir) {
  p = createHistFacetPlot(meltFacetPlotData(data$subclonal.fraction, samplename), paste(samplename, "Fraction Tumour Cells"), "Fraction Tumour Cells", "Count", binwidth=0.05)
  p = p + ggplot2::xlim(0,1.5)
  createPng(p, file.path(outpath, paste(samplename, "_fractionOfCells.png", sep="")), width=1500, height=500)
}

#' Plot allele frequency histogram
#' 
#' @param data A DPClust input table
#' @param samplename Name of the sample
#' @param outdir Directory where the figure is to be stored
#' @author sd11
#' @export
plot_vaf_hist = function(data, samplename, outdir) {
  p = createHistFacetPlot(meltFacetPlotData(data$mutCount/(data$mutCount+data$WTCount), samplename), paste(samplename, "alleleFrequency"), "Allele Frequency", "Count", binwidth=0.01)
  createPng(p, file.path(outpath, paste(samplename, "_alleleFrequency.png", sep="")), width=1500, height=500)
}
