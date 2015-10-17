#' Function to convert mutation burdens into mutation copy number
#' @param burden A vector containing mutation burdens
#' @param totalCopyNumber A vector with total tumour copynumber
#' @param cellularity Float with the cellularity of this sample
#' @param normalCopyNumber A vector with the total copy number of the normal
#' @return A vector with mutation copy number
#' @author dw9
#' @export
mutationBurdenToMutationCopyNumber = function(burden, totalCopyNumber, cellularity, normalCopyNumber=rep(2, length(burden))) {
  mutCopyNumber = burden/cellularity*(cellularity*totalCopyNumber+normalCopyNumber*(1-cellularity))
  mutCopyNumber[is.nan(mutCopyNumber)] = 0
  return(mutCopyNumber)
}

#' Function to convert mutation copy number to mutation burden
#' @param copyNumber A vector containing mutation copy number
#' @param totalCopyNumber A vector with total tumour copynumber
#' @param cellularity Float with the cellularity of this sample
#' @param normalCopyNumber A vector with the total copy number of the normal
#' @return A vector with mutation burdens
#' @author dw9
#' @export
mutationCopyNumberToMutationBurden = function(copyNumber, totalCopyNumber, cellularity, normalCopyNumber=rep(2, length(copyNumber))){
  burden = copyNumber*cellularity/(cellularity*totalCopyNumber+normalCopyNumber*(1-cellularity))
  burden[is.nan(burden) | (burden<0.000001)] = 0.000001
  burden[burden>0.999999] = 0.999999
  return(burden)	
}