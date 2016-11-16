##############################################
# Copy number preparations
##############################################
#' This function takes a Battenberg subclones.txt file, pulls out the A copy number profile,
#' then it calculates the total CN and a weighted ploidy before it subsets the input data
#' by columns: "chr", "startpos", "endpos", "nMaj1_A", "nMin1_A", "frac1_A", 
#' "nMaj2_A", "nMin2_A", "frac2_A", "SDfrac_A". 
#' 
#' Then logic is applied for classifying copy number aberrations into
#' a series of classes. There are two main groups: Clonal (starting with c) and
#' Subclonal (starting with s).
#' 
#' A CNA can be:
#' * HD : Homozygous deletion
#' * LOH : Loss of heterozygosity
#' * Amp : Amplification
#' * Gain : Gain
#' * NoCNV : Normal copy number without aberrations
#' * Loss : Loss
#' 
#' The copy number segments, classification, weighted ploidy and samplename are saved
#' to the specified output file
#' 
#' @param samplename, A String containing the samplename
#' @param bb_subclones_file, A String containing the full path to a Battenberg subclones.txt file
#' @param gender, Specify either the string male or female
#' @param outfile, A String with the full path to where the output should be written
#' @author tjm
#' @export
collate_bb_subclones = function(samplename, bb_subclones_file, gender, outfile) {
  
  if(gender == 'male' | gender == 'Male') {
    is_male = T
  } else if(gender == 'female' | gender == 'Female') {
    is_male = F
  } else {
    stop("Unknown gender supplied, exit.")
  }
  
  allsegs = NULL
  cndata = read.table(bb_subclones_file, header=T, stringsAsFactors=F)
  cndata = cbind(samplename, cndata)
  # Remove the X chromosome when sample is male as Battenberg cannot infer copy number there
  if (is_male) {
    cndata = cndata[!cndata$chr=="X",]
  }
  
  # Remove all segments witout a call
  if (any(is.na(cndata$nMin1_A) | is.na(cndata$nMaj1_A))) {
    cndata = cndata[!(is.na(cndata$nMin1_A) | is.na(cndata$nMaj1_A)),]
  }
  
  # Obtain the total copy number for each segment
  totalcn = rep(0, dim(cndata)[1])
  for (k in 1:nrow(cndata)) {
    # Clonal copy number
    if (cndata$frac1_A[k] == 1) {
      totalcn[k] = cndata$nMaj1_A[k] + cndata$nMin1_A[k]
      # Subclonal copy number - which is represented as a mixture of two states
    } else if (cndata$frac1_A[k]!="Inf" & cndata$frac1_A[k]!="-Inf") {
      totalcn[k] = (cndata$nMaj1_A[k] + cndata$nMin1_A[k]) * cndata$frac1_A[k] + (cndata$nMaj2_A[k] + cndata$nMin2_A[k]) * cndata$frac2_A[k]
    }
    cndataweighted = totalcn * (cndata$endpos - cndata$startpos) / (sum(as.numeric(cndata$endpos-cndata$startpos)))
  }
  # Ploidy weighted by segment size
  ploidy = round(sum(cndataweighted)/2)*2
  
  # If there is no column with a tumour name, add it in temporarily
  if (! "Tumour_Name" %in% colnames(cndata)) {
    cndata = data.frame(Tumour_Name=samplename, cndata)
  }
  
  allsegs = data.frame(cndata[, c("Tumour_Name", "chr", "startpos", 
                                  "endpos", "nMaj1_A", "nMin1_A", 
                                  "frac1_A", "nMaj2_A", "nMin2_A", 
                                  "frac2_A", "SDfrac_A")], tumour_ploidy=ploidy)
  
  # Now classify all segments into a category
  tot = dim(allsegs)[1]
  # if you have clonal LOH, subclonal LOH is not counted!!!  Losses not counted
  allsegsa <- NULL
  CNA <- NULL
  for (i in 1:dim(allsegs)[1]) {
    # Clonal copy number
    if (is.na(allsegs$nMaj2_A[i])) {
      if (allsegs$nMin1_A[i] == 0 & allsegs$nMaj1_A[i] == 0) {
        allsegsa <- rbind(allsegsa, allsegs[i,c(1:7,11:12)])
        CNA <- c(CNA, "cHD")
      }
      if (xor(allsegs$nMin1_A[i] == 0, allsegs$nMaj1_A[i] == 0)) {
        allsegsa <- rbind(allsegsa, allsegs[i,c(1:7,11:12)])
        CNA <- c(CNA, "cLOH")
      }
      if ((allsegs$nMaj1_A[i] > allsegs$tumour_ploidy[i]/2) |
            (allsegs$nMin1_A[i] > allsegs$tumour_ploidy[i]/2)) {
        if ((allsegs$nMaj1_A[i] > allsegs$tumour_ploidy[i]*2) |
              (allsegs$nMin1_A[i] > allsegs$tumour_ploidy[i]*2)) {
          allsegsa <- rbind(allsegsa, allsegs[i,c(1:7,11:12)])
          CNA <- c(CNA, "cAmp")
        }
        else {
          allsegsa <- rbind(allsegsa, allsegs[i,c(1:7,11:12)])
          CNA <- c(CNA, "cGain")
        }
      }
      if ((allsegs$nMaj1_A[i] == allsegs$tumour_ploidy[i]/2) &
            (allsegs$nMin1_A[i] == allsegs$tumour_ploidy[i]/2)) {
        allsegsa <- rbind(allsegsa, allsegs[i,c(1:7,11:12)])
        CNA <- c(CNA, "NoCNV")
      }
      if ((allsegs$nMaj1_A[i] < allsegs$tumour_ploidy[i]/2) |
            (allsegs$nMin1_A[i] < allsegs$tumour_ploidy[i]/2)) {
        allsegsa <- rbind(allsegsa, allsegs[i,c(1:7,11:12)])
        CNA <- c(CNA, "cLoss")
      }
    }
    
    # Subclonal copy number
    if (!is.na(allsegs$nMaj2_A[i])) {
      if (allsegs$nMin1_A[i] == 0 & allsegs$nMaj1_A[i] == 0) {
        CNA <- c(CNA, "sHD")
        allsegsa <- rbind(allsegsa, allsegs[i,c(1:7,11:12)])
      }
      if ((allsegs$nMin1_A[i] == 0 | allsegs$nMaj1_A[i] == 0) &
            xor(allsegs$nMin2_A[i] == 0, allsegs$nMaj2_A[i] == 0)) {
        CNA <- c(CNA, "cLOH")
        tmp <- allsegs[i,c(1:4,8:12)]
        names(tmp) <- names(allsegs[i,c(1:7,11:12)])
        allsegsa <- rbind(allsegsa, tmp[,names(allsegs[i,c(1:7,11:12)])])
      }
      else if (xor(allsegs$nMin1_A[i] == 0, allsegs$nMaj1_A[i] == 0)) {
        CNA <- c(CNA, "sLOH")
        allsegsa <- rbind(allsegsa, allsegs[i,c(1:7,11:12)])
      }
      if ((allsegs$nMaj1_A[i] > allsegs$tumour_ploidy[i]/2 |
             allsegs$nMin1_A[i] > allsegs$tumour_ploidy[i]/2) &
            (allsegs$nMaj2_A[i] > allsegs$tumour_ploidy[i]/2 |
               allsegs$nMin2_A[i] > allsegs$tumour_ploidy[i]/2)) {
        if ((allsegs$nMaj1_A[i] > allsegs$tumour_ploidy[i]*2 |
               allsegs$nMin1_A[i] > allsegs$tumour_ploidy[i]*2) &
              (allsegs$nMaj2_A[i] > allsegs$tumour_ploidy[i]*2 |
                 allsegs$nMin2_A[i] > allsegs$tumour_ploidy[i]*2)) {
          CNA <- c(CNA, "cAmp")
          allsegsa <- rbind(allsegsa, allsegs[i,c(1:7,11:12)])
        }
        else {
          CNA <- c(CNA, "cGain")
          allsegsa <- rbind(allsegsa, allsegs[i,c(1:7,11:12)])
        }
      }
      if ((allsegs$nMin2_A[i] > allsegs$tumour_ploidy[i]/2 & allsegs$nMin1_A[i] != allsegs$nMin2_A[i])|
            (allsegs$nMaj2_A[i] > allsegs$tumour_ploidy[i]/2 & allsegs$nMaj1_A[i] != allsegs$nMaj2_A[i])) {
        if ((allsegs$nMin2_A[i] > allsegs$tumour_ploidy[i]*2 & allsegs$nMin1_A[i] != allsegs$nMin2_A[i]) |
              (allsegs$nMaj2_A[i] > allsegs$tumour_ploidy[i]*2 & allsegs$nMaj1_A[i] != allsegs$nMaj2_A[i])) {
          CNA <- c(CNA, "sAmp")
          tmp <- allsegs[i,c(1:4,8:12)]
          names(tmp) <- names(allsegs[i,c(1:7,11:12)])
          allsegsa <- rbind(allsegsa, tmp[,names(allsegs[i,c(1:7,11:12)])])
        }
        else {
          CNA <- c(CNA, "sGain")
          tmp <- allsegs[i,c(1:4,8:12)]
          names(tmp) <- names(allsegs[i,c(1:7,11:12)])
          allsegsa <- rbind(allsegsa, tmp[,names(allsegs[i,c(1:7,11:12)])])
        }
      }
      if ((allsegs$nMaj1_A[i] < allsegs$tumour_ploidy[i]/2 |
             allsegs$nMin1_A[i] < allsegs$tumour_ploidy[i]/2) &
            (allsegs$nMaj2_A[i] < allsegs$tumour_ploidy[i]/2 |
               allsegs$nMin2_A[i] < allsegs$tumour_ploidy[i]/2)) {
        CNA <- c(CNA, "cLoss")
        tmp <- allsegs[i,c(1:4,8:12)]
        names(tmp) <- names(allsegs[i,c(1:7,11:12)])
        allsegsa <- rbind(allsegsa, tmp[,names(allsegs[i,c(1:7,11:12)])])
      }
      if (((allsegs$nMaj1_A[i] < allsegs$tumour_ploidy[i]/2 & allsegs$nMaj1_A[i] != allsegs$nMaj2_A[i]) |
             (allsegs$nMin1_A[i] < allsegs$tumour_ploidy[i]/2 & allsegs$nMin1_A[i] != allsegs$nMin2_A[i]))) {
        allsegsa <- rbind(allsegsa, allsegs[i,c(1:7,11:12)])
        CNA <- c(CNA, "sLoss")      
      }
    }
    if(i %% 100 ==0){
      print(paste(i,"/",tot))
    }
  }  
  
  allsegsa <- cbind(allsegsa, CNA)
  allsegsa$frac1_A[allsegsa$CNA == "cGain"|allsegsa$CNA == "cAmp"|allsegsa$CNA == "cLOH"|allsegsa$CNA == "cHD"|allsegsa$CNA == "cLoss"] = 1
  
  # Switch the columns
  allsegsa = data.frame(allsegsa[,-1], Tumour_Name=samplename)
  
  write.table(allsegsa, file=outfile, sep="\t", quote=F, row.names=F)
}
#' This function extends existing categories by a list of sizes
#' @param infile string that points to a collated Battenberg subclones file
#' @param outfile string where the output will be written.
#' @param size_categories named vector with size options to consider for each CNA category
#' @author sd11
extend_bb_cn_categories = function(infile, outfile, size_categories=c("3"=10^3, "4"=10^4, "5"=10^5, "6"=10^6, "7"=10^7, "8"=10^8, "9"=10^9, "10"=10^10)) {
  dat = read.table(infile, header=T, stringsAsFactors=F)
  
  # Assign a segment to the first category that is larger than it
  dat$category = sapply(1:nrow(dat), 
                        function(i, dat, size_categories) { 
                          seg_len = dat$endpos[i] - dat$startpos[i]
                          cat_index = which(!(size_categories < seg_len))[1]
                          return(paste(dat$CNA[i], "_", names(cat_index), sep=""))
                        }, dat=dat, size_categories=size_categories)
  
  write.table(dat, file=outfile, sep="\t", row.names=F, quote=F)
}

#' Count how often each of the categories is present in this sample
#' @param infile string pointing to a collated Battenberg subclones file with extended CNA categories
#' @param outfile string pointing to where the catalog should be written
#' @param size_categories named vector with the size options to consider. This should be the same vector given to extend_bb_cn_categories
#' @author sd11
get_bb_cn_catalog = function(infile, outfile, size_categories=c("3"=10^3, "4"=10^4, "5"=10^5, "6"=10^6, "7"=10^7, "8"=10^8, "9"=10^9, "10"=10^10)) {
  dat = read.table(infile, header=T, stringsAsFactors=F)
  # Match the CN types with the size categories to create all catalog entries
  CN_types = c("NoCNV", "cAmp", "cGain", "cLOH", "cLoss", "sAmp", "sGain", "sLOH", "sLoss")
  CN_classes = sapply(1:length(size_categories), 
                      function(i, CN_classes, size_categories) { paste(CN_types, names(size_categories)[i], sep="_") }, 
                      CN_classes=CN_classes, size_categories=size_categories)
  # Count how often each category is available
  catalog = do.call(rbind, lapply(1:length(CN_classes), 
                                  function(i, CN_classes, dat) { data.frame(category=CN_classes[i], count=sum(dat$category==CN_classes[i])) }, 
                                  CN_classes=CN_classes, dat=dat))
  
  write.table(catalog, file=outfile, quote=F, sep="\t", row.names=F)
}
