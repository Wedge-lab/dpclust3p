############################################
# Copy number convert scripts
############################################

#' Transform ASCAT output into Battenberg
#' 
#' This function takes the ascat segments, acf and ploidy files to create
#' a minimum Battenberg subclones and rho_and_psi file for use in pre-processing.
#' These files only contain the essentials required by this package.
#' @param outfile.prefix String prefix for the output files. Internally _subclones.txt and _rho_and_psi.txt will be added.
#' @param segments.file String pointing to a segments.txt ASCAT output file
#' @param acf.file String pointing to a acf.txt ASCAT output file
#' @param ploidy.file String pointing to a ploidy.txt ASCAT output file
#' @author sd11
#' @export
ascatToBattenberg = function(outfile.prefix, segments.file, acf.file, ploidy.file) {
  # Construct a minimum Battenberg file with the copy number segments
  d = read.table(segments.file, header=T, stringsAsFactors=F)
  subclones = d[,2:6]
  colnames(subclones)[4] = "nMaj1_A"
  colnames(subclones)[5] = "nMin1_A"
  subclones$frac1_A = 1
  subclones$nMaj2_A = NA
  subclones$nMin2_A = NA
  subclones$frac2_A = NA
  subclones$SDfrac_A = NA
  write.table(subclones, paste(outfile.prefix, "_subclones.txt", sep=""), quote=F, sep="\t")
  
  # Now construct a minimum rho/psi file
  cellularity = read.table(acf.file, header=T)[1,1]
  ploidy = read.table(ploidy.file, header=T)[1,1]
  
  purity_ploidy = array(NA, c(3,5))
  colnames(purity_ploidy) = c("rho", "psi", "ploidy", "distance", "is.best")
  rownames(purity_ploidy) = c("ASCAT", "FRAC_GENOME", "REF_SEG")
  purity_ploidy["FRAC_GENOME", "rho"] = cellularity
  purity_ploidy["FRAC_GENOME", "psi"] = ploidy
  write.table(purity_ploidy, paste(outfile.prefix, "_rho_and_psi.txt", sep=""), quote=F, sep="\t")
}

#' Transform ASCAT NGS output into Battenberg
#' 
#' This function takes the ascat copynumber.caveman.csv and samplestatistics files to create
#' a minimum Battenberg subclones and rho_and_psi file for use in pre-processing.
#' These files only contain the essentials required by this package.
#' @param outfile.prefix String prefix for the output files. Internally _subclones.txt and _rho_and_psi.txt will be added.
#' @param copynumber.caveman.file String pointing to an ASCAT NGS copynumber.caveman.csv file
#' @param samplestatistics.file String point to an ASCAT NGS samplestatistics.csv file
#' @author sd11
#' @export
ascatNgsToBattenberg = function(outfile.prefix, copynumber.caveman.file, samplestatistics.file) {
  # Construct a minimum Battenberg file with the copy number segments
  d = read.table(copynumber.caveman.file, sep=",", header=F, stringsAsFactors=F)
  colnames(d) = c("count", "chr", "startpos", "endpos", "normal_total", "normal_minor", "tumour_total", "tumour_minor")
  subclones = d[,2:4]
  subclones$nMaj1_A = d$tumour_total-d$tumour_minor
  subclones$nMin1_A = d$tumour_minor
  subclones$frac1_A = 1
  subclones$nMaj2_A = NA
  subclones$nMin2_A = NA
  subclones$frac2_A = NA
  write.table(subclones, paste(outfile.prefix, "_subclones.txt", sep=""), quote=F, sep="\t")
  
  # Now construct a minimum rho/psi file
  samplestats = read.table(samplestatistics.file, header=F, stringsAsFactors=F)
  
  purity_ploidy = array(NA, c(3,5))
  colnames(purity_ploidy) = c("rho", "psi", "ploidy", "distance", "is.best")
  rownames(purity_ploidy) = c("ASCAT", "FRAC_GENOME", "REF_SEG")
  purity_ploidy["FRAC_GENOME", "rho"] = samplestats[samplestats$V1=="rho",2]
  purity_ploidy["FRAC_GENOME", "psi"] = samplestats[samplestats$V1=="Ploidy",2]
  write.table(purity_ploidy, paste(outfile.prefix, "_rho_and_psi.txt", sep=""), quote=F, sep="\t")
}

#' Check if the Y chromosome is present and the donor is male. If this statement is TRUE we need to 
#' there is no estimate of the Y chromosome, which BB currently provides as two copies of X
#' @param subclone.data A read-in Battenberg subclones.txt file as a data.frame
#' @return The input data.frame with a Y chromosome entry added and possible X chromosome adjustment
#' @author sd11
#' @noRd
addYchromToBattenberg = function(subclone.data) {
  # Take the largest segment on X chromosome
  xlargest = which.max(subclone.data[subclone.data$chr=="X", ]$endpos - subclone.data[subclone.data$chr=="X", ]$startpos)
  xlargest = subclone.data[subclone.data$chr=="X", ][xlargest,]
  
  # Get the total availability of X
  xnmaj = xlargest$nMaj1_A * xlargest$frac1_A + ifelse(is.na(xlargest$frac2_A), 0, xlargest$nMaj2_A * xlargest$frac2_A)
  xnmin = xlargest$nMin1_A * xlargest$frac1_A + ifelse(is.na(xlargest$frac2_A), 0, xlargest$nMin2_A * xlargest$frac2_A)
  
  # If there are no two copies we assume that Y was lost and make no changes
  if (xnmaj + xnmin >=2 & xnmin >= 1) {
    # We have at least 1 copy of either allele. One must therefore be the Y chromosome
    
    # Remove a copy of the X chromosome
    ynMaj_1 = subclone.data[subclone.data$chr=="X", c("nMin1_A")]
    yfrac_1 = subclone.data[subclone.data$chr=="X", c("frac1_A")]
    subclone.data[subclone.data$chr=="X", c("nMin1_A")] = 0
    
    # if there was subclonal copy number we need to remove/adapt that too
    if (!is.na(xlargest$frac2_A)) {
      ynMaj_2 = subclone.data[subclone.data$chr=="X", c("nMin2_A")]
      yfrac_2 = subclone.data[subclone.data$chr=="X", c("frac2_A")]
      subclone.data[subclone.data$chr=="X", c("nMin2_A")] = 0
    } else {
      ynMaj_2 = NA
      yfrac_2 = NA
    }
  } else {
    # No 2 copies, we assume that Y has been lost
    ynMaj_1 = 0
    yfrac_1 = 1
    ynMaj_2 = NA
    yfrac_2 = NA
  }
  
  # Add Y
  subclone.data = rbind(subclone.data, 
                        data.frame(chr="Y", startpos=0, endpos=59373566, BAF=NA, pval=NA, LogR=NA, ntot=NA, 
                                   nMaj1_A=ynMaj_1, nMin1_A=0, frac1_A=yfrac_1, nMaj2_A=ynMaj_2, nMin2_A=NA, frac2_A=yfrac_2, SDfrac_A=NA, SDfrac_A_BS=NA, frac1_A_0.025=NA, frac1_A_0.975=NA, 
                                   nMaj1_B=NA, nMin1_B=NA, frac1_B=NA, nMaj2_B=NA, nMin2_B=NA, frac2_B=NA, SDfrac_B=NA, SDfrac_B_BS=NA, frac1_B_0.025=NA, frac1_B_0.975=NA, 
                                   nMaj1_C=NA, nMin1_C=NA, frac1_C=NA, nMaj2_C=NA, nMin2_C=NA, frac2_C=NA, SDfrac_C=NA, SDfrac_C_BS=NA, frac1_C_0.025=NA, frac1_C_0.975=NA, 
                                   nMaj1_D=NA, nMin1_D=NA, frac1_D=NA, nMaj2_D=NA, nMin2_D=NA, frac2_D=NA, SDfrac_D=NA, SDfrac_D_BS=NA, frac1_D_0.025=NA, frac1_D_0.975=NA, 
                                   nMaj1_E=NA, nMin1_E=NA, frac1_E=NA, nMaj2_E=NA, nMin2_E=NA, frac2_E=NA, SDfrac_E=NA, SDfrac_E_BS=NA, frac1_E_0.025=NA, frac1_E_0.975=NA, 
                                   nMaj1_F=NA, nMin1_F=NA, frac1_F=NA, nMaj2_F=NA, nMin2_F=NA, frac2_F=NA, SDfrac_F=NA, SDfrac_F_BS=NA, frac1_F_0.025=NA, frac1_F_0.975=NA)
  )
  return(subclone.data)
}


##############################################
# Copy number preparations
##############################################
#' This function takes a Battenberg styled data.frame, pulls out the A copy number profile,
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
#' @param samplename A String containing the samplename
#' @param cndata A data.frame in Battenberg subclones format
#' @param gender Specify either the string male or female
#' @param outfile A String with the full path to where the output should be written
#' @param ploidy An optional parameter that is used to call gains/losses against. If higher than the ploidy a segment is considered a gain, lower a loss. When left as NA the sample ploidy is used (Default NA).
#' @param exclude_sex_chroms Optional parameter whether to exclude sex chromosomes (Default FALSE)
#' @author tjm
#' @export
collate_bb_subclones = function(samplename, cndata, gender, outfile, ploidy=NA, exclude_sex_chroms=F) {
  
  if(gender == 'male' | gender == 'Male') {
    is_male = T
  } else if(gender == 'female' | gender == 'Female') {
    is_male = F
  } else {
    stop("Unknown gender supplied, exit.")
  }
  
  allsegs = NULL
  cndata = cbind(samplename, cndata)
  # Remove the sex chromosomes if required
  if (exclude_sex_chroms) {
    cndata = cndata[!(cndata$chr=="X" | cndata$chr=="Y"),]
  }
  
  # Remove all segments without a call
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
  
  # Ploidy weighted by segment size - if not provided as input
  if (is.na(ploidy)) {
    ploidy = round(sum(cndataweighted)/2)*2
  }
  
  # If there is no column with a tumour name, add it in temporarily
  if (! "Tumour_Name" %in% colnames(cndata)) {
    cndata = data.frame(Tumour_Name=samplename, cndata)
  }
  
  allsegs = data.frame(cndata[, c("Tumour_Name", "chr", "startpos", 
                                  "endpos", "nMaj1_A", "nMin1_A", 
                                  "frac1_A", "nMaj2_A", "nMin2_A", 
                                  "frac2_A", "SDfrac_A")], tumour_ploidy=ploidy)
  
  ######################################################
  # Now classify all segments into a category
  ######################################################
  
  # Columns to be selected
  select_state_1 = c(1:7,11:12)
  select_state_2 = c(1:4,8:12)
  
  tot = dim(allsegs)[1]
  # if you have clonal LOH, subclonal LOH is not counted!!!  Losses not counted
  allsegsa <- NULL
  CNA <- NULL
  for (i in 1:dim(allsegs)[1]) {
    
    is_hd = allsegs$nMin1_A[i] == 0 & allsegs$nMaj1_A[i] == 0
    is_loh_1 = xor(allsegs$nMin1_A[i] == 0, allsegs$nMaj1_A[i] == 0)
    if (!is.na(allsegs$nMaj2_A[i])) {
      is_loh_2 = xor(allsegs$nMin2_A[i] == 0, allsegs$nMaj2_A[i] == 0)
    }
    

    ######################################################
    # Determine the logic category parameters
    ######################################################
    
    # X/Y are a special case in males
    is_sex_chrom = (allsegs$chr[i]=="X" | allsegs$chr[i]=="Y")
    if (is_male & is_sex_chrom) {
      sex_chromosome_maj_exp = ploidy/2
      sex_chromosome_maj_amp = ploidy*2
      
      # check that min==0, if not, then switch
      if (allsegs$nMaj1_A[i]==0 & allsegs$nMin1_A[i] > 0) {
        allsegs$nMaj1_A[i] = allsegs$nMin1_A[i]
        allsegs$nMin1_A[i] = 0
        
        if (!is.na(allsegs$nMaj2_A[i])) {
          allsegs$nMaj2_A[i] = allsegs$nMin2_A[i]
          allsegs$nMin2_A[i] = 0
        }
      }
      
      is_normal = allsegs$nMaj1_A[i] == sex_chromosome_maj_exp
      is_gained_1 = allsegs$nMaj1_A[i] > sex_chromosome_maj_exp
      is_amplified_1 = allsegs$nMaj1_A[i] > sex_chromosome_maj_amp
      is_lost_1 = allsegs$nMaj1_A[i] < sex_chromosome_maj_exp
      
      # Only required when subclonal copy number
      if (!is.na(allsegs$nMaj2_A[i])) {
        is_gained_2 = allsegs$nMaj2_A[i] > sex_chromosome_maj_exp
        is_amplified_2 = allsegs$nMaj2_A[i] > sex_chromosome_maj_amp
        is_lost_2 = allsegs$nMaj2_A[i] < sex_chromosome_maj_exp
        
        is_gained_min = F #(allsegs$nMin2_A[i] > ploidy/2 & allsegs$nMin1_A[i] != allsegs$nMin2_A[i])
        is_gained_maj = (allsegs$nMaj2_A[i] > sex_chromosome_maj_exp & allsegs$nMaj1_A[i] != allsegs$nMaj2_A[i])
        
        is_amplified_min = F #(allsegs$nMin2_A[i] > ploidy*2 & allsegs$nMin1_A[i] != allsegs$nMin2_A[i])
        is_amplified_maj = (allsegs$nMaj2_A[i] > sex_chromosome_maj_amp & allsegs$nMaj1_A[i] != allsegs$nMaj2_A[i])

        is_lost_min = F #(allsegs$nMin1_A[i] < ploidy/2 & allsegs$nMin1_A[i] != allsegs$nMin2_A[i])        
        is_lost_maj = (allsegs$nMaj1_A[i] < sex_chromosome_maj_exp & allsegs$nMaj1_A[i] != allsegs$nMaj2_A[i])
      }
      
    } else {
     
      is_normal = (allsegs$nMaj1_A[i] == ploidy/2) & (allsegs$nMin1_A[i] == ploidy/2)
      is_gained_1 = (allsegs$nMaj1_A[i] > ploidy/2 | allsegs$nMin1_A[i] > ploidy/2)
      is_amplified_1 = (allsegs$nMaj1_A[i] > ploidy*2) | (allsegs$nMin1_A[i] > ploidy*2)
      is_lost_1 = (allsegs$nMaj1_A[i] < ploidy/2) | (allsegs$nMin1_A[i] < ploidy/2)
      
      # Only required when subclonal copy number
      if (!is.na(allsegs$nMaj2_A[i])) {
        is_gained_2 = (allsegs$nMaj2_A[i] > ploidy/2 | allsegs$nMin2_A[i] > ploidy/2)
        is_amplified_2 = (allsegs$nMaj2_A[i] > ploidy*2 | allsegs$nMin2_A[i] > ploidy*2)
        is_lost_2 = (allsegs$nMaj2_A[i] < ploidy/2 | allsegs$nMin2_A[i] < ploidy/2)
        
        is_gained_min = (allsegs$nMin2_A[i] > ploidy/2 & allsegs$nMin1_A[i] != allsegs$nMin2_A[i])
        is_gained_maj = (allsegs$nMaj2_A[i] > ploidy/2 & allsegs$nMaj1_A[i] != allsegs$nMaj2_A[i])
        
        is_amplified_min = (allsegs$nMin2_A[i] > ploidy*2 & allsegs$nMin1_A[i] != allsegs$nMin2_A[i])
        is_amplified_maj = (allsegs$nMaj2_A[i] > ploidy*2 & allsegs$nMaj1_A[i] != allsegs$nMaj2_A[i])
        
        is_lost_maj = (allsegs$nMaj1_A[i] < ploidy/2 & allsegs$nMaj1_A[i] != allsegs$nMaj2_A[i])
        is_lost_min = (allsegs$nMin1_A[i] < ploidy/2 & allsegs$nMin1_A[i] != allsegs$nMin2_A[i])
      }
    }

    ######################################################
    # Actual logic to perform the classification
    ######################################################
    # Clonal copy number
    if (is.na(allsegs$nMaj2_A[i])) {
      if (is_hd & !is_sex_chrom) {
        allsegsa <- rbind(allsegsa, allsegs[i,select_state_1])
        CNA <- c(CNA, "cHD")
      }
      if (is_loh_1 & !is_sex_chrom) {
        allsegsa <- rbind(allsegsa, allsegs[i,select_state_1])
        CNA <- c(CNA, "cLOH")
      }
      if (is_gained_1) {
        if (is_amplified_1) {
          allsegsa <- rbind(allsegsa, allsegs[i,select_state_1])
          CNA <- c(CNA, "cAmp")
        }
        else {
          allsegsa <- rbind(allsegsa, allsegs[i,select_state_1])
          CNA <- c(CNA, "cGain")
        }
      }
      if (is_normal) {
        allsegsa <- rbind(allsegsa, allsegs[i,select_state_1])
        CNA <- c(CNA, "NoCNV")
      }
      if (is_lost_1) {
        allsegsa <- rbind(allsegsa, allsegs[i,select_state_1])
        CNA <- c(CNA, "cLoss")
      }
    }
    
    # Subclonal copy number
    if (!is.na(allsegs$nMaj2_A[i])) {
      if (is_hd & !is_sex_chrom) {
        CNA <- c(CNA, "sHD")
        allsegsa <- rbind(allsegsa, allsegs[i,select_state_1])
      }
      if (is_loh_1 & is_loh_2 & !is_sex_chrom) {
        CNA <- c(CNA, "cLOH")
        tmp <- allsegs[i,select_state_2]
        names(tmp) <- names(allsegs[i,select_state_1])
        allsegsa <- rbind(allsegsa, tmp[,names(allsegs[i,select_state_1])])
      }
      else if (is_loh_1 & !is_sex_chrom) {
        CNA <- c(CNA, "sLOH")
        allsegsa <- rbind(allsegsa, allsegs[i,select_state_1])
      }
      if (is_gained_1 & is_gained_2) {
        if (is_amplified_1 & is_amplified_2) {
          CNA <- c(CNA, "cAmp")
          allsegsa <- rbind(allsegsa, allsegs[i,select_state_1])
        }
        else {
          CNA <- c(CNA, "cGain")
          allsegsa <- rbind(allsegsa, allsegs[i,select_state_1])
        }
      }
      if (is_gained_min | is_gained_maj) {
        if (is_amplified_min | is_amplified_maj) {
          CNA <- c(CNA, "sAmp")
          tmp <- allsegs[i,select_state_2]
          names(tmp) <- names(allsegs[i,select_state_1])
          allsegsa <- rbind(allsegsa, tmp[,names(allsegs[i,select_state_1])])
        }
        else {
          CNA <- c(CNA, "sGain")
          tmp <- allsegs[i,select_state_2]
          names(tmp) <- names(allsegs[i,select_state_1])
          allsegsa <- rbind(allsegsa, tmp[,names(allsegs[i,select_state_1])])
        }
      }
      if (is_lost_1 & is_lost_2) {
        CNA <- c(CNA, "cLoss")
        tmp <- allsegs[i,select_state_2]
        names(tmp) <- names(allsegs[i,select_state_1])
        allsegsa <- rbind(allsegsa, tmp[,names(allsegs[i,select_state_1])])
      }
      if ((is_lost_maj | is_lost_min)) {
        allsegsa <- rbind(allsegsa, allsegs[i,select_state_1])
        CNA <- c(CNA, "sLoss")      
      }
    }

    if(i %% 100 ==0){
      print(paste(i,"/",tot))
    }
  }  
  
  allsegsa <- cbind(allsegsa, CNA)
  allsegsa$frac1_A[allsegsa$CNA == "cGain" | allsegsa$CNA == "cAmp" | allsegsa$CNA == "cLOH" | allsegsa$CNA == "cHD" | allsegsa$CNA == "cLoss"] = 1
  
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
