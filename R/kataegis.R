#' Identify Kataegis events in SNV mutation data
#'
#' @param samplename Samplename used when writing output files
#' @param dpInfile A DP input file
#' @param outdir Directory where output will be written
#' @param gamma_param Gamma parameter to be used for segmentation
#' @param kmin Kmin parameter to be used for segmentation
#' @param kataegis.threshold Intermutation distance, if NA will be set to: (a) 900000 > SNVs: 100,
#' (b) 500000 > SNVs: 250, (c) 100000 > SNVs, (d) otherwise 1000 (Default NA)
#' @param minMuts Minimum number of mutations within kataegis.threshold to call Kataegis
#' @param logScale Transform intermutation distance to logscale
#' @param makePlots Make a figure for each Kataegis event
#' @param removeFromDPin Remove SNVs identified as part of a Kataegis event from the DP input file
#' @author dw9
#' @export
identifyKataegis<-function(samplename, dpInfile, outdir=".", gamma_param=25, kmin=2, kataegis.threshold=NA, minMuts=6, logScale=F, makePlots=F, removeFromDPin=F){
  data = read.table(dpInfile, header=T, stringsAsFactors=F)
  # It is possible to add move columns here: The first 2 columns in chrpos must contain chr and pos. If there are additional columns, these will be included in the kataegis_mutations output
  chrpos = data[,c(1,3)]

  # No threshold supplied, therefore use these heuristics
  if (is.na(kataegis.threshold)) {
    no.muts = nrow(data)
    print(paste("#muts=",no.muts,sep=""))
    if(no.muts>900000){
      kataegis.threshold = 100
    }else if(no.muts>500000){
      kataegis.threshold = 250
    }else if(no.muts>100000){
      kataegis.threshold = 500
    }else{
      kataegis.threshold = 1000
    }
  }

  print(samplename)
  segmentedIntermutationDistances = array(NA,c(0,3))
  kat.regions.all = array(NA,c(0,5))
  for(chr in unique(chrpos[,1])) { #c(1:22,"X","Y")){
    print(paste("chr ",chr,sep=""))
    positions = chrpos[chrpos[,1]==chr,2]
    if(length(positions)>1){
      positions = sort(positions)
      if(logScale){
        interDist = log(positions[-1]-positions[-(length(positions))])
      }else{
        interDist = positions[-1]-positions[-(length(positions))]
      }
      saved.positions = positions
      #use midpoints for plotting
      positions = (positions[-1]+positions[-(length(positions))])/2
      sdev <- getMad(interDist,k=25)
      #res= selectFastPcf(interDist,kmin,gamma_param*sdev,T)
      res= exactPcf(interDist,kmin,gamma_param*sdev,T)
      segInterDist = res$yhat
      if(makePlots){
        makeKataegisPlot(outdir, samplename, chr, interDist, segInterDist, gamma_param, kmin)
      }
      segmentedIntermutationDistances = rbind(segmentedIntermutationDistances,cbind(chr,positions,segInterDist))

      #use positions, rather than midpoints, to define extent of kataegis
      positions = saved.positions
      if(logScale){
        katLoci = segInterDist<=log(kataegis.threshold)
      }else{
        katLoci = segInterDist<=kataegis.threshold
      }
      if(sum(katLoci)>0){
        start.regions = which(c(katLoci, F) & !c(F, katLoci))
        end.regions = which(!c(katLoci, F) & c(F, katLoci))
        if(length(start.regions)>0){
          for(r in 1:length(start.regions)){
            kat.regions.all = rbind(kat.regions.all,
              c(chr,positions[start.regions[r]],positions[end.regions[r]],positions[end.regions[r]]-positions[start.regions[r]],end.regions[r]-start.regions[r]+1))
          }
        }
      }
    }
  }

  if(logScale){
    outfile_suffix = "_logScale.txt"
  }else{
    outfile_suffix = ".txt"
  }

  colnames(kat.regions.all) = c("chr", "start", "end", "size", "no_of_variants")
  kat.regions.all = kat.regions.all[as.numeric(kat.regions.all[,5])>=minMuts,]
  # single row causes problems
  if(length(kat.regions.all)==5){
    temp = array(kat.regions.all,c(1,5))
    colnames(temp) = names(kat.regions.all)
    kat.regions.all = temp
  }
  write.table(kat.regions.all,paste(outdir,"/",samplename,"_kataegis_regions_fromSegmentation_gamma",gamma_param,"_threshold",kataegis.threshold,"_kmin",kmin,"_minMuts",minMuts,outfile_suffix,sep=""),sep="\t",quote=F,row.names=F,col.names=T)

  colnames(segmentedIntermutationDistances)=c("chr","pos","intermutation.distances")
  write.table(segmentedIntermutationDistances,paste(outdir,"/",samplename,"_kataegis_segmentedIntermutationDistances_gamma",gamma_param,"_threshold",kataegis.threshold,"_kmin",kmin,"_minMuts",minMuts,outfile_suffix,sep=""),sep="\t",quote=F,row.names=F)

  kat.muts = chrpos[0,]
  if(nrow(kat.regions.all)>0){
    for(k in 1:nrow(kat.regions.all)){
      #print(paste(kat.regions.all[k,5],sum(chrpos[,1]==kat.regions.all[k,1] & chrpos[,2]>=as.numeric(kat.regions.all[k,2]) & chrpos[,2]<=as.numeric(kat.regions.all[k,3])),sep=","))
      kat.muts = rbind(kat.muts,chrpos[chrpos[,1]==kat.regions.all[k,1] & chrpos[,2]>=as.numeric(kat.regions.all[k,2]) & chrpos[,2]<=as.numeric(kat.regions.all[k,3]),])
    }
  }
  write.table(kat.muts,paste(outdir,"/",samplename,"_kataegis_mutationsWithin_Kataegis_gamma",gamma_param,"_threshold",kataegis.threshold,"_kmin",kmin,"_minMuts",minMuts,".txt",sep=""),sep="\t",quote=F,row.names=F)
  if(nrow(kat.muts)!=sum(as.numeric(kat.regions.all[,5]))){
    stop("number of kataegis variants is inconsistent")
  }

  #' Remove the Kataegis SNVs from the DP input file
  if (removeFromDPin) {
    # Filter the DP input file
    chrpos_kat = paste(kat.muts[,1], kat.muts[,2], sep="_")
    chrpos_data = paste(data[,1], data[,3], sep="_")
    select = !(chrpos_data %in% chrpos_kat)
    data = data[select,]
    write.table(data, file=dpInfile, quote=F, sep="\t", row.names=F)
  }
}

#' Make Kataegis plot
#' @noRd
makeKataegisPlot = function(outdir, samplename, chr, interDist, segInterDist, gamma_param, kmin) {
  if(logScale){
    png(filename = paste(outdir,"/",samplename,"_chr",chr,"_kataegis_segmentedIntermutationDistance_logScale_gamma",gamma_param,"_kmin",kmin,".png",sep=""), width = 2000, height = 1000, res = 200)
    par(mar = c(5,5,5,0.5), cex = 0.6, cex.main=3, cex.axis = 2, cex.lab = 2)
    plot(c(min(positions)/1000000,max(positions)/1000000),c(0,ceiling(max(interDist))),pch=".",type = "n",
         main = paste("Chromosome ", chr, sep=""), xlab = "Position (Mb)", ylab = "log(intermutation distance) (segmented)")
    points(positions/1000000,interDist,pch=20,col="red",cex=0.5)
    points(positions/1000000,segInterDist,pch=20,cex=2,col="green")
  }else{
    png(filename = paste(outdir,"/",samplename,"_chr",chr,"_kataegis_segmentedIntermutationDistance_gamma",gamma_param,"_kmin",kmin,".png",sep=""), width = 2000, height = 1000, res = 200)
    par(mar = c(5,5,5,0.5), cex = 0.6, cex.main=3, cex.axis = 2, cex.lab = 2)
    plot(c(min(positions)/1000000,max(positions)/1000000),c(0,ceiling(log(max(interDist)))),pch=".",type = "n",
         main = paste("Chromosome ", chr, sep=""), xlab = "Position (Mb)", ylab = "log(intermutation distance) (segmented)")
    points(positions/1000000,log(interDist),pch=20,col="red",cex=0.5)
    points(positions/1000000,log(segInterDist),pch=20,cex=2,col="green")
  }
}
