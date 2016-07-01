###
### CGH pipeline source file
###

### Functions from CGH_useful_function.R

##Useful CGH manipulation functions

plot.copynumber.calls <- function(case_LogR_ratioCalled) {
  
png(paste0(sampleNames(case_LogR_ratioCalled)[1],".calls.png"), height=720, width=1784)
par(mar=c(5, 6, 4, 2) + 0.1)
chr.total <- table(chromosomes(case_LogR_ratioCalled))
chrs      <- unique(chromosomes(case_LogR_ratioCalled))

call.colours <- cbind(calls(case_LogR_ratioCalled), NA)

for(i in 1:nrow(call.colours)) {
  if(call.colours[i,1]=="0") {
    call.colours[i,2] <- "darkgrey"
  }
  if(call.colours[i,1]=="1") {
    call.colours[i,2] <- "darkgreen"
  }
  if(call.colours[i,1]=="2") {
    call.colours[i,2] <- "green"
  }
  if(call.colours[i,1]=="-1") {
    call.colours[i,2] <- "darkred"
  }
  if(call.colours[i,1]=="-2") {
    call.colours[i,2] <- "red"
  }
}

plot(copynumber(case_LogR_ratioCalled),
     pch=20,
     cex=0.02,
     ylim=c(-1,1),
     col=call.colours[,2],
     xaxt="n",
     xlab="Chromosomes",
     ylab="LogR.ratio",
     cex.axis=2,
     cex.lab=2,
     cex.main=2,
     main=sampleNames(case_LogR_ratioCalled)[1])

prev.chr <- NULL
abline(v=0, lty="dotted")

for (chr in min(chrs):max(chrs)) {
  
  current.chr <- sum(prev.chr, chr.total[chr])
  text.offset <- chr.total[chr]/2
  abline(v=current.chr, lty="dotted")
  text(current.chr-text.offset, -1, paste0(chr), cex=1)
  prev.chr <- current.chr
  
}

abline(h=0, lty="dotted")

dev.off()

}

regions.as.table <- function(region, as.df=TRUE) {
  ##Takes a CGH compressed object and converts it to a dataframe
  ##If as.df is set to FALSE it returns a matrix
  ##Gives a weird warning message which doesn't impact the results
  output <- cbind(chromosomes(region), bpstart(region), bpend(region), regions(region))
  colnames(output)[1:3] <- c("Chr", "Start", "End")
  if(as.df) {
    output <- as.data.frame(output)
  }
  return(output)
}

subset.chromosome <- function(CGHobject, chr=NULL) {
  ##Takes a CGHobject and subsets a chromosome from it
  ##chr must be set, NULL returns error
  output <- CGHobject[chromosomes(CGHobject)==chr]
  return(output)
}

std_CGHplot <- function(x, y, y.lim=c(-1,1), char=".", colour="black") {
  ##plots x,y using normal parameters for looking at LogR data and calls
  plot(x, y, ylim=y.lim, pch=char, col=colour)
}

###Quality control functions

##mBAF validation functions for G,A & L

mBAF_based_gain_QC <- function(regions.table, mBAF.table, min.baf=0.58) {
  ##takes the calls in the form of a dataframe (col.4 as the start column)
  ##queries this against a dataframe in the same format containing average bafs
  ##changes all gains which have a baf less than min.baf back to zero (no call)
  output <- list()
  length(output) <- 2
  
  for (i in 4:ncol(regions.table)) {
  index <- which(regions.table[,i]==1)
  bad.baf <- which(mBAF.table[index,i] < min.baf)
  
  if(length(bad.baf)>0) {
  regions.table[index[bad.baf],i] <- 0
  }
  
  }
  
  output[[1]] <- regions.table
  output[[2]] <- mBAF.table
  names(output) <- c("regions.filtered", "mBaf.filtered")
  return(output)
  
}

mBAF_based_loss_QC <- function(regions.table, mBAF.table, min.baf=0.64) {
  ##takes the calls in the form of a dataframe (col.4 as the start column)
  ##queries this against a dataframe in the same format containing average bafs
  ##changes all losses which have a baf less than min.baf back to zero (no call)
  output <- list()
  length(output) <- 2
  
  for (i in 4:ncol(regions.table)) {
    index <- which(regions.table[,i]==-1)
    bad.baf <- which(mBAF.table[index,i] < min.baf)
    
    if(length(bad.baf)>0) {
      regions.table[index[bad.baf],i] <- 0
    }
    
  }
  
  output[[1]] <- regions.table
  output[[2]] <- mBAF.table
  names(output) <- c("regions.filtered", "mBaf.filtered")
  return(output)
  
}

mBAF_based_amp_QC <- function(regions.table, mBAF.table, min.baf=0.595) {
  ##takes the calls in the form of a dataframe (col.4 as the start column)
  ##queries this against a dataframe in the same format containing average bafs
  ##changes all amps which have a baf less than min.baf back to zero (no call)
  output <- list()
  length(output) <- 2
  
  for (i in 4:ncol(regions.table)) {
    index <- which(regions.table[,i]==2)
    bad.baf <- which(mBAF.table[index,i] < min.baf)
    
    if(length(bad.baf)>0) {
      regions.table[index[bad.baf],i] <- 0
    }
    
  }
  
  output[[1]] <- regions.table
  output[[2]] <- mBAF.table
  names(output) <- c("regions.filtered", "mBaf.filtered")
  return(output)
  
}

### CGH pipeline main functions

##Taken from CGH_pipeline_20141006.R

movav_normalise <- function(sample, control, weight.lim=3) {
  ##Normalises samples LogR values based on the moving average of all normals
  ##Both objects must be same length
  popvar <- function(x, na.rm = FALSE) {
    # Calculates population variance instead of sample variance (which is the
    # default of the var() function in R).
    #
    # Args:
    # x: a vector of the population data.
    # na.rm: a logical value indicating whether NA values should be stripped
    # before the computation proceeds.
    #
    # Returns:
    # The population variance.
    if(na.rm) {
      x <- x[!is.na(x)]
    } else if(any(is.na(x))) {
      return(NA)
    }
    mean((x-mean(x))^2)
  }
  
  r <- NULL
  weights <- seq(-weight.lim,weight.lim, by=0.1)
  n <- 0
  
  for (i in weights) {
    n          <- n+1
    norm       <- sample - (control*i)
    r[n] <- popvar(norm, na.rm=TRUE)
  }
  
  results <- data.frame(weights, r)
  
  final_weight <- results[which(results$r==min(results$r)),1]
  
  output <- sample - (control*final_weight)
  
  cat("final weight is:",final_weight,"\n")
  
  return(output)
}

CGHregions_GC_edit <- function(input, critfound=2) {
  #The main CGHregion function extracted to be used without the av.error calculation and for c to be determined
  CGHdata <- input
  kolnam  <- sampleNames(CGHdata)
  ncolm   <- ncol(CGHdata)
  chromo  <- chromosomes(CGHdata)
  bppos   <- bpstart(CGHdata)
  numbr   <- nrow(CGHdata)
  CGHdata <- data.frame(chromosomes(CGHdata), bpstart(CGHdata), calls(CGHdata))
  
  normstate   <- 0
  levels      <- c(-2, -1, 0, 1, 2)
  
  cspr        <- CGHregions:::.deterreg(CGHdata=CGHdata, critfound, ncolm, normstate, levels)
  
  regionsfound <- cspr[[2]]
  countnomono  <- cspr[[3]]
  
  print(paste("c = ",critfound,", nr of regions: ", nrow(regionsfound), sep=""))
  
  res     <- apply(regionsfound, 1, CGHregions:::.whichsign2, ctdat=countnomono, levels=levels)
  prof    <- t(sapply(res, function(x) {as.vector(x[[1]], mode="numeric")}))
  nclone  <- apply(regionsfound,1,CGHregions:::.ntd,ct=countnomono)
  aved    <- signif(as.vector(lapply(res, function(x) {x[[2]]}), mode="numeric")/nclone, digits=3)
  bp      <- t(apply(regionsfound, 1, CGHregions:::.findbp, bppos = bppos))
  chrreg  <- apply(regionsfound, 1, CGHregions:::.findchr, chr = chromo)
  towrite <- cbind(regionsfound, bp, chrreg, nclone, aved, prof)
  rownames(towrite) <- c()
  od      <- order(towrite[,1])
  towrite <- towrite[od, -(1:2)]
  kolnamnew <- c("bp start", "bp end", "chromosome", "nclone", "Ave Dist", kolnam)
  colnames(towrite) <- kolnamnew
  
  annotation  <- data.frame(Chromosome=towrite[,3], Start=towrite[,1], End=towrite[,2], Nclone=towrite[,4], AveDist=towrite[,5])
  metadata    <- data.frame(  labelDescription=c("Chromosomal position",
                                                 "Basepair position start",
                                                 "Basepair position end",
                                                 "Number of clones in region",
                                                 "Average distance"),
                              row.names=c("Chromosome",
                                          "Start",
                                          "End",
                                          "Nclone",
                                          "AveDist")
  )
  
  dimLabels   <- c("featureNames", "featureColumns")
  annotation  <- new("AnnotatedDataFrame", data=annotation, dimLabels=dimLabels, varMetadata=metadata)
  result      <- new("cghRegions", regions=as.matrix(towrite[,6:ncol(towrite)]), featureData=annotation)
  
  return(result)
}

plot_post_CGHregions_filtering <- function(test_data, col.start=4, save=FALSE, title="post_CGHregions_filtering") {
  ##plots a dataframe of calls similarly to the plot function in CGHregions
  ##colours represent average calls in these regions - 0 is not necessarily no changes, they could average out
  ##col.start is important for telling the function where the first column of the first sample is
  ##it then assumes every column to the end is another array
  weight <- rep(c(0,0.25), times=ceiling(nrow(test_data)/2))[1:nrow(test_data)]
  
  test_data <- data.frame(test_data, test_data[,1]+weight)
  
  test_data <- data.frame(test_data, rowMeans(test_data[,col.start:(ncol(test_data)-1)]))
  
  colours <- NULL
  
  reds <- c("firebrick4", "firebrick3", "firebrick2", "firebrick1")
  reds <- rev(reds)
  greens <- c("chartreuse4", "chartreuse3", "chartreuse2", "chartreuse1")
  
  for(i in 1:nrow(test_data)) {
    
    if(test_data[i,ncol(test_data)]< -0.75) {colours[i] <- reds[1]
    } else if (test_data[i,ncol(test_data)]>= -0.75 & test_data[i,ncol(test_data)]< -0.5) {colours[i] <- reds[2]
    } else if (test_data[i,ncol(test_data)]>= -0.5 & test_data[i,ncol(test_data)]< -0.25) {colours[i] <- reds[3]
    } else if (test_data[i,ncol(test_data)]>= -0.25 & test_data[i,ncol(test_data)]< 0) {colours[i] <- reds[4]
    } else if (test_data[i,ncol(test_data)]==0) {colours[i] <- "black"
    } else if (test_data[i,ncol(test_data)]<= 0.25 & test_data[i,ncol(test_data)]> 0) {colours[i] <- greens[1]
    } else if (test_data[i,ncol(test_data)]<= 0.5 & test_data[i,ncol(test_data)]> 0.25) {colours[i] <- greens[2]
    } else if (test_data[i,ncol(test_data)]<= 0.75 & test_data[i,ncol(test_data)]> 0.5) {colours[i] <- greens[3]
    } else if (test_data[i,ncol(test_data)]>= 0.75) {colours[i] <- greens[4]}}
  
  if(save){
    pdf(paste0(title,".pdf"))
  }
  
  for (i in 1:nrow(test_data)) {
    
    x <- seq(test_data[i,2], test_data[i,3], by=100000)
    plot(x, rep(test_data[i,(ncol(test_data)-1)], times=length(x)), pch=20, ylim=c(0,24), xlim=c(0,250000000), xlab="", ylab="", xaxt="n", yaxt="n", col=colours[i])
    par(new="TRUE")
  }
  
  axis(2, at=seq(0,24, by=1))
  axis(1, at=seq(0,250000000, by=10000000))
  title(main=paste0(title))
  if(save){
    dev.off()
  }
}

mean_mBAF_per_region <- function(regions_df, mBAF_df, chr_num=23, info_loss_threshold=0.9) {
  ##calculates, for a patient, the average mbaf in this region in every array
  ##chr_num must be set to the number of chromosomes you are actually analysing  
  mean_mBAF_per_region <- NULL
  
  for(chr in 1:chr_num) {
    mBAF_df_chr_specific    <- mBAF_df[mBAF_df$Chr==chr,]
    regions_df_chr_specific <- regions_df[regions_df[,"Chr"]==chr,]
    for(i in 1:nrow(regions_df_chr_specific)){
      region_mBAFs <- mBAF_df_chr_specific[which(mBAF_df_chr_specific$Position>=regions_df_chr_specific[,"Start"][i] & mBAF_df_chr_specific$Position<=regions_df_chr_specific[,"End"][i]),]
      
      mean_mBAF_indiv_region <- NULL
      
      for(k in 4:ncol(region_mBAFs)){
        
        mean_mBAF_indiv_region[k] <- mean(na.omit(region_mBAFs[,k]))
        
        if(mean_mBAF_indiv_region[k]=="NaN"){mean_mBAF_indiv_region[k]<-1}
        if(sum(is.na(region_mBAFs[,k])) / sum(table(is.na(region_mBAFs[,k]))) >= info_loss_threshold){
          mean_mBAF_indiv_region[k]<-1
          print(paste0("Region with ",(1-info_loss_threshold)*100,"% or less SNP information is reverted to 1"))
        }
        
      }
      
      mean_mBAF_indiv_region[1] <- chr
      mean_mBAF_indiv_region[2] <- regions_df_chr_specific[,"Start"][i]
      mean_mBAF_indiv_region[3] <- regions_df_chr_specific[,"End"][i]
      
      mean_mBAF_per_region <- rbind(mean_mBAF_per_region, mean_mBAF_indiv_region)
      
    }
  }
  
  rownames(mean_mBAF_per_region) <- NULL
  return(mean_mBAF_per_region)
}