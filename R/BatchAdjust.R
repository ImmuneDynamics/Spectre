#' BatchAdjust
#' 
#' @export

####         
#
# BatchAdjust.R   ## 
#
# main function : BatchAdjust()
#


# Install flowCore:
#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("flowCore")

library("flowCore") # for read.FCS

# Define asinh factor:
#g_asinh_b <- 1/5; # global asinh factor; transform is asinh(x*b); corresponds to a_b in cytofkit cytof_exprsExtract()
g_asinh_b <- 1; # global asinh factor; transform is asinh(x*b); corresponds to a_b in cytofkit cytof_exprsExtract()

# logToFile
# Maintain a log.
logToFile <- function(logfilename, logmessage, timestamp=FALSE, echo=TRUE, overwrite=FALSE){
  # Should we add a timestamp?
  if(timestamp){
    logmessage <- sprintf("%s  %s", Sys.time(), logmessage);
  }
  if(overwrite==TRUE){
    # Over write existing.
    logfile <- file(logfilename, open="wt");
  }else{
    # Open the file for appending, in text mode.
    logfile <- file(logfilename, open="at");
  }
  writeLines(logmessage, con=logfile);
  close(logfile);
  if(echo){
    print(logmessage, quote=FALSE);
  }
}

# get80thPercentile
# Return the requested (80th) percentile of the vector vec.
#  perc <- .8; # percentile. median would be .5. Highest % zeros in any sample is 76.35%
get80thPercentile <- function(vec, perc=.8){
  npoints <- length(vec);
  vec_sorted <- sort(vec);
  perci <- ceiling(perc * npoints);
  perc_value <- vec_sorted[perci];
  return(perc_value);
}


# listBatchesPresent
# Find all anchor files in basedir,
# parse file names to return a vector of batch numbers.
#    xxx[batchKeyword][##][anchorKeyword]xxx.fcs
# example: 011118_Barcode_7_anchor stim.fcs:  anchorKeyword="anchor stim", batchKeyword="Barcode_"
# example: Set10_CTT0.fcs:  anchorKeyword="CTT0"; batchKeyword="Set";
listBatchesPresent <- function(basedir, batchKeyword="Barcode_", anchorKeyword="anchor stim"){
  batches_present <- c();
  basedir_escapeSpace <- gsub(pattern=" ", replacement="\\ ", x=basedir, fixed=TRUE);
  grepForAnchor_escapeSpace <- gsub(pattern=" ", replacement="\\ ", x=anchorKeyword, fixed=TRUE);
  
  ls_cmd <- sprintf("ls -1 %s/*%s*.fcs", basedir_escapeSpace, grepForAnchor_escapeSpace);
  anchors_list <- sort(system(ls_cmd, intern=TRUE, ignore.stdout = FALSE, ignore.stderr = TRUE, wait = TRUE));
  
  underscore_anchorKeyword <- sprintf("_%s", anchorKeyword);
  for(ananchor in anchors_list){
    first_part <- unlist(strsplit(ananchor, underscore_anchorKeyword, fixed=TRUE));
    next_part <- unlist(strsplit(first_part[1], batchKeyword, fixed=TRUE));
    this_batch_num <- as.numeric(next_part[2]); # numeric
    batches_present <- c(batches_present, this_batch_num);
  }
  return(sort(batches_present));
}


# getMinEventCount
# Get the minimum number of events across anchor files.
getMinEventCount <- function(anchorKeyword="anchor stim", basedir="/Users/ron-home/projects/DATA/SLE_Malaria/Bead_Normalized_Debarcoded_Singlet_JG_unzipped"){
  
  #T0 <- Sys.time();
  whichlines <- NULL;
  basedir_escapeSpace <- gsub(pattern=" ", replacement="\\ ", x=basedir, fixed=TRUE);
  grepForAnchor_escapeSpace <- gsub(pattern=" ", replacement="\\ ", x=anchorKeyword, fixed=TRUE);
  ls_cmd <- sprintf("ls -1 %s/*%s*.fcs", basedir_escapeSpace, grepForAnchor_escapeSpace);
  anchors_list <- sort(system(ls_cmd, intern=TRUE, ignore.stdout = FALSE, ignore.stderr = TRUE, wait = TRUE));
  anchor_counter <- 0;
  minCount <- Inf;
  for(ananchor in anchors_list){
    anchor_counter <- anchor_counter + 1;
    thisFCSobject <- read.FCS(ananchor, transformation=NULL, which.lines=whichlines);
    #thisFCSobject <- read.FCS(ananchor, transformation=NULL, which.lines=whichlines, column.pattern="Time");
    #this_data <- exprs(thisFCSobject);
    Nrows <- nrow(thisFCSobject);
    if(Nrows < minCount){
      minCount <- Nrows;
    }
    #print(sprintf("%i %i %s", anchor_counter, nrow(thisFCSobject), basename(ananchor)),q=F);
  }
  #T1 <- Sys.time();
  #print(T1-T0);
  return(minCount);
}


# get_cols_to_norm
# List all columns in the most recently created fcs file.
# This is only used if no channelsToAdjust file is specified.
# Try to remove non-data channels:   "Time" "Event_length" "Center" "Offset" "Width" "Residual" 
# return cols_to_norm 
get_cols_to_norm <- function(basedir, anchorKeyword=c()){
  cols_to_skip <- c("Time","Event_length","Center","Offset","Width","Residual");
  if(is.null(anchorKeyword)){
    anchorKeyword <- "";
  }
  basedir_escapeSpace <- gsub(pattern=" ", replacement="\\ ", x=basedir, fixed=TRUE);
  grepForAnchor_escapeSpace <- gsub(pattern=" ", replacement="\\ ", x=anchorKeyword, fixed=TRUE);
  
  ls_cmd <- sprintf("ls -1t %s/*%s*.fcs", basedir_escapeSpace, grepForAnchor_escapeSpace);
  anchors_list <- sort(system(ls_cmd, intern=TRUE, ignore.stdout = FALSE, ignore.stderr = TRUE, wait = TRUE));
  
  if(length(anchors_list) == 0){
    stop("Found no FCS files, nothing to adjust.");
  }
  cols_to_norm <- c();
  anchor_file <- anchors_list[1];
  anchor_object <- read.FCS(anchor_file, which.lines=10);
  anchor_data <- exprs(anchor_object);
  chi <- 0;
  chis_to_drop <- c();
  for(channelname in colnames(anchor_data)){
    chi <- chi + 1;
    for(skip in cols_to_skip){
      #if( length(grep(pattern=skip, x=channelname, ignore.case=TRUE)) > 0 ){}
      if( length(grep(pattern=skip, x=channelname, fixed=TRUE)) > 0 ){
        chis_to_drop <- c(chis_to_drop, chi);
        next;
      }
    }
  }
  cols_to_norm <- colnames(anchor_data)[-chis_to_drop];
  print(sprintf("No channels file. Using channels in most recent FCS file: %s", anchor_file), q=F);
  return(cols_to_norm);
}

# getBatchNumFromFilename
# Parse filename, return batch number.
getBatchNumFromFilename <- function(anchorFileName, batchKeyword="Barcode_", anchorKeyword="anchor stim"){
  underscore_anchorKeyword <- sprintf("_%s", anchorKeyword);
  first_part <- unlist(strsplit(anchorFileName, underscore_anchorKeyword, fixed=TRUE));
  next_part <- unlist(strsplit(first_part[1], batchKeyword, fixed=TRUE));
  this_batch_num <- as.numeric(next_part[2]); # numeric
  return(this_batch_num);
}

# getNZ
# Return non-zero elements of vector.
getNZ <- function(vec){
  wnz <- which(vec > 0);
  return(vec[wnz]);
}

# parseP
# Drop the trailing p, return a number 1-100: "80p" -> 80 ; "50p" -> 50
parseP <- function(str="80p"){
  if(length(grep("p$", str)) == 1){
    num <- as.numeric(sub("p$", "", str));
    if(is.finite(num) && num>=1 && num <=100){
      return(num);
    }
  }
  print(sprintf("parseP: %s doesn't fit pattern", str));
  stop();
}


# getValueMappings
# Return mappingFunctionsList for quantile normalization.
# mappingFunctionsList[[batch]][[acol]]
getValueMappings <- function(anchorKeyword, batchKeyword, basedir, minCount, batches, cols_to_norm, transformation=TRUE, outputfile, nz_only=FALSE){
  
  mt0 <- Sys.time();
  whichlines <- minCount;
  basedir_escapeSpace <- gsub(pattern=" ", replacement="\\ ", x=basedir, fixed=TRUE);
  grepForAnchor_escapeSpace <- gsub(pattern=" ", replacement="\\ ", x=anchorKeyword, fixed=TRUE);
  ls_cmd <- sprintf("ls -1 %s/*%s*.fcs", basedir_escapeSpace, grepForAnchor_escapeSpace);
  anchors_list <- sort(system(ls_cmd, intern=TRUE, ignore.stdout = FALSE, ignore.stderr = TRUE, wait = TRUE));
  anchor_counter <- 0;
  
  # Build a 2-level list: anchorDataListList[[col_to_norm]][[batch]]
  # Top level is indexed by column name (channel) to be adjusted: anchorDataList <- anchorDataListList[[col_to_norm]]
  # Second level is indexed by batch: anchorData_batchX <- unlist(anchorDataList[[batch]])
  # Initialize list.
  anchorDataListList <- list();
  pZeros <- list(); # percentage of zeros.
  # Read events for an anchor.
  # Load into lists
  
  Nbatches <- length(batches);
  for(ananchor in anchors_list){
    thisBatchNum <- getBatchNumFromFilename(basename(ananchor), batchKeyword=batchKeyword, anchorKeyword=anchorKeyword); # note this is numeric
    anchor_counter <- anchor_counter + 1;
    logToFile(outputfile, sprintf("Start loading events batch: %i (%i of %i)", thisBatchNum, anchor_counter, Nbatches), timestamp=TRUE);
    thisFCSobject <- read.FCS(ananchor, transformation=NULL, which.lines=whichlines);
    this_data <- exprs(thisFCSobject);
    for(acol in colnames(this_data)){
      if(!(acol %in% cols_to_norm)){
        next;
      }
      if(transformation){
        anchorDataListList[[acol]][[as.character(thisBatchNum)]] <- asinh(this_data[,acol] * g_asinh_b);
      } else{
        anchorDataListList[[acol]][[as.character(thisBatchNum)]] <- this_data[,acol];
      }
      pZeros[[acol]][[as.character(thisBatchNum)]] <- 100 * sum(anchorDataListList[[acol]][[as.character(thisBatchNum)]] <= 0) /
        length(anchorDataListList[[acol]][[as.character(thisBatchNum)]]);
      
      if(nz_only){
        wgz <- which(anchorDataListList[[acol]][[as.character(thisBatchNum)]] > 0);
        anchorDataListList[[acol]][[as.character(thisBatchNum)]] <- anchorDataListList[[acol]][[as.character(thisBatchNum)]][wgz];
      }
    }
  }
  
  
  # Create the reference quantile.
  nqpoints <- 100000;
  qtype <- 8; # Quantile algorithm. type=7 is default. type=8 recommended by Hyndman and Fan (1996)
  refq <- list();
  for(acol in cols_to_norm){
    # WAS: Data for all batches in this channel.
    #refq[[acol]] <- quantile(unlist(anchorDataListList[[acol]], use.names=FALSE), probs=seq(0,1,length.out=nqpoints), names=FALSE, type=qtype);
    # Now, just use batch 1 as target.
    refq[[acol]] <- quantile(anchorDataListList[[acol]][["1"]], probs=seq(0,1,length.out=nqpoints), names=FALSE, type=qtype);
  }
  
  # Create a mapping function for each batch for each column to norm.
  # This will map values in for the anchor to the reference values.
  mappingFunctionsList <- list();
  # mappingFunctionsList[[batch]][[acol]]
  for(abatch in batches){
    thisBatchChar <- as.character(abatch);
    thisBatchFunctionsList <- list(); # indexed by column name to adjust
    for(acol in cols_to_norm){
      ## Option to help at end?
      #maxRefValue <- max(refq[[acol]]);
      #if(max(anchorDataListList[[acol]][[thisBatchChar]]) < maxRefValue){
      #   qx <- quantile(c(anchorDataListList[[acol]][[thisBatchChar]], maxRefValue), probs=seq(0,1,length.out=nqpoints), names=FALSE, type=qtype);
      #} else{
      #   qx <- quantile(anchorDataListList[[acol]][[thisBatchChar]], probs=seq(0,1,length.out=nqpoints), names=FALSE, type=qtype);
      #}
      qx <- quantile(anchorDataListList[[acol]][[thisBatchChar]], probs=seq(0,1,length.out=nqpoints), names=FALSE, type=qtype);
      spf <- splinefun(x=qx, y=refq[[acol]], method="monoH.FC", ties=min);
      #spf <- approxfun(x=qx, y=refq[[acol]], rule=2, ties=min, yleft=0); # clips
      thisBatchFunctionsList[[acol]] <- spf;
    }
    mappingFunctionsList[[thisBatchChar]] <- thisBatchFunctionsList;
    logToFile(outputfile, sprintf("Done mappingFunctionsList[[%s]]", thisBatchChar), timestamp=TRUE, echo=FALSE);
  }
  
  save(mappingFunctionsList, file=sprintf("%s/mappingFunctionsList.Rdata", dirname(outputfile)));
  
  mt1 <- Sys.time();
  logToFile(outputfile, "getValueMappings duration:", timestamp=TRUE, echo=FALSE);
  logToFile(outputfile, format(mt1-mt0), timestamp=FALSE, echo=FALSE);
  return(mappingFunctionsList);
}

# Ratios. Flip and negate values less than one. 
# Just for visualization.
flipLT1 <- function(vec){
  wlt1 <- which(vec < 1);
  if(length(wlt1) > 0){
    flipped <- -1 * (1/vec[wlt1]);
    vec[wlt1] <- flipped;
  }
  return(vec);
}

# Plot scaling factors, one plot per channel, each with a bar for each batch.
barplot_scalingFactors <- function(scalingFactorsList, postdir){
  # scalingFactorsList[[batch]][[acol]] 
  
  # Format for plotting: index by channel
  # listByCh[[ch]][[batch]] 
  listByCh <- list();
  for(bname in names(scalingFactorsList)){
    if(bname == "1"){
      next;
    }
    sfThisBatchByCh <- scalingFactorsList[[bname]];
    for(ch in names(sfThisBatchByCh)){
      if(is.null(listByCh[[ch]])){
        listByCh[[ch]] <- list();
      }
      listByCh[[ch]][[bname]] <- sfThisBatchByCh[[ch]];
    }
  }
  # plot for each channel
  Nplots <- length(listByCh);
  plotCols <- 9;
  plotRows <- ceiling( Nplots / plotCols )
  pngwidth<-4500; pngheight<-2000;
  pngwidth<-plotCols*500; pngheight<-plotRows*400;
  png(filename=sprintf("%s/ScalingFactors.png", postdir) , width=pngwidth, height=pngheight);
  #png(filename=sprintf("%s/ScalingFactorsFlipped.png", postdir) , width=pngwidth, height=pngheight);
  layout(matrix(1:(plotCols*plotRows), ncol=plotCols, byrow=T));
  par(mar=c(2.5,2.5,2.5,1)+0.1); # 'c(bottom, left, top, right)' default is 'c(5, 4, 4, 2) + 0.1'.
  
  par(cex=1);
  for(ch in names(listByCh)){
    barplot(unlist(listByCh[[ch]]), main=ch, col="limegreen");
    #barplot(flipLT1(unlist(listByCh[[ch]])), main=ch, col="limegreen");
    abline(h=1, lty=2, col="grey", lwd=2);
    #abline(h=-1, lty=2, col="grey", lwd=2);
    
  }
  dev.off();
}

# getScalingFactors
# Return	scalingFactorsList 
# scalingFactorsList[[batch]][[acol]] 
# method = 80p | hybrid | SD | sd 
getScalingFactors <- function(anchorKeyword, batchKeyword, basedir, minCount, batches, cols_to_norm, transformation=TRUE, nz_only=FALSE, outputfile, method="80p"){
  
  mt0 <- Sys.time();
  #whichlines <- minCount;
  whichlines <- NULL;
  basedir_escapeSpace <- gsub(pattern=" ", replacement="\\ ", x=basedir, fixed=TRUE);
  grepForAnchor_escapeSpace <- gsub(pattern=" ", replacement="\\ ", x=anchorKeyword, fixed=TRUE);
  ls_cmd <- sprintf("ls -1 %s/*%s*.fcs", basedir_escapeSpace, grepForAnchor_escapeSpace);
  anchors_list <- sort(system(ls_cmd, intern=TRUE, ignore.stdout = FALSE, ignore.stderr = TRUE, wait = TRUE));
  anchor_counter <- 0;
  
  # Build a 2-level list: anchorDataListList[[col_to_norm]][[batch]]
  # Top level is indexed by column name (channel) to be adjusted: anchorDataList <- anchorDataListList[[col_to_norm]]
  # Second level is indexed by batch: anchorData_batchX <- unlist(anchorDataList[[batch]])
  # Initialize list.
  anchorDataListList <- list();
  pZeros <- list(); # percentage of zeros.
  # Read events for an anchor. 
  # Load into lists
  Nbatches <- length(batches);
  for(ananchor in anchors_list){
    thisBatchNum <- getBatchNumFromFilename(basename(ananchor), batchKeyword=batchKeyword, anchorKeyword=anchorKeyword); # note this is numeric
    anchor_counter <- anchor_counter + 1;
    logToFile(outputfile, sprintf("Start loading events batch: %i (%i of %i)", thisBatchNum, anchor_counter, Nbatches), timestamp=TRUE);
    thisFCSobject <- read.FCS(ananchor, transformation=NULL, which.lines=whichlines);
    this_data <- exprs(thisFCSobject);
    if(!is.null(minCount)){
      if(minCount > nrow(this_data)){
        stop("minCount too high");
      }
      this_data <- this_data[sample.int(n=nrow(this_data), size=minCount),];
    }
    for(acol in colnames(this_data)){
      if(!(acol %in% cols_to_norm)){
        next;
      }
      if(transformation){
        #anchorDataListList[[acol]][[as.character(thisBatchNum)]] <- sort(asinh(this_data[,acol] * g_asinh_b));
        anchorDataListList[[acol]][[as.character(thisBatchNum)]] <- (asinh(this_data[,acol] * g_asinh_b));
      } else{
        #anchorDataListList[[acol]][[as.character(thisBatchNum)]] <- sort(this_data[,acol]);
        anchorDataListList[[acol]][[as.character(thisBatchNum)]] <- (this_data[,acol]);
      }
      #nZeros[[acol]][[as.character(thisBatchNum)]] <- sum(anchorDataListList[[acol]][[as.character(thisBatchNum)]] == 0);
      pZeros[[acol]][[as.character(thisBatchNum)]] <- 100 * sum(anchorDataListList[[acol]][[as.character(thisBatchNum)]] <= 0) / 
        length(anchorDataListList[[acol]][[as.character(thisBatchNum)]]);
    }
  }
  
  # Create the reference.
  #allSD <- list();
  #allMean <- list();
  #allMedian <- list();
  maxFrac0ForMedianThreshold <- 0.39; # hybrid
  if(method == "SD" || method == "sd"){
    maxFrac0ForMedianThreshold <- -1; # do SD only
  }
  else { # method == "80p" || method == "50p" || method == "hybrid" ...
    # For % or median scaling.
    if(method == "hybrid"){
      perc <- .8; # percentile. median would be .5. Highest % zeros in any sample is 76.35%
    } else{ # method == "50p", or "80p"...
      maxFrac0ForMedianThreshold <- 1;   # do 80p only (or median)
      perc <- parseP(method)/100;
      #print(sprintf("method %s using perc: %.2f", method, perc), q=F);
    }
  }
  FractionZerosPooled <- list();
  for(acol in cols_to_norm){
    zmax <- max(unlist(pZeros[[acol]]));
    zmin <- min(unlist(pZeros[[acol]]));
    zsd <- sd(unlist(pZeros[[acol]]));
    pooled_vec <- unlist(anchorDataListList[[acol]], use.names=FALSE);
    #allSD[[acol]] <- sd(pooled_vec);
    #allMean[[acol]] <- mean(pooled_vec);
    #allMedian[[acol]] <- median(pooled_vec);
    zeros_ref <- sum(pooled_vec <= 0);
    length_ref <- length(pooled_vec);
    pzeros_ref <- 100*zeros_ref/length_ref;
    fzeros_ref <- zeros_ref/length_ref;
    FractionZerosPooled[[acol]] <- fzeros_ref;
    #print(sprintf("acol: %s  zsd: %.2f   zmin: %2.f%%  zmax: %.2f%%   pooled0s: %.2f%%   %i of %i", acol, zsd, zmin, zmax, pzeros_ref, zeros_ref, length_ref));
    if( FractionZerosPooled[[acol]] > maxFrac0ForMedianThreshold ){
      #print(sprintf("acol: %s using SD scaling zsd: %.2f   zmin: %2.f%%  zmax: %.2f%%   pooled0s: %.2f%%   %i of %i", acol, zsd, zmin, zmax, pzeros_ref, zeros_ref, length_ref));
    } else{
      #print(sprintf("acol: %s using percentile scaling zsd: %.2f   zmin: %2.f%%  zmax: %.2f%%   pooled0s: %.2f%%   %i of %i", acol, zsd, zmin, zmax, pzeros_ref, zeros_ref, length_ref));
    }
    
  }
  
  
  # Create a scaling factor for each batch for each column to norm.
  # This will scale values in for the anchor to the reference values.
  scalingFactorsList <- list();
  # scalingFactorsList[[batch]][[acol]] 
  for(abatch in batches){
    thisBatchChar <- as.character(abatch);
    thisBatchScalingFactors <- list(); # indexed by column name to adjust
    
    for(acol in cols_to_norm){
      thisvec <- anchorDataListList[[acol]][[thisBatchChar]];
      refvec <- anchorDataListList[[acol]][["1"]];
      if(nz_only){
        thisvec <- getNZ(thisvec);
        refvec <- getNZ(refvec);
      }
      # Make the target based on the fraction of zeros in this channel.
      if( FractionZerosPooled[[acol]] > maxFrac0ForMedianThreshold ){
        # use SD scaling
        scf <- sd(refvec) / sd(thisvec);
      } else{
        # use percentile 80 or median scaling
        refvalue <- get80thPercentile(refvec, perc);
        thisvalue <- get80thPercentile(thisvec, perc);
        if(is.null(refvalue) || refvalue == 0){
          zmaxplus <- ceiling(max(unlist(pZeros[[acol]])) + .01);
          print(sprintf("***  Try increasing percentile scaling to %ip, or exclude channel %s.  ***", zmaxplus, acol), q=F);
          stop(sprintf("zero scaling factor: batch %s channel %s", thisBatchChar, acol));
        }
        if(is.null(thisvalue) || thisvalue == 0){
          zmaxplus <- ceiling(max(unlist(pZeros[[acol]])) + .01);
          print(sprintf("***  Try increasing percentile scaling to %ip, or exclude channel %s.  ***", zmaxplus, acol), q=F);
          stop(sprintf("undefined scaling factor: batch %s channel %s", thisBatchChar, acol));
        }
        scf <- refvalue / thisvalue;
      }
      
      thisBatchScalingFactors[[acol]] <- scf;
    }
    scalingFactorsList[[thisBatchChar]] <- thisBatchScalingFactors;
    logToFile(outputfile, sprintf("Done scalingFactorsList[[%s]]", thisBatchChar), timestamp=TRUE, echo=FALSE);
  }
  
  save(scalingFactorsList, file=sprintf("%s/scalingFactorsList.Rdata", dirname(outputfile)));
  
  barplot_scalingFactors(scalingFactorsList, postdir=dirname(outputfile));
  mt1 <- Sys.time();
  logToFile(outputfile, "getScalingFactors duration:", timestamp=TRUE, echo=FALSE);
  logToFile(outputfile, format(mt1-mt0), timestamp=FALSE, echo=FALSE);
  return(scalingFactorsList); 
}

# BatchAdjust
# Main function for batch adjustment.
# method = 95p | hybrid | SD | sd | quantile | QN
# addExt: Add an extension to the output file name to distinguish from original. eg addExt="_BN"
#         The default addExt=c() makes no change to the output filename.
BatchAdjust <- function(
  basedir=".",
  outdir=".",
  channelsFile = "ChannelsToAdjust.txt",
  batchKeyword="Barcode_",
  anchorKeyword = "anchor stim",
  nz_only=FALSE,
  method="95p",
  transformation=FALSE,
  addExt=c(),
  plotDiagnostics=TRUE){
  
  whichlines <- NULL;
  timestamp <- format(Sys.time(), format="%Y.%m.%d.%H%M%S");
  # Escape any spaces in file names.
  # Also enforce that batchKeyword doesn't contain spaces...
  basedir_escapeSpace <- gsub(pattern=" ", replacement="\\ ", x=basedir, fixed=TRUE);
  outdir_escapeSpace <- gsub(pattern=" ", replacement="\\ ", x=outdir, fixed=TRUE);
  anchorKeyword_escapeSpace <- gsub(pattern=" ", replacement="\\ ", x=anchorKeyword, fixed=TRUE);
  mkdir_cmd <- sprintf(paste("mkdir ", "-p ", "%s"), outdir_escapeSpace);
  mkret <- system(mkdir_cmd, intern=TRUE, ignore.stdout = FALSE, ignore.stderr = TRUE, wait = TRUE);
  outputfile <- sprintf("%s/LOG_BatchAdjust.%s.txt", outdir, timestamp);
  
  logToFile(logfilename=outputfile, logmessage="BatchAdjust.R", timestamp=TRUE, echo=TRUE, overwrite=FALSE);
  logToFile(outputfile, sprintf("basedir:%s",basedir));
  logToFile(outputfile, sprintf("outdir:%s",outdir));
  logToFile(outputfile, sprintf("channelsFile:%s",channelsFile));
  logToFile(outputfile, sprintf("batchKeyword:%s",batchKeyword));
  logToFile(outputfile, sprintf("anchorKeyword:%s",anchorKeyword));
  if(transformation){
    logToFile(outputfile, sprintf("transformation:%s","TRUE"));
  }else{
    logToFile(outputfile, sprintf("transformation:%s","FALSE"));
  }
  if(nz_only){
    logToFile(outputfile, sprintf("nz_only:%s","TRUE"));
  }else{
    logToFile(outputfile, sprintf("nz_only:%s","FALSE"), echo=FALSE);
  }
  logToFile(outputfile, sprintf("method:%s", method));
  
  
  t0 <- Sys.time(); # Time full duration.
  
  
  # Parameter checking.
  if( (is.null(addExt)||addExt=="") && (basedir==outdir)){
    print(basedir);
    print(outdir);
    stop("basedir=outdir will overwrite source files. Please specifiy basedir and/or outdir, or set addExt.");
  }
  # Enforce that batchKeyword doesn't contain spaces...
  if( length(grep(pattern=" ", x=batchKeyword, fixed=TRUE)) > 0 ){
    stop("batchKeyword (%s) must not contain spaces.");
  }
  if(is.null(channelsFile)){
    cols_to_norm <- get_cols_to_norm(basedir=basedir, anchorKeyword=anchorKeyword);
    logToFile(outputfile, "No channelsFile provided.");
  } else{
    logToFile(outputfile, sprintf("Reading channels file: %s", channelsFile));
    cols_to_norm <- unique(as.character(read.table(file=channelsFile, header=F, sep="\n", quote="", as.is=TRUE)[,1]));
  }
  N_ref_channels <- length(cols_to_norm);
  logToFile(outputfile, sprintf("Adjusting %i channels.", N_ref_channels));
  
  batchesToAdjust <- listBatchesPresent(basedir, batchKeyword=batchKeyword, anchorKeyword=anchorKeyword);
  if(! ( 1 %in% batchesToAdjust) ){
    print("Batch 1 anchor not found.", q=F);
    print("Check file naming requirements?", q=F);
    stop("The reference anchor must exist and have the batch number 1.");
  }
  logToFile(outputfile, "batchesToAdjust:");
  logToFile(outputfile, paste(batchesToAdjust, collapse=" "));
  
  
  batch_counter <- 0;
  file_counter <- 0;
  
  nbatches <- length(batchesToAdjust);
  
  minCountAnchors <-  NULL; # get all events. Doesn't need to be the same number for all.
  #minCountAnchors <-  310085; # for anchor stim files
  #if(is.null(minCountAnchors)){
  #   logToFile(outputfile, "minCount: Using all anchor events.");
  #} else{
  #   logToFile(outputfile, sprintf("minCount: Using %i events per anchor.", minCountAnchors));
  #}
  
  if(method == "quantile" || method == "QN" || method == "Quantile"){
    mappingFunctionsList <- getValueMappings(anchorKeyword=anchorKeyword, batchKeyword=batchKeyword, basedir=basedir, minCount=minCountAnchors, batches=batchesToAdjust, cols_to_norm=cols_to_norm, transformation=transformation, outputfile=outputfile, nz_only=nz_only);
    # mappingFunctionsList[[batch]][[colname]]
  } else{
    scalingFactorsList <- getScalingFactors(anchorKeyword=anchorKeyword, batchKeyword=batchKeyword, basedir=basedir, minCount=minCountAnchors, batches=batchesToAdjust, cols_to_norm=cols_to_norm, transformation=transformation, nz_only=nz_only, outputfile=outputfile, method=method);
    # scalingFactorsList[[batch]][[colname]] 
  }
  
  print("Adjusting...");
  # For each batch...
  for(thisbatch in batchesToAdjust){
    td0 <- Sys.time(); # time this batch
    if(method == "quantile" || method == "QN"){
      #thisBatchMappingFunctions <- mappingFunctionsList[[as.character(thisbatch)]];
      thisBatchAdjustments <- mappingFunctionsList[[as.character(thisbatch)]];
    } else{
      #thisBatchScalingFactors <- scalingFactorsList[[as.character(thisbatch)]];
      thisBatchAdjustments <- scalingFactorsList[[as.character(thisbatch)]];
    }
    
    batch_counter <- batch_counter + 1;
    logToFile(outputfile, sprintf("Batch %i of %i", batch_counter, nbatches));
    logToFile(outputfile, as.character(thisbatch));
    
    # Apply to each file.
    #  Non-Anchor filenaming:   xxx[batchKeyword][##]_xxx.fcs
    ls_cmd <- sprintf("ls -1 %s/*%s%i_*.fcs", basedir_escapeSpace, batchKeyword, thisbatch);
    fcs_files <- sort(system(ls_cmd, intern=TRUE, ignore.stdout = FALSE, ignore.stderr = TRUE, wait = TRUE));
    for(fcsfile in fcs_files){
      file_counter <- file_counter + 1;
      logToFile(outputfile, sprintf("file %i", file_counter));
      
      # Make the adjustment per channel.
      tf0 <- Sys.time();
      thisFCSobject <- read.FCS(fcsfile, transformation=NULL, which.lines=whichlines);
      this_data <- exprs(thisFCSobject);
      these_parameters <- parameters(thisFCSobject);
      
      if(transformation){
        this_data <- asinh(this_data*g_asinh_b);
      }
      for(colname in colnames(this_data)){
        if(!(colname %in% cols_to_norm)){
          next;
        }
        if(method == "quantile" || method == "QN"){
          # Apply mapping function.
          mapFun <- thisBatchAdjustments[[colname]];
          this_data[,colname] <- mapFun(this_data[,colname]);
        } else{
          # Apply scaling factor.
          scalingFactor <- thisBatchAdjustments[[colname]];
          this_data[,colname] <- scalingFactor * (this_data[,colname]);
        }
        w0 <- which(this_data[,colname] < 0);
        if(length(w0) > 0){
          this_data[w0,colname] <- 0;
        }
      }
      
      if(transformation){ # undo the transform: inverse of asinh is sinh
        # asinh(ck/5) == ckt
        # 5*sinh(ckt) == ck
        this_data <- sinh(this_data) / g_asinh_b;
      }
      #exprs(thisFCSobject) <- this_data;
      # Need to set the ranges:
      new_params <- pData(these_parameters);
      new_params[,"minRange"] <- floor(apply(this_data, MARGIN=2, FUN=min));
      new_params[,"maxRange"] <- ceiling(apply(this_data, MARGIN=2, FUN=max));
      new_params[,"range"] <- new_params[,"maxRange"] - new_params[,"minRange"] + 1;
      pData(these_parameters) <- new_params;
      desc <-  thisFCSobject@description;
      for(paramrowi in 1:nrow(new_params)){
        rangeStr <- sprintf("$P%iR", paramrowi);
        desc[[rangeStr]] <- new_params[paramrowi, "maxRange"] + 1;
      }
      newFCSobject <- flowFrame(exprs=this_data, parameters=these_parameters, description=desc);
      #newFCSobject <- flowFrame(exprs=this_data, parameters=these_parameters)
      
      #addExt <- "_BN";
      # Add an extension to the output file name to distinguish.
      #  addExt <- c(); # or don't
      if(is.null(addExt)){
        outfilename <- sprintf("%s/%s", outdir, basename(fcsfile));
      }else{
        replfcs <- sprintf("%s.fcs", addExt);
        #basenameW_BNext <- gsub(pattern="\\.fcs$", replacement="_BN.fcs", x=basename(fcsfile), fixed=F);
        basenameW_BNext <- gsub(pattern="\\.fcs$", replacement=replfcs, x=basename(fcsfile), fixed=F);
        outfilename <- sprintf("%s/%s", outdir, basenameW_BNext);
      }
      #write.FCS(thisFCSobject, filename=outfilename);
      write.FCS(newFCSobject, filename=outfilename);
      tf1 <- Sys.time();
      #print("Per file read,norm,write:");
      #print(tf1-tf0);
    }
    td1 <- Sys.time();
    logToFile(outputfile, "Per batch:", echo=FALSE); # ~1-2 minutes 
    logToFile(outputfile, format(td1-td0), echo=FALSE);
  } # for each batch number
  
  t1 <- Sys.time();
  logToFile(outputfile, "Finished. Duration:", timestamp=TRUE);
  logToFile(outputfile, format(t1-t0)); 
  
  if(plotDiagnostics){
    
    t00 <- Sys.time();
    print("Testing total variance reduction:", q=F);
    totalVar_allEvents(predir=basedir, postdir=outdir, batchKeyword=batchKeyword, channelsFile=channelsFile, anchorKeyword=anchorKeyword, colorPre="blue", colorPost="wheat");
    print("Plotting pre/post distributions...", q=F);
    call_plotAllAPrePost1ch(plotnz=TRUE, xlim=c(0,8), postdir=outdir, anchorKeyword=anchorKeyword, batchKeyword=batchKeyword, predir=basedir, colorPre="blue", colorPost="wheat", addExt=addExt, channelsFile=channelsFile);
    t11 <- Sys.time();
    logToFile(outputfile, "Finished diagnostic plots, duration:", timestamp=TRUE);
    logToFile(outputfile, format(t11-t00)); 
    
  }
} # BatchAdjust









#####################################################
#  Plot distributions pre/post batch adjustment.
# 

# Plot all anchors for one channel.
plotAllAPrePost1ch <- function(ch="CD3", xlim=c(0,8), plotnz=TRUE, postdir=c(), anchorKeyword="anchor stim", batchKeyword="Barcode_", predir=c(), colorPre="blue", colorPost="wheat", addExt=c()){
  grepForAnchor_escapeSpace <- gsub(pattern=" ", replacement="\\ ", x=anchorKeyword, fixed=TRUE);
  predir_escapeSpace <- gsub(pattern=" ", replacement="\\ ", x=predir, fixed=TRUE);
  
  chname <- get_ch_name(ch, predir);
  print(sprintf("%s", chname), q=F);
  anchorKeyword_underscore <- gsub(pattern=" ", replacement="_", x=anchorKeyword, fixed=TRUE);
  
  basedir <- postdir;
  basedir_escapeSpace <- gsub(pattern=" ", replacement="\\ ", x=basedir, fixed=TRUE);
  pngoutdir_escaped <- sprintf("%s/DistributionPlots", basedir_escapeSpace);
  mkdir_cmd <- sprintf(paste("mkdir ", "-p ", "%s"), pngoutdir_escaped);
  mkret <- system(mkdir_cmd, intern=TRUE, ignore.stdout = FALSE, ignore.stderr = TRUE, wait = TRUE);
  pngoutdir <- sprintf("%s/DistributionPlots", basedir);
  pngname <- sprintf("%s/%s.png", pngoutdir, chname);
  
  ls_cmd <- sprintf("ls -1 %s/*%s*.fcs", basedir_escapeSpace, grepForAnchor_escapeSpace);
  anchors_list <- system(ls_cmd, intern=TRUE, ignore.stdout = FALSE, ignore.stderr = TRUE, wait = TRUE);
  anchors_list <- sortAnchorsByBatch(anchors_list,  anchorKeyword=anchorKeyword, batchKeyword=batchKeyword);
  N_anchors <- length(anchors_list);
  
  pngwidth<-900; 
  pngheight <- N_anchors * 260;
  png(filename=pngname, width=pngwidth, height=pngheight);
  layout(matrix(1:(2*N_anchors), ncol=1));
  
  for(basedir in c(predir, postdir)){
    for(fname in anchors_list){
      if(basedir == predir){
        col <- colorPre;
        if(is.null(addExt)){
          #fname <- sprintf("%s/%s", predir_escapeSpace, basename(fname));
          fname <- sprintf("%s/%s", predir, basename(fname));
        } else{
          addedPat <- sprintf("%s.fcs$", addExt);
          #fname <- sprintf("%s/%s", predir_escapeSpace, gsub(pattern=addedPat, replacement=".fcs", x=basename(fname), fixed=F));
          fname <- sprintf("%s/%s", predir, gsub(pattern=addedPat, replacement=".fcs", x=basename(fname), fixed=F));
        }
      } else{
        col <- colorPost;
      }
      batchNum <- getBatchNumFromFilename(fname, batchKeyword=batchKeyword, anchorKeyword=anchorKeyword);
      fce <- fc_read(fname, trans=TRUE, fullnames=TRUE);
      plot1(fce, ch=ch, xlim=xlim, plotnz=plotnz, chname=chname, barcolor=col, basedir=basedir);
      #title(ylab=batchNum, las=1, cex.lab=3, line=0, col.lab=1);
      mtext(text=batchNum, side=2, cex=2, line=0, las=1);
    }
  }
  dev.off();
}

get_ch_name <- function(ch, basedir=c()){
  if(is.null(basedir)){
    stop("get_ch_name must be supplied with basedir, a directory to find a sample fcs file.")
  }
  basedir_escapeSpace <- gsub(pattern=" ", replacement="\\ ", x=basedir, fixed=TRUE);
  ls_cmd <- sprintf("ls -1t %s/*.fcs", basedir_escapeSpace);
  anchors_list <- system(ls_cmd, intern=TRUE, ignore.stdout = FALSE, ignore.stderr = TRUE, wait = TRUE);   
  anchor_file <- anchors_list[1];
  fc <- read.FCS(anchor_file, which.lines=10);
  chname <- grep(ch,pData(parameters(fc))$desc, value=T);
  if(length(chname) == 0){
    chi <- grep(ch,colnames(fc), value=F);
    chname <- pData(parameters(fc))$desc[chi];
  }
  if(length(chname) > 1){
    chname <- chname[1];
  }
  return(chname);
}


call_plotAllAPrePost1ch <- function(plotnz=TRUE, xlim=c(0,8), postdir=c(), anchorKeyword="anchor stim", batchKeyword="Barcode_", predir=c(), colorPre="lightblue", colorPost="wheat", addExt=c(), channelsFile = "ChannelsToAdjust_example.txt"){
  #print(sprintf("call_plotAllAPrePost1ch %s", postdir), q=F);
  cols_to_norm <- unique(as.character(read.table(file=channelsFile, header=F, sep="\n", quote="", as.is=TRUE)[,1]));
  for(ch in cols_to_norm){
    plotAllAPrePost1ch(ch=ch, xlim=xlim, plotnz=plotnz, postdir=postdir, batchKeyword=batchKeyword, anchorKeyword=anchorKeyword, predir=predir, colorPre=colorPre, colorPost=colorPost, addExt=addExt);
  }
}


# fc_read() Load data from .fcs file.
# fullnames=TRUE -> 145Nd_IFNg eg  more informative
# fullnames=FALSE -> Nd145Di eg  matches cols_to_norm
# trans=TRUE -> asinh
fc_read <- function(fname="011118_Barcode_7_anchor stim_CD3+ CD19+.fcs", trans=TRUE, which.lines=NULL, fullnames=TRUE){
  fc <- read.FCS(fname, transformation=NULL, which.lines=which.lines, truncate_max_range=FALSE); # flowFrame
  fce <- exprs(fc); # matrix
  # flowFrame is an AnnotatedDataFrame from Biobase. To access useful names, have to use this:
  #p <- parameters(fc);
  #pData(p)$desc
  if(fullnames){
    colnames(fce) <- pData(parameters(fc))$desc;
  }
  if(trans){
    fce <- asinh(fce*g_asinh_b);
  }
  return(fce);
}


# plot1, fce is loaded by fc_read via flowCore. Channel names (ch and chname) are affected by fullnames=T|F
# chname is the column name to plot. 
plot1 <- function(fce, ch="In113Di", xlim=c(0,8), plotnz=FALSE, chname="113In_CD57", bordercolor=1, barcolor="wheat", brks=c(), basedir=c()){
  if(is.null(chname)){
    chname <- get_ch_name(ch, basedir);
  }
  if(plotnz){
    wnz <- which(fce[,chname] > 0);
    vec <- fce[wnz,chname];
  } else{
    vec <- fce[,chname];
  }
  #barcolor <- "lightgreen";
  #barcolor <- "wheat";
  #barcolor <- "violet";
  #barcolor <- "lightblue";
  if(is.null(brks)){
    #brks <- 300;
    #brks <-  seq(0,9,length.out=400);
    #brks <-  seq(0,9,length.out=300);
    #brks <-  seq(0,10,length.out=300);
    #brks <-  seq(0,10,length.out=250);
    #brks <-  seq(0,11,length.out=250);
    brks <-  seq(0,12.5,length.out=250);
  }
  if(is.null(xlim)){
    hist(vec, breaks=brks, xlab="", main="", border=bordercolor, col=barcolor, yaxt="n");
  } else{
    hist(vec, breaks=brks, xlab="", ylab="", main="", xlim=xlim, border=bordercolor, col=barcolor, yaxt="n", xaxt="n");
  }
  axis(side=1, labels=FALSE);
}


sortAnchorsByBatch <- function(anchor_list,  anchorKeyword="anchor stim", batchKeyword="Barcode_"){
  Nanchors <- length(anchor_list);
  batchNums <- rep(0, Nanchors);
  for(ai in 1:Nanchors){
    fname <- anchor_list[ai];
    batchNums[ai] <- getBatchNumFromFilename(fname, batchKeyword=batchKeyword, anchorKeyword=anchorKeyword);
  }
  oi <- order(as.numeric(batchNums));
  return(anchor_list[oi]);
}



############################################################################################
# Testing mean cytokine levels, in all cell events,  not in subpopulations.


# Return a vector of channel names, long names from the description field desc.
get_cols_to_norm_long_names <- function(channelsFile = "ChannelsToAdjust_example.txt", basedir){
  cols_to_norm <- unique(as.character(read.table(file=channelsFile, header=F, sep="\n", quote="", as.is=TRUE)[,1]));
  long_names <- c();
  for(cn in cols_to_norm){
    long_names <- c(long_names, get_ch_name(cn, basedir=basedir));
  }
  return(long_names);
}

# Compute the trace of the covariance matrix.
# Input matricies are  Rows:markers   X   Cols:barcode samples
varDiffTotalPrePostT <- function(tbl_pre, tbl_post){
  
  ####   Transpose   ####
  tbl_pre <- t(tbl_pre);
  tbl_post <- t(tbl_post);
  #######################
  
  cpre <- cov(tbl_pre);
  cpost <- cov(tbl_post);
  epre <- eigen(cpre, symmetric=TRUE, only.values=TRUE);
  epost <- eigen(cpost, symmetric=TRUE, only.values=TRUE);
  trace_pre <- sum(epre$values);
  trace_post <- sum(epost$values);
  return(trace_pre - trace_post);
}

# Swap columns of pre/post: 1to1, 2to2, etc
# Calculate varDiffTotal for all possible permutations. Return vector.
permuteColsPrePost <- function(mat_pre, mat_post){
  if(!(all(colnames(mat_pre) == colnames(mat_post)))){
    stop("pre/post colnames don't match");
  }
  mat <- cbind(mat_pre, mat_post);
  N_bcs <- ncol(mat_pre);
  tfvec <- rep(TRUE, N_bcs);
  count <- 0;
  testres <- rep(0, (2^N_bcs));
  for(itr in 0:((2^N_bcs)-1)){
    tfvec1 <- as.logical(intToBits(itr)[1:N_bcs]);
    tfvec <- c(tfvec1, !tfvec1)
    perm_mat_pre <- mat[,tfvec];
    perm_mat_post <- mat[,!tfvec];
    count <- count + 1;
    testres[count] <- varDiffTotalPrePostT(perm_mat_pre, perm_mat_post);
  }
  return(testres);
}


############################################################################################

# boxplot (top left)
# rows: cytokines, cols: anchor samples
# addPoints: add a point to the boxplot for each sample. (may want to adjust point size)
boxSummaryValsPrePost <- function(mat_pre, mat_post, colorPre="blue", colorPost="wheat", addPoints=FALSE){
  # make a list
  valueslist <- list();
  for(rowi in nrow(mat_pre):1){
    vname <- sprintf("post %s", rownames(mat_post)[rowi]);
    valueslist[[vname]] <- mat_post[rowi,];
    vname <- sprintf("pre %s", rownames(mat_pre)[rowi]);
    valueslist[[vname]] <- mat_pre[rowi,];
  }
  boxplot(valueslist, boxwex=.5, xaxt="n", yaxt="n", cex.axis=1, las=2, range=0, lwd=1, boxlwd=1.25, horizontal=TRUE, col=rep(c(colorPost, colorPre), 2*nrow(mat_pre)), labels=FALSE);
  axis(side=1, las=0); # 1=below
  
  if(addPoints){
    if(nrow(mat_pre) > 15){
      ptcex <- .3;
    } else{
      ptcex <- .8;
    }
    for(vi in 1:length(valueslist)){
      points(y=rep(vi, length(valueslist[[vi]])), x=valueslist[[vi]], cex=ptcex, lwd=1, col="darkred");
      points(y=rep(vi, length(valueslist[[vi]])), x=valueslist[[vi]], cex=ptcex/2, lwd=1, col="yellow");
    }
  }
}

# Variance barplot, top right
varBarPlot <- function(mat_pre, mat_post, colorPre="blue", colorPost="wheat"){
  if(nrow(mat_pre) > 15){
    cexnames <- .7;
  } else{
    cexnames <- 1;
  }
  vpre <- apply(mat_pre, MARGIN=1, FUN=var);
  vpost <- apply(mat_post, MARGIN=1, FUN=var);
  mat <- rbind( rev(vpost), rev(vpre));
  barplot(height=mat, horiz=TRUE, beside=TRUE, col=rep(c(colorPost, colorPre), 2*nrow(mat_pre)), las=1, cex.names=cexnames);
  #title(main="Variance (in mean signal)");
  title(main="Variance");
}

get_bar_code_from_filename <- function(fname, batchKeyword="Plate"){
  parts <- unlist(strsplit( x=gsub(pattern=".fcs$", replacement="", x=fname), split=batchKeyword));
  subparts <- unlist(strsplit( x=parts[2], split="_"));
  return(subparts[1]);
}

# Return a matrix of values (mean or maybe var).
# Rows are cytokines (or any marker).
# Cols are anchor samples.
# Does it matter if trans=FALSE or TRUE?
# Similar to above, but not subpopulations. All cell events.
get_summaries_per_channel <- function(basedir, cols_to_use, batchKeyword="Plate", anchorKeyword = "Sample2"){
  basedir_escapeSpace <- gsub(pattern=" ", replacement="\\ ", x=basedir, fixed=TRUE);
  ls_cmd <- sprintf("ls -1 %s/*.fcs", basedir_escapeSpace);
  fcsfiles <- sort(system(ls_cmd, intern=TRUE, ignore.stdout = FALSE, ignore.stderr = TRUE, wait = TRUE));
  #togrep <- "anchor stim";
  togrep <- anchorKeyword;
  filelist <- sort(grep(togrep, fcsfiles, fixed=TRUE, value=TRUE)); 
  outmat <- c();
  for(fname in filelist){
    #print(fname, q=F);
    bc <- get_bar_code_from_filename(fname, batchKeyword);
    #fce <- fc_read(fname=fname, trans=FALSE, which.lines=NULL, fullnames=TRUE);
    fce <- fc_read(fname=fname, trans=TRUE, which.lines=NULL, fullnames=TRUE);
    tf2use <- colnames(fce) %in% cols_to_use;
    vals <- apply(fce[,tf2use], MARGIN=2, FUN="mean"); # or var
    outmat <- cbind(outmat, vals);
    colnames(outmat)[ncol(outmat)] <- bc;
  }
  return(outmat);
}

# Three plots: 
#  1. horizontal boxplot, pre and post box for each channel.
#  2. Variance of plot 1. 
#  3. Null distribution and p-value
# Total variance (reduction) in mean cytokine levels for all cell events
# Permutation test for significance swapping pre/post.
totalVar_allEvents <- function(predir, postdir, batchKeyword, channelsFile, anchorKeyword, colorPre="blue", colorPost="wheat"){
  t00 <- Sys.time();
  cols_to_use <- get_cols_to_norm_long_names(channelsFile=channelsFile, basedir=predir);
  
  mat_pre <- get_summaries_per_channel(predir, cols_to_use, batchKeyword, anchorKeyword);
  mat_post <- get_summaries_per_channel(postdir, cols_to_use, batchKeyword, anchorKeyword);
  if(!(all(colnames(mat_pre) == colnames(mat_post)))){
    stop("totalVar_allEvents: File names in pre and post directories don't match.");
  }
  test_real <- varDiffTotalPrePostT(mat_pre, mat_post);
  perm_tests <- permuteColsPrePost(mat_pre, mat_post);
  N_as_or_more_extreme <- sum(perm_tests >= test_real);
  pv <- N_as_or_more_extreme / length(perm_tests);
  
  pngwidth<-1600; pngheight<-1200;
  png(filename=sprintf("%s/PrePostVariance.png", postdir) , width=pngwidth, height=pngheight);
  layout(matrix(c(1,2,3,3), nrow=2, byrow=TRUE));
  title <- sprintf("%s", "Mean signal intensity per replicate");
  
  # plot 1:
  par(cex=1.5);
  par(mar=c(2,1,2,0)); # 'c(bottom, left, top, right)' default is 'c(5, 4, 4, 2) + 0.1'.
  boxSummaryValsPrePost(mat_pre, mat_post, colorPre, colorPost);
  title(main=title);
  legend("topright", legend=c("Pre", "Post"), col=c(colorPre, colorPost), lwd=9);
  legend("topright", legend=c("Pre", "Post"), col=c(colorPre, colorPost), lwd=9);
  
  # plot 2: variance bars
  #par(mar=c(2,0,2,1)); # 'c(bottom, left, top, right)' default is 'c(5, 4, 4, 2) + 0.1'.
  par(mar=c(2,7.5,2,1)); # 'c(bottom, left, top, right)' default is 'c(5, 4, 4, 2) + 0.1'.
  varBarPlot(mat_pre, mat_post, colorPre, colorPost);
  
  # plot 3: Bottom NULL dist hist
  par(cex=1.5);
  par(mar=c(5,4.5,4,1)+0.1); # 'c(bottom, left, top, right)' default is 'c(5, 4, 4, 2) + 0.1'.
  #hist(perm_tests, breaks=70, main="", xlab="Test statistic null distribution");
  #hist(perm_tests, breaks=70, main="", xlab="Change in total variance null distribution");
  hist(perm_tests, breaks=70, main="", xlab="Change in total variance", cex.lab=1.5);
  abline(v=test_real, col=2, lwd=5);
  #title(main=sprintf("p = %.05f", pv), line=-1);
  title(main=sprintf("p = %.05f", pv), line=0);
  dev.off();
  print(sprintf("p = %.05f", pv), q=F);
  
  t11 <- Sys.time();
  t11-t00;
}
