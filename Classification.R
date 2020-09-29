#!/export/data/apps/R-3.1.3/bin/Rscript

#get directory of script
opts <- commandArgs(trailingOnly = FALSE)

file.arg.name <- "--file="
script.name <- sub(file.arg.name, "", opts[grep(file.arg.name, opts)])
scriptdir <- dirname(script.name)

m <- match("--args", opts, 0L)
if (m) {
  opts <- opts[-seq_len(m)]
} else {
  opts <- character()
}

print(paste0("options: '", paste(opts, collapse="', '"), "'"))
usage <- "Rscript Classification.R <NKI_1M file> <sample_type> = breast/ovarian; <cls> = b1.191/b1.371/b1/b2 ; <correct_platform> = TRUE/FALSE ; <missing2centroid> = TRUE/FALSE <outcls> = output directory";


if (length(opts) != 6) {
  print(usage); 
} else {
file <- opts[1];
sample_type <- opts[2];
cls <- opts[3];
correct_platform <- as.logical(opts[4]);
missing2centroid <- as.logical(opts[5]);
outcls <- as.character(opts[6])

outbase <- paste0(outcls, "/", sub("\\.txt","", basename(file)))
outbelong <- paste0(cls,'-',sample_type,'-',correct_platform,'-',missing2centroid)

cls_session.txt <- paste0(outbase, ".session.", outbelong,".txt")
if (file.exists(cls_session.txt)) {
  stop(paste("Session file already present, remove", cls_session.txt, "and rerun to overwrite."))
}
sink(file=cls_session.txt)

print(paste0("scriptdir: ",scriptdir))

#require(multicore) // disabled to circumvent setting up outside of server environment. 
# no really heavy lifting done, segmentation could be run in parallel to decrease waiting time 
#(without: ~ 1 min on 1 laptop core)
# segmentation methods
require(cghseg)
# Range methods
require(GenomicRanges)
# interpolation and padding
require(zoo)
# shrunken centroid classification
require(pamr)
# BRCA1 and BRCA2 breast cancer classifiers
require(nkiBRCA)

####################################################################################################
### functions (not in separate file because uses arguments above)
####################################################################################################

# doCghSeg #########################################################################################
####################################################################################################
# adapted from chris_cghsegmentation_functions.R
# -------------------------------------------------------------------
# Copyright 2011 Christiaan Klijn <c.klijn@nki.nl>
# Project: CGH segmentation
# Description: Functions to segment your CGH data
# -------------------------------------------------------------------

doCghSeg <- function(chromNum, allKC, chrom, maploc) {
  tumnames <- colnames(allKC)[-1:-2]
  profiles <- as.matrix(allKC[chrom == chromNum, -c(1,2)])
  subchrom <- chrom[chrom==chromNum]
  submaploc <- maploc[chrom==chromNum]
		
	# check profiles are ordered (errors might be introduced by liftover)
  stopifnot(order(submaploc)==1:length(submaploc)) 

  # Segmentation
  n <- nrow(profiles)
  CGHd <- new("CGHdata",Y=profiles)
  CGHo <- new("CGHoptions")
  alpha(CGHo) <- 0.1 ## look for 1 to a maximum of 500 segments per tumor,
  select(CGHo) <- "mBIC" ## this define the way the number of segment is selected ( in between 1 and Kmax)

  CGHr <- uniseg(CGHd,CGHo) ##The segments with their mean can be found in CGHr@mu.

  allstart <- unlist(lapply(CGHr@mu, function (x) {return(x$begin)}))
  allend <- unlist(lapply(CGHr@mu, function (x) {return(x$end)}))
  allmeans <- unlist(lapply(CGHr@mu, function (x) {return(x$mean)}))

  allID <-rep(gsub('X','', names(CGHr@mu)),
    unlist(lapply(CGHr@mu, nrow))) 
  
  segments <- GRanges(
    seqnames = Rle(paste('chr', chromNum, sep=''), length(allstart)),
    ranges=IRanges(start=submaploc[allstart], end=submaploc[allend]),
    seg.mean=allmeans,
    ID=allID) # Make a range object to contain the start and the end probe for each segment

  return(segments)
}

# extractCghSeg####################################################################################
###################################################################################################
##    Philip Schouten 2013 <philip.schouten@gmail.com>
##    make a chrom, maploc, sample .... sample dataframe from doCghSeg output

extractCghSeg <- function(cghSegObj, chrom, maploc, coln) {

  if( is.null(names(cghSegObj))) {
    stop('cghSegObj does not have chromosome names')
  }

  startInd <- lapply(unique(as.character(chrom)), function(x)  match(start(cghSegObj[[x]]), maploc[chrom==x]))
  names(startInd) <- unique(as.character(chrom))
  endInd <- lapply(unique(as.character(chrom)), function(x) match(end(cghSegObj[[x]]), maploc[chrom==x]))
  names(endInd) <- unique(as.character(chrom))
  segMeans <- lapply(cghSegObj, function(x) values(x)$seg.mean)
  patients <- lapply(cghSegObj, function(x) values(x)$ID)

	# unpack segments
  tmp1 <- lapply(unique(as.character(chrom)), function(x) rep(segMeans[[x]], times=((endInd[[x]]+1)-startInd[[x]])))

	# catch error
  if(!length(unlist(tmp1)) == length(chrom)*length(unique(unlist(patients)))) {
    # duplicated maplocs in nki cause this problem
    # only chrom 3 has this problem:
    # lapply(1:22, function(x) duplicated(maploc[chrom==x]))
    # change nki chrom 3 maploc 249
    stop('missing segMeans for some probes')
  }

	names(tmp1) <- unique(as.character(chrom))
  # convert to matrix
  tmp2 <- lapply(unique(as.character(chrom)), function(x) matrix(tmp1[[x]], ncol= length(unique(patients[[x]]))))

	# reduce to dataframe
  final <- do.call(rbind, tmp2)

	# convert to KCSmart dataframe
  final2 <- data.frame(chrom=chrom, maploc=maploc, final, stringsAsFactors=F)
  colnames(final2) <- coln

  return(final2)
}

# fixMissing2Centroid#################################################################################################
######################################################################################################################
##    Philip Schouten 2013 <philip.schouten@gmail.com>
##    find rows with all values missing and put them to classifier centroid
##    mean. Finds single NAs and replace with linear approximation through the
##    zoo package. Or interpolate all missings with zoo package.

## the input definitions for the various classifiers are:
## Breast BRCA1 191 probe: cls='b1.191', sample_type='breast'
## Breast BRCA1 371 probe: cls='b1.371', sample_type='breast'
## Breast BRCA2 703 probe: cls='b2', sample_type='breast'
## Ovarian: cls=NA, sample_type='ovarian', fillM=NA 
## For ovarian cls is irrelevant because the processing is the same. If probes need to be set to centroid mean an
## ovarian cancer specific BRCA1 and BRCA2 centroid object should be created and the processing should be done as
## for the breast cancer case below. Since the high correlation between copy number segments it was deemed not necessary since
## the classifier is based on segmented data. Filling with na.approx retains information of the surrounding non-missing probes 
## by linear interpolation

fixMissing2Centroid <- function(cls=c('b1.371', 'b1.191','b2',NA), dt, fillM=c('zoo', 'ct', NA), sample_type=c('breast','ovarian')) {
    # first create data.frame; this is equal for both tumor types
            
    # the classifier was trained on hg18 data
    plf <- read.delim(paste0(scriptdir,"/platformnki_hg18.txt"))

    # merge the input data with the platform file. The platform file contains duplicated chrom and maploc positions
    # which correspond to triple spotted array positions which are marked by a different order value. In the merge we use
    # all.x to retain the dimensions of the matrix.
    comb <- merge(plf, dt, by.x=c('chrom','maploc'), by.y=c('chrom', 'maploc'), all.x=T)
    comb <- comb[order(comb$chrom, comb$maploc), ]
    # remove chromosome Y values, not required for classification
    comb[comb$chrom == 24,-1:-8] <- 0        
    # we might have non-mappable positions; all values in a row are NA or random missings in a sample; some values in a row are NA
    allmissing <- apply(comb[, -1:-8], 1, function(x) all(is.na(x)))
    somemissing <- apply(comb[,-1:-8], 1, function(x) any(is.na(x))) & ! allmissing
    pos_some_missing <- which(is.na(as.matrix(comb[somemissing,-1:-8])))

    if (sample_type == 'breast')  {
      # for breast cancer samples probes that contain only missings (nonmappable) should be set to the
      # to mean of the sum of the class centroids to prevent influence on the classification (Schouten et al BCRT 2012)
      # here might also be implemented to interpolate the missing probes. 
          
      # get centroids of the various classifiers, get the number of columns and get the classifier probes
      if (cls=='b1.371') {
        ct <- b1.371.ct[,9:10]
        }
    
      if (cls=='b1.191') {
        ct <- b1.191.ct[,9:10]
        }
    
      if (cls=='b2') {
        ct <- b2.704.ct[,9:10]

        # here segmentation with NA approximated, since the segmentation functions do not work with missings
        # since the probes to be changed afterwards are defined earlier filling the missings here is no problem
        comb[,-1:-8] <- na.approx(comb[,-1:-8], na.rm=F)
        comb[,-1:-8] <- na.locf(comb[,-1:-8], fromLast=T,na.rm=F)
        comb[,-1:-8] <- na.locf(comb[,-1:-8])

        sg <- lapply(1:24, function(x) doCghSeg(x, allKC=comb[, colnames(comb) %in% colnames(dt)],
            chrom=plf$chrom, maploc=plf$maploc))

        names(sg) <- 1:24

        comb[,-1:-8] <- extractCghSeg(sg, chrom=plf$chrom, maploc=plf$maploc, coln=colnames(dt))[,-1:-2]
        }

      classprobe <- ! (ct[1] == ct[2])

      #create a filling object
      fl <- matrix(rep((ct[,1]+ ct[,2])/2,   ncol(comb)-8), ncol=(ncol(comb)-8))

      # fill individual missing probes with classifier centroids in case ct is chosen
      if(fillM=='ct') {
        tm <- as.matrix(comb[,-1:-8])
        tm[allmissing,] <- fl[allmissing,]
        #	 change somemissing to interpolated, they're randomly missing and therefore should not be filled with centroid avg.
        # set the segmented data back to NA for BRCA2 (for BRCA1 does not do anything, some_missing are not yet filled, as
        # is caused by segmentation for the BRCA2 classifier)
        tm[somemissing,][pos_some_missing] <- NA
        comb[,-1:-8] <- tm
        # fill random missings with interpolation
        comb[,-1:-8] <- na.approx(comb[,-1:-8], na.rm=F)
        # push backward and forward to remove missings at the start and end
        comb[,-1:-8] <- na.locf(comb[,-1:-8], fromLast=T,na.rm=F)
        comb[,-1:-8] <- na.locf(comb[,-1:-8])
      }	else if (fillM=='zoo') {					
        # fill with linear interpolation in case 'zoo' is chosen.	
        comb[,-1:-8] <- na.approx(comb[,-1:-8], na.rm=F)
        # push backward and forward to remove missings at start and end
        comb[,-1:-8] <- na.locf(comb[,-1:-8], fromLast=T,na.rm=F)
        comb[,-1:-8] <- na.locf(comb[,-1:-8])
      }
	     
    cat(paste(sum(allmissing & classprobe), ' classifier probes set to centroid for all samples\n', sep=""))
    cat(paste(sum(classprobe & somemissing ), ' classifier probes approximated or put to centroid in some samples \n', sep='')) 
    }	
			
	if (sample_type == 'ovarian') {
    # linear interpolation for ovarian cancer due to high correlation between probes.
    comb[,-1:-8] <- na.approx(comb[,-1:-8], na.rm=F)
    comb[,-1:-8] <- na.locf(comb[,-1:-8], fromLast=T)
    comb[,-1:-8] <- na.locf(comb[,-1:-8])

    # segmentation
    sg <- lapply(1:24, function(x) doCghSeg(x, allKC=comb[, colnames(comb) %in% colnames(dt)],
    chrom=plf$chrom, maploc=plf$maploc))
    names(sg) <- 1:24
    comb[,-1:-8] <- extractCghSeg(sg, chrom=plf$chrom, maploc=plf$maploc, coln=colnames(dt))[,-1:-2]
    } 

  return(comb[, colnames(comb) %in% colnames(dt)])
}

# fillMat ##############################################################################################################
########################################################################################################################
##    Philip Schouten 2020 <philip.schouten@gmail.com>
# fill missing by linear interpolation and subsequently fill trailing and starting NAs by bringing backward and forward
# the last known 

fillMat <- function(dt) {
    # The breast cancer classifier was developed on hg18 in 2008 by Joosse
    plf <- read.delim(paste0(scriptdir,"/platformnki_hg18.txt"))

    # merge the input data with the platform file. The platform file contains duplicated chrom and maploc positions
    # which correspond to triple spotted array positions which are marked by a different order value. In the merge we use
    # all.x to retain the dimensions of the matrix.
    comb<- merge(plf, dt, by.x=c('chrom','maploc'), by.y=c('chrom', 'maploc'), all.x=T)
    comb <- comb[order(comb$chrom, comb$maploc), ]
    # linear interpolation for ovarian cancer due to high correlation between probes.
    comb[,-1:-8] <- na.approx(comb[,-1:-8], na.rm=F)
    comb[,-1:-8] <- na.locf(comb[,-1:-8], fromLast=T, na.rm=F)
    comb[,-1:-8] <- na.locf(comb[,-1:-8])
    return(comb[, colnames(comb) %in% colnames(dt)])
}

# correctDataset  ######################################################################################################
########################################################################################################################
##    Philip Schouten 2020 <philip.schouten@gmail.com>
# adjusted from Rubayte to include  breast and ovarian cancer correction
# apply centering and scaling correction of breast cancer and ovarian cancer data.
# Aim: correct centering and scaling differences that are assumed to 
# arise between different copy number profiling methods/preprocessing methods.
# the correction in this file is towards the BAC array data, and therefor applies 
# to the breast cancer classifiers as provided in the nkiBRCA package.

correctDataset <- function(dt, sample_type=c('breast','ovarian'),filedir) {
  #Author: Philip Schouten

  
  # load unsegmented breast or ovarian cancer ratios, depending on sample_type
  # breast  cancer unsegmented data:
  if (sample_type=='breast') {
      load(paste0(scriptdir, '/AnnOncB1PaperRatios.RDa'))
      refdata <- AnnOncB1PaperRatios
  }
  # ovarian cancer unsegmented data:
  if (sample_type=='ovarian') {
    load(paste0(scriptdir, '/ov.ratios.RDa'))
    refdata <- ov.ratios
  }


  # load platform file
  # platform file. These are the locations of the BACs. New data needs to be mapped
  # to chr1:22 and chrX described in here (3248 probes). N.B. this drops chrY however
  # Y should not be taken into account in the correction as in the sequencing data
  # has extreme values. Average of the ratios within
  # the BAC clone. E.g. if BAC clone is chr1 100000 - 200000, find all ratios within
  # this interval, take the mean, and use this as the ratio for this BAC location. 
  plf <- read.delim(paste0(scriptdir, '/platformnki_hg18.txt'), sep='\t', stringsAsFactors=F)
  
  # this is the data to correct
  newdata <- (dt[1:3248,-1:-2])
  meannew <- apply(newdata,1,mean) 
 
  # this is the reference data
  meanref <- apply(refdata[1:3248,-1:-2],1,mean)

  # for breast cancer correction based on ratios use:
  # - AnnOncB1PaperRatios
  # for breast cancer correction based on segments (not recommended) use:
  #- AnnOncB1PaperSegments
  # - for ovarian cancer correction based on ratios use: 
  # ov.ratios (ovr)
  # for ovarian cancer correction based on segments use:
  # ov.segments (ovs)


  # fit simple linear model for centering and scaling the new dataset to the 
  # platform used for training.

  # I'd prefer to fit the model on unsegmented data, correct and then segment the data
  # afterwards. If no unsegmented data is available it's possible to perform correction
  # with the segmented data. Always check afterwards that the average profiles better resemble
  # each other and for location specific differences
  ft <- glm(sort(meanref) ~ sort(meannew))

  # correct dataset
  correctednew <- (newdata  * coef(ft)[2]) +  coef(ft)[1]
  print(coef(ft))


  dt[1:3248, -1:-2] <- correctednew

  return(dt)
}

##############################################################################
## Classify 
##############################################################################

print(file)
print(sample_type)
print(cls)
print(correct_platform)
print(missing2centroid)
stopifnot(file.exists(file))
tmp <- read.delim(file, stringsAsFactors=F)
colnames(tmp)[3] <- 'chrom'

# the breast cancer classifier was built on hg 18 in 2008 (Joosse)
plf <- read.delim(paste0(scriptdir,"/platformnki_hg18.txt"))

comb <- merge(plf,tmp, by='Order', all.x=T)
comb <- comb[order(comb$chrom.x, comb$maploc), ]

kc <- data.frame(chrom=plf$chrom, maploc=plf$maploc, comb[,-1:-17])

print('created KC dataframe')

# platform correction, segmentation and classification can't handle NAs. Due to correlation between neighbouring locations
# linear interpolation is reasonable, alternative is to set the missing probe to the  average of the class centroids 
# (no effect of the probe for classifcation). For breast cancer this is handles in fixMissing2Centroid. For ovarian cancer
# in principle is handled by linear interpolation. 
# we save a matrix of missing values so breast cancer can be corrected using existing functions and for qc for ovarian cancer
  print('start pipeline')
  missing_mat <- is.na(as.matrix(kc))

  kc <- fillMat(kc)

  if (correct_platform) {
    kc <- correctDataset(kc, sample_type, filedir=scriptdir)
  }
  print('corrected dataset')

  # correct_platform needs to occur before changing the centroids of missing class otherwise the values will be 
  # adjusted by platform correction
  if (missing2centroid) {
    # reset the initial missing values
    kc <- as.matrix(kc)
    # for the BRCA2 breast cancer classifier kc will be overwritten and thus no ratios file will exist. 
    # therefore we save before resetting the missing values. In the ratios file the missings are interpolated and not set to centroid.
    # This is not fixable, as the centroid is set after segmentation. N.B. default is to let interpolation and segmentation 
    # and not to set to centroid.
    if (cls =='b2' && sample_type=='breast') {
      ratb2 <- kc
    }
    kc[missing_mat] <- NA
    kc <- data.frame(kc, stringsAsFactors=F)
    kc <- fixMissing2Centroid(cls=cls, dt=kc, fillM='ct', sample_type=sample_type)
    print('fixed missing to centroid') 
  }

  # not optional ; output requires sg
  #  if (segment) {
  # since missings are filled above, the only part of the script that is used is segmentation
  if (cls =='b2' && sample_type=='breast' && missing2centroid) {
    sg <- kc
  } else {
    sg <- fixMissing2Centroid(sample_type='ovarian', dt=kc)
  }
  print('segmented')

  # classifiers 
  if (cls == 'b1.191' && sample_type=='breast') {
    pred <- apply(kc[,-1:-2], 2, b1191)  
  } else if (cls == 'b1.371' && sample_type=='breast') {
    pred <- apply(kc[,-1:-2], 2, b1371)  
  } else if (cls == 'b2' && sample_type=='breast') {
    pred <- apply(kc[,-1:-2], 2, b2704)  
  } else if (cls == 'b1' &&  sample_type=='ovarian') {
    load(paste0(scriptdir,'/ov.B1.RDa'))
    pred <- with(ov.B1, pamr.predict(m, newx=as.matrix(sg[1:3248,-1:-2]), 
      threshold=delta[sel],type='posterior'))[,2]
  } else if (cls == 'b2' &&  sample_type=='ovarian') {
    load(paste0(scriptdir,'/ov.B2.RDa'))
    pred <- with(ov.B2, pamr.predict(m, newx=as.matrix(sg[1:3248,-1:-2]), 
      threshold=delta[sel],type='posterior'))[,2]
  } else {
    stop(paste(cls, "is not implemented for", sample_type))
  }
  print('created predictions')

if (exists('ratb2')) {
  ratios <- as.data.frame(ratb2[,-1:-2]);
} else {
  ratios <- as.data.frame(kc[,-1:-2]);
}
segments <- as.data.frame(sg[,-1:-2]);

write.table(ratios, file = paste0(outbase, ".ratios-", outbelong, ".tsv"), quote = FALSE, row.names = FALSE, sep = "\t")
write.table(segments, file = paste0(outbase, ".segments-", outbelong ,".tsv"), quote = FALSE, row.names = FALSE, sep = "\t")

pred.out <- data.frame(sample_id = colnames(ratios), class_probability=round(pred,3))
write.table(pred.out, file = paste0(outbase, ".pred-", outbelong ,".tsv"), quote = FALSE, row.names = FALSE, sep = "\t")

print('wrote output files')

print(sessionInfo())
print(paste0('opts:', opts))
sink()
}
