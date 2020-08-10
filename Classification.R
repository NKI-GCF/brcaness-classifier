#!/export/data/apps/R-3.1.3/bin/Rscript

# ask Roel whether parallel is already incorporated, turned of for processing with windows
#require(multicore)

# Range methods
require(GenomicRanges)
# interpolation and padding
require(zoo)
# still required?
require(WriteXLS)
# shrunken centroid classification
require(pamr)
# BRCA1 and BRCA2 breast cancer classifiers
require(nkiBRCA)

#get directory of script
opts <- commandArgs(trailingOnly = FALSE)
file.arg.name <- "--file="
script.name <- sub(file.arg.name, "", opts[grep(file.arg.name, opts)])
scriptdir <- dirname(script.name)
print(paste0("scriptdir: ",scriptdir))
pipeline <- dirname(scriptdir)
print(paste0('pipeline:',pipeline))

m <- match("--args", opts, 0L)
if (m) {
  opts <- opts[-seq_len(m)]
} else {
  opts <- character()
}

print(paste0("options: '", paste(opts, collapse="', '"), "'"))
usage <- "Rscript cmd_EL.R <NKI_1M file> <sample_type> <cls> <variation_pipeline> <correct_platform> <missing2centroid";



# chris_cghsegmentation_functions.R
# -------------------------------------------------------------------
# Copyright 2011 Christiaan Klijn <c.klijn@nki.nl>
# Project: CGH segmentation
# Description: Functions to segment your CGH data
# -------------------------------------------------------------------

doCghSeg <- function(chromNum, allKC, chrom, maploc) {

  require(cghseg)
  require(GenomicRanges)

  tumnames <- colnames(allKC)[-1:-2]
  profiles <- as.matrix(allKC[chrom == chromNum, -c(1,2)])
  subchrom <- chrom[chrom==chromNum]
  submaploc <- maploc[chrom==chromNum]

  print(paste('chrom', chromNum,'is ordered:', !any(!order(submaploc)==1:length(submaploc))))

  # Segmentation
  n <- nrow(profiles)
  CGHd <- new("CGHdata",Y=profiles)
  CGHo <- new("CGHoptions")
  alpha(CGHo) <- 0.1 ## look for 1 to a maximum of 500 segments per tumor,
  select(CGHo) <- "mBIC" ## this define the way the number of segment is selected ( in between 1 and Kmax)
  #wavenorm(CGHo)<-'position' ## Wave correction using position

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
    ID=allID)

  # Make a range object to contain the start and the end probe for each segment

  return(segments)
}


#######################################
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

  tmp1 <- lapply(unique(as.character(chrom)), function(x) rep(segMeans[[x]], times=((endInd[[x]]+1)-startInd[[x]])))

  if(!length(unlist(tmp1)) == length(chrom)*length(unique(unlist(patients)))) {
    # duplicated maplocs in nki cause this problem
    # only chrom 3 has this problem:
    # lapply(1:22, function(x) duplicated(maploc[chrom==x]))
    # change nki chrom 3 maploc 249
    stop('missing segMeans for some probes')
  }

  names(tmp1) <- unique(as.character(chrom))
  tmp2 <- lapply(unique(as.character(chrom)), function(x) matrix(tmp1[[x]], ncol= length(unique(patients[[x]]))))


  final <- do.call(rbind, tmp2)

  final2 <- data.frame(chrom=chrom, maploc=maploc, final, stringsAsFactors=F)
  colnames(final2) <- coln

  return(final2)
}


#######################################
##    Philip Schouten 2013 <philip.schouten@gmail.com>
##    find rows with all values missing and put them to classifier centroid
##    mean. Finds single NAs and replace with linear approximation through the
##    zoo package or to the centroid with ct.

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
            
        # De classifier (in 2008 door Simon Joosse ontwikkeld) is nog op 18 gebouwd.
    plf <- read.delim(paste0(pipeline,"/bac/platformnki_hg18.txt"))

        # merge the input data with the platform file. The platform file contains duplicated chrom and maploc positions
        # which correspond to triple spotted array positions which are marked by a different order value. In the merge we use
        # all.x to retain the dimensions of the matrix.
    comb<- merge(plf, dt, by.x=c('chrom','maploc'), by.y=c('chrom', 'maploc'), all.x=T)
    comb <- comb[order(comb$chrom, comb$maploc), ]
    comb[comb$chrom == 24,-1:-8] <- 0        
    # we might have non-mappable positions; all values in a row are NA or random missings in a sample; some values in a row are NA
    allmissing <- apply(comb[, -1:-8], 1, function(x) all(is.na(x)))
    somemissing <- apply(comb[,-1:-8], 1, function(x) any(is.na(x))) & ! allmissing
    pos_some_missing <- which(is.na(as.matrix(comb[somemissing,-1:-8])))
  
    # set NA's on the Y (=24) chrom, to 0, it lacks classifier probes anyway

    
    # note for later: it would make sense to implement platform correction here. Think of how chromosome Y needs to be handled.
    # if implemented here the missing classifier probes will be correctly changed to the average of the class centroids

  
    if (sample_type == 'breast')  {
      # for breast cancer samples probes that contain only missings (nonmappable) should be set to the
      # to mean of the sum of the class centroids to prevent influence on the classification (Schouten et al BCRT 2012)
      # here might also be implemented to interpolate the missing probes. 
          
      # get centroids of the various classifiers, get the number of columns and get the classifier probes
      if (cls=='b1.371') {
        ct <- read.delim(paste0(pipeline,"/bac/B1371centroids.txt"), stringsAsFactors=F, header=F)
        n <- 2
        classprobe <- ! (ct[1] == ct[2])
        }
    
      if (cls=='b1.191') {
        ct <- read.delim(paste0(pipeline,"/bac/B1centroids.txt"), stringsAsFactors=F, header=F)
        n <- 3
        classprobe <- ! (ct[2] == ct[3])
        }
    
      if (cls=='b2') {
        ct <- read.delim(paste0(pipeline,"/bac/B2centroids.txt"), stringsAsFactors=F, header=F)
        n <- 3
        classprobe <- ! (ct[2] == ct[3])

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

 	#create a filling object
 	fl <- matrix(rep((ct[,n-1]+ ct[,n])/2,   ncol(comb)-8), ncol=(ncol(comb)-8))
			
	# fill individual missing probes with classifier centroids in case ct is chosen
	if(fillM=='ct') {
        tm <- as.matrix(comb[,-1:-8])
        tm[allmissing,] <- fl[allmissing,]
        # change somemissing to interpolated, they're randomly missing and therefore should not be filled with centroid avg.
        # set the segmented data back to NA for BRCA2 (for BRCA1 does not do anything)
        tm[somemissing,][pos_some_missing] <- NA
        comb[,-1:-8] <- tm
        # fill random missings with interpolation
        comb[,-1:-8] <- na.approx(comb[,-1:-8], na.rm=F)
        comb[,-1:-8] <- na.locf(comb[,-1:-8], fromLast=T,na.rm=F)
        comb[,-1:-8] <- na.locf(comb[,-1:-8])
    } else if (fillM=='zoo') {					
    # fill with linear interpolation in case 'zoo' is chosen.	
        comb[,-1:-8] <- na.approx(comb[,-1:-8], na.rm=F)
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

# fill missing by linear interpolation and subsequently fill trailing and starting NAs by bringing backward and forward
# the last known 
fillMat <- function(dt) {
        # De classifier (in 2008 door Simon Joosse ontwikkeld) is nog op 18 gebouwd.
    plf <- read.delim(paste0(pipeline,"/bac/platformnki_hg18.txt"))

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

# adjusted from Rubayte to include  breast and ovarian cancer correction
correctDataset <- function(dt, sample_type=c('breast','ovarian'),filedir) {
  #Author: Philip Schouten
  # Aim: correct centering and scaling differences that are assumed to 
  # arise between different copy number profiling methods/preprocessing methods.
  # the correction in this file is towards the BAC array data, and therefor applies 
  # to the breast cancer classifiers as provided in the nkiBRCA package.
  
  # load unsegmented breast or ovarian cancer ratios, depending on sample_type
  # breast  cancer unsegmented data:
  if (sample_type=='breast') {
      load(paste0(filedir, '/AnnOncB1PaperRatios.RDa'))
      ovr <- AnnOncB1PaperRatios
  }
  # ovarian cancer unsegmented data:
  if (sample_type=='ovarian') {
    load(paste0(filedir, '/ov.ratios.RDa'))
    ovr <- ov.ratios
  }


  # load platform file
  # platform file. These are the locations of the BACs. New data needs to be mapped
  # to chr1:22 and chrX described in here (3248 probes). Average of the ratios within
  # the BAC clone. E.g. if BAC clone is chr1 100000 - 200000, find all ratios within
  # this interval, take the mean, and use this as the ratio for this BAC location. 
  plf <- read.delim(paste0(filedir, '/platformnki.txt'), sep='\t', stringsAsFactors=F)
 

 
  # create example ratios
  fakeratios <- (dt[1:3248,-1:-2])
  meanfake <- apply(fakeratios,1,mean) 

 
  # mean ratios and segmented ratios
  meanbr <- apply(ovr[1:3248,-1:-2],1,mean)

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
  ft <- glm(sort(meanbr) ~ sort(meanfake))

  # correct dataset
  correctedfake <- (fakeratios  * coef(ft)[2]) +  coef(ft)[1]
  print(coef(ft))


  dt[1:3248, -1:-2] <- correctedfake

  return(dt)
}


if (length(opts) < 4) {
  print(usage);
} else {
  file <- opts[1];
  sample_type <- opts[2];
  cls <- opts[3];
  variation_pipeline <- as.logical(opts[4]);
  correct_platform <- as.logical(opts[5]);
  missing2centroid <- as.logical(opts[6]);
  print(file)
  print(sample_type)
  print(cls)
  print(variation_pipeline)
  print(correct_platform)
  print(missing2centroid)
   stopifnot(file.exists(file))
  tmp <- read.delim(file, stringsAsFactors=F)
  print(head(colnames(tmp)))
  colnames(tmp)[3] <- 'chrom'
  # De classifier (in 2008 door Simon Joosse ontwikkelt) is nog op 18 gebouwd.
  plf <- read.delim(paste0(pipeline,"/bac/platformnki_hg18.txt"))

  comb <- merge(plf,tmp, by='Order', all.x=T)
  comb <- comb[order(comb$chrom.x, comb$maploc), ]

  print(head(colnames(comb)),20)

  kc <- data.frame(chrom=plf$chrom, maploc=plf$maploc, comb[,-1:-17])
  print(colnames(comb[1:18]))
  print(head(colnames(kc),4))

  # fixed pipelines (these could also be called with arguments in the process.sh and the correct ratio/segment output
  # of variations below (line 285)):
  if (! variation_pipeline) {
    if (sample_type=='breast' & cls =='b1.191') {
      # flow 1 breast BRCA1 191 without platform correction
      kc <- fixMissing2Centroid('b1.191', dt=kc, fillM='ct', sample_type='breast')
      print(any(is.na(as.matrix(kc))))
      # since kc contains no missings the fixMissing2Centroid(sample_type='ovarian' only does segmentation)
      sg <- fixMissing2Centroid(dt=kc,  sample_type='ovarian')
      pred <- apply(kc[,-1:-2], 2, b1191)
      }

    if (sample_type=='breast' & cls =='b1.371')  {
      # flow 2 breast BRCA1 371 without platform correction
      kc <- fixMissing2Centroid('b1.371', dt=kc,fillM='ct', sample_type='breast')
      print(any(is.na(as.matrix(kc))))
      # since kc contains no missings the fixMissing2Centroid(sample_type='ovarian' only does segmentation)
      sg <- fixMissing2Centroid(dt=kc,  sample_type='ovarian')
      pred <- apply(kc[,-1:-2], 2, b1371)
      }

    if (sample_type=='breast' & cls =='b2')  {
      # flow 3 breast BRCA2 without platform correction
      sg <- fixMissing2Centroid('b2', dt=kc,fillM='ct', sample_type='breast')
      print(any(is.na(as.matrix(kc))))
      # fill missings by interpolation (this needs to be done here because it overwrites kc)
      kc <- fillMat(kc)
      pred <- apply(sg[,-1:-2], 2, b2704)
      }  

    if (sample_type=='ovarian' & cls=='b1') {
      # flow 4 ovarian BRCA1 with platform correction
      # fill missings by linear interpolation
      kc <- fillMat(kc[1:3248,])
      # platformcorrection ovarian
      kc <- correctDataset(kc[1:3248,], sample_type='ovarian', filedir=paste0(pipeline,'/brcaClassifier'))
      # now only uses segmentation from fixMissing2Centroid
      print(any(is.na(kc)))
      print(head(kc[kc$chrom == 22,]))
      sg <- fixMissing2Centroid(cls=cls, sample_type='ovarian', dt=kc)
      load(paste0(pipeline,'/brcaClassifier/ov.B1.RDa'))
      pred <- with(ov.B1, pamr.predict(m, newx=as.matrix(sg[1:3248,-1:-2]), 
        threshold=delta[sel],type='posterior'))[,2]
      }

      if (sample_type=='ovarian' & cls=='b2') {
      # flow 5 ovarian BRCA2 with platform correction
      # fill missings by linear interpolation
      kc <- fillMat(kc[1:3248,])
      # platformcorrection ovarian
      kc <- correctDataset(kc[1:3248,], sample_type='ovarian', filedir=paste0(pipeline,'/brcaClassifier'))
      # now only uses segmentation from fixMissing2Centroid
      sg <- fixMissing2Centroid(cls=cls,sample_type='ovarian', dt=kc)
      load(paste0(pipeline,'/brcaClassifier/ov.B2.RDa'))      
      pred <- with(ov.B2, pamr.predict(m, newx=as.matrix(sg[1:3248,-1:-2]), 
        threshold=delta[sel], type='posterior'))[,2]
      }
  }
   
  # return variations of the ratios and segments for plotting/classifation outside of the fixed pipelines

  # platform correction, segmentation and classification can't handle NAs. Due to correlation between neighbouring locations
  # linear interpolation is reasonable, alternative is to set the missing probe to the  average of the class centroids 
  # (no effect of the probe for classifcation). For breast cancer this is handles in fixMissing2Centroid. For ovarian cancer
  # in principle is handled by linear interpolation. 
  # we save a matrix of missing values so breast cancer can be corrected using existing functions and for qc for ovarian cancer
  if (variation_pipeline) {  
    missing_mat <- is.na(as.matrix(kc))

    kc <- fillMat(kc)

    if (correct_platform) {
      kc <- correctDataset(kc, sample_type, filedir=paste0(pipeline,'/brcaClassifier'))
    }

    # correct_platform needs to occur before changing the centroids of missing class otherwise the values will be 
    # adjusted by platform correction
    if (missing2centroid) {
      # reset the initial missing values
      kc <- as.matrix(kc)
      kc[missing_mat] <- NA
      kc <- data.frame(kc, stringsAsFactors=F)
      kc <- fixMissing2Centroid(cls=cls, dt=kc, fillM='ct', sample_type=sample_type)
    }

    # not optional ; output requires sg
    #  if (segment) {
    # since missings are filled above, the only part of the script that is used is segmentation
    if (cls =='b2' && sample_type=='breast' && missing2centroid) {
      print(head(kc[,1:5]))
      sg <- kc
      print('no')
    } else {
      print(any(is.na(kc)))
      print(head(kc[kc$chrom == 22,]))
      sg <- fixMissing2Centroid(sample_type='ovarian', dt=kc)
      print('yes')
    }
  }

  # classifiers could be run here to also produce pred.

  if (exists("pred")) {
      print(head(pred))
      print(dim(pred))
  }

  # since we do classification in the pipeline external classification is not necessary anymore,
  # therefor statistics page is irrelevant
  ratios <- as.data.frame(kc[,-1:-2]);
  segments <- as.data.frame(sg[,-1:-2]);
  # cnames <- c(" ", "Hybridisation quality", "Profile quality", "Mean standard deviation",
  #   "Mean intensity channel 1", "Mean intensity channel 2")
  # str(cnames)
  # statistics <- as.data.frame(matrix(data = 1, nrow = ncol(segments), ncol = length(cnames)))
  # statistics[1] <- colnames(segments)
  # colnames(statistics) <- cnames

  # add to filename the settings of the pipeline that produced the output.

  print(dim(ratios))
  print(head(ratios[,1:5]))
  print(dim(segments))
  print(head(segments[,1:5]))
  print(colnames(kc)[1:6])
  
  print(paste0(gsub("\\.txt","", basename(file)), ".ratios.tsv"))

  print(paste0('table of all equal ratios:', table(apply(ratios[,-1:-2],1,function(x) all(x==x[1])))))
  print(paste0('table of all equal segments:', table(apply(segments[,-1:-2],1,function(x) all(x==x[1])))))


  #WriteXLS(c("ratios", "segments", "statistics"), paste0(gsub("\\.txt","", basename(file)), ".xls"))
  write.table(ratios, file = paste0(gsub("\\.txt","", basename(file)), ".ratios-", cls,'-',sample_type,'-',variation_pipeline,'-',correct_platform,'-',missing2centroid, ".tsv"), quote = FALSE, row.names = FALSE, sep = "\t")
  write.table(segments, file = paste0(gsub("\\.txt","", basename(file)), ".segments-", cls,'-',sample_type,'-',variation_pipeline,'-',correct_platform,'-',missing2centroid ,".tsv"), quote = FALSE, row.names = FALSE, sep = "\t")
  #write.table(statistics, file = paste0(gsub("\\.txt","", basename(file)), ".statistics.tsv"), quote = FALSE, row.names = FALSE, sep = "\t")

  if (exists("pred"))  {
    pred.out <- data.frame(sample_id = colnames(ratios), class_probability=round(pred,3))
    print(dim(pred.out))
    write.table(pred.out, file = paste0(gsub("\\.txt","", basename(file)), ".pred.", "-", cls,'-',sample_type,'-',variation_pipeline,'-',correct_platform,'-',missing2centroid ,".tsv"), quote = FALSE, row.names = FALSE, sep = "\t")
  }

  sink(file=paste0(gsub("\\.txt","", basename(file)), ".session.", sample_type,"-", cls,".txt"))
  print(sessionInfo())
  print(paste0('opts:', opts))
  sink()
}