#!/usr/bin/Rscript
opts <- commandArgs(trailingOnly = TRUE);

stopifnot(length(opts) == 5)

bw <- opts[1]
mapq <- opts[2]
seqlen <- opts[3]
blacklist = opts[4]
outbase <- opts[5]

if (outbase != "") {
  log_file <- file(paste0("/output/log/create_", outbase, "_sessioninfo.txt"), open="wt")
  sink(log_file)
  sink(log_file, type="message")
}
  source("/app/SeqCNV.R")
  require("IRanges")
  require("gtools")
  rd_file <- paste0("rd_", bw, ".Rda")
  if (file.exists(rd_file)) {
    load(rd_file)
  } else {

  mappa <- paste0("/app/GRCh38_", seqlen, "bp-q", mapq, "-", bw, "k.bed.gz")
  print(paste("Using mappabilitybility file", mappa))
  stopifnot(file.exists(mappa))

  blacklist <- paste0("/app/", blacklist)
  print(paste("Using mappabilitybility file", blacklist))
  stopifnot(file.exists(blacklist))

  mappability <- read.delim(gzfile(mappa), header=FALSE, col.names=c("chr","start","end","mappos", "mappa"), skip=1)
  black <- read.delim(blacklist, header=F, col.names=c("chr", "start", "end"))

  f <- findCovFiles(pattern=paste0(".*counts-", bw, "000-q", mapq, ".txt"))
  stopifnot(length(f) > 0)
  gc <- read.delim(paste0("/app/gccontent-", bw, "k.txt"))[,c(1,2,3,5)]
  colnames(gc) <- c("chr","start","end","gc")

  data <- loadCovData(f, gc=gc, mappa=mappability, black=black, exclude="MT", datacol=4)
  print("finished loading data")

  colnames(data$cov) <- gsub(paste0("-counts-", bw, "000-q", mapq, ".txt"),"",f)
  #fix the names
  sampnames <- paste(sub("(_[CGAT]{6,}).*","\\1", sub("^.*/","", colnames(data$cov)), perl=T), paste0(round(colSums(data$cov) / 1e6,1),"M"), sep="_")
  colnames(data$cov) <- sampnames

  usepoints <-  !(data$anno$black | data$anno$chr %in% c("X", "Y", "MT"))

  ratios <- sapply(1:ncol(data$cov), function(s) {
    plots <- if(outbase != "") {paste0("qc", bw, "K/", colnames(data$cov)[s], ".png")} else {NULL}
    tng(data.frame(count=data$cov[,s], gc=data$anno$gc, mappa=data$anno$mappa), use=usepoints, plots=plots)
  })
  print("finished creating ratios")

  colnames(ratios) <- sampnames
  rd <- list(ratios=ratios[!data$anno$black & !is.na(data$anno$mappa) & data$anno$mappa > .2, ], anno=data$anno[!data$anno$black & !is.na(data$anno$mappa) & data$anno$mappa > .2, ])

  # make ranges for later bac overlap

    save(rd, file=rd_file)
  }
  res <- RangedData(IRanges(start=rd$anno$start, end=rd$anno$end), space=rd$anno$chr, rd$ratios)

  rd$ratios_orig <- rd$ratios
  # correct linear fit for effects
  for (i in 1:ncol(rd$ratios)) {
    rd$ratios[,i] <- rd$ratios[,i] - median(rd$ratios[,i], na.rm=TRUE)
  }

  for (chr in c(1:22, "X")) {
    pdfname <- paste0("/output/chr", chr ,"-varmappa-", bw, "K.pdf")
    if(!file.exists(pdfname)) {
      print(pdfname)
      pdf(pdfname, width=0, height=0, paper="A4r")
      par(mfrow=c(1,1))
      for(i in seq.int(1,ncol(rd$ratios),4)) {
        for(p in i:min((i+3),ncol(rd$ratios))) {
          plotCNV(rd, sample=p, pch=".", main=paste("chromosome", chr, colnames(rd$ratios)[p]), ylim=c(-4,4), chr=chr)
        }
        if(dev.capabilities()$locator) readline()
      }
      dev.off()
    }
  }

  # create Nexus data
  dn <- paste0("/output/singlechannel_mappa.2_ratios", bw, "000")
  dir.create(dn, showWarnings=F)
  fns <- paste(colnames(rd$ratios), "txt", sep=".")
  writeAllNexusNormalized(rd, filenames = fns, path=dn)

  pdfname <- paste0("/output/all-varmappa-", bw, "K.pdf")
  if(!file.exists(pdfname)) {

    print(pdfname)
    pdf(pdfname, width=0, height=0, paper="A4r")
    par(mfrow=c(1,1))
    for(i in seq.int(1,ncol(rd$ratios),4)) {
      for(p in i:min((i+3),ncol(rd$ratios))) {
        plotCNV(rd[rd$anno$chr %in% c(1:22, "X")], sample=p, pch=".", main=colnames(rd$ratios)[p], ylim=c(-4,4))
      }
      if(dev.capabilities()$locator) readline()
    }
    dev.off()
  }

  stopifnot(outbase != "");

  #load bac nki platform
  read.delim("/app/platformnki_lifthg38-nochr.txt", header=T)->bac

  # fix factor order
  bac$om <- droplevels(bac$om)
  bac$om <- factor(droplevels(bac$om), levels=c(1:22, "X", "Y"))
  bac <- bac[bac$Start > 0,]
  bacr <- RangedData(IRanges(start=bac$Start, end=bac$End), space=bac$om, id=bac$Order, clone=bac$Clone)

  #make 1M
  bacr1m <- resize(ranges(bacr), width=1e6, fix="center")

  #summarize (mean) data in rd using regions in bac
  sn <- colnames(res)
  res$bacid <- rep(NA_integer_, nrow(res))
  overlaps <- as.matrix(findOverlaps(ranges(res), bacr1m, minoverlap=10000L, select="all"))
  res$bacid[overlaps[,1]] <- bacr$id[overlaps[,2]]

  resdf <- as.data.frame(res)

  tobac <- sapply(sn,function(n) tapply(resdf[,n], resdf$bacid, function(x) mean(x[is.finite(x)])))

  #final table
  tab <- data.frame(
  bac,
  as.data.frame(bacr1m),
  tobac[match(as.character(bac$Order), rownames(tobac)),])
  dir <- getwd()
  write.table(tab, file=paste0(outbase, ".txt"), sep="\t", quote=F, row.names=F)

sessionInfo()
sink(type="message")
sink()

