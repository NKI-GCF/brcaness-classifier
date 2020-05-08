#!/usr/bin/Rscript
sink("create_NKI_1m_sessioninfo.txt")
sink("create_NKI_1m_sessioninfo.txt", type="message")

  source(paste0("/app/SeqCNV/SeqCNV.R"), chdir=TRUE)

  mappa65 <- paste0("/app/mappability/GRCh38_65bp-q15-20k.bed.gz")
  print(paste("Using mappabilitybility file", mappa65))
  stopifnot(file.exists(mappa65))

  blacklist <- paste0("/app/mappability/GRCh38-blacklist-nochr.bed")
  print(paste("Using mappabilitybility file", blacklist))
  stopifnot(file.exists(blacklist))

  mappability <- read.delim(gzfile(mappa65), header=FALSE, col.names=c("chr","start","end","mappos", "mappa"))
  black <- read.delim(blacklist, header=F, col.names=c("chr", "start", "end"))

  f <- findCovFiles(pattern=paste0(".*counts-20000-q15.txt"))
  stopifnot(length(f) > 0)
  gc <- read.delim(paste0("/tmp/gccontent-20000.txt"))[,c(1,2,3,5)]
  colnames(gc) <- c("chr","start","end","gc")

  data <- loadCovData(f, gc=gc, mappa=mappability, black=black, exclude="MT", datacol=4)
  print("finished loading data")

  colnames(data$cov) <- gsub(paste0("_L00[1-8]-s-counts-20000-q15.txt"),"",f)
  #fix the names
  sampnames <- paste(sub("(_[CGAT]{6,}).*","\\1", sub("^.*/","", colnames(data$cov)), perl=T), paste0(round(colSums(data$cov) / 1e6,1),"M"), sep="_")
  colnames(data$cov) <- sampnames

  usepoints <-  !(data$anno$black | data$anno$chr %in% c("X", "Y", "MT"))

  ratios <- sapply(1:ncol(data$cov), function(s) {
    plots <- paste0("qc20K/", colnames(data$cov)[s], ".png")
    tng(data.frame(count=data$cov[,s], gc=data$anno$gc, mappa=data$anno$mappa), use=usepoints, plot=plots)
  })
  print("finished creating ratios")

  colnames(ratios) <- sampnames
  rd <- list(ratios=ratios[!data$anno$black & !is.na(data$anno$mappa) & data$anno$mappa > .2, ], anno=data$anno[!data$anno$black & !is.na(data$anno$mappa) & data$anno$mappa > .2, ])

  #load bac nki platform
  read.delim(paste0("/app/bac/platformnki_lifthg38-nochr.txt"), header=T)->bac

  # fix factor order
  bac$om <- droplevels(bac$om)
  bac$om <- factor(droplevels(bac$om), levels=c(1:22, "X", "Y"))
  bac <- bac[bac$Start > 0,]
  bacr <- RangedData(IRanges(start=bac$Start, end=bac$End), space=bac$om, id=bac$Order, clone=bac$Clone)

  #make 1M
  bacr1m <- resize(ranges(bacr), width=1e6, fix="center")

  #summarize (mean) data in rd using regions in bac
  res <- RangedData(IRanges(start=rd$anno$start, end=rd$anno$end), space=rd$anno$chr, rd$ratios)
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
  write.table(tab, file=paste0("NKI_1M.txt"), sep="\t", quote=F, row.names=F)

sessioInfo()
sink(type="message")
sink()

