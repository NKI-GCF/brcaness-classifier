#!/usr/bin/Rscript
initial.options <- commandArgs(trailingOnly = FALSE)
file.arg.name <- "--file="
script.name <- sub(file.arg.name, "", initial.options[grep(file.arg.name, initial.options)])
SCRIPTDIR <- dirname(script.name)

# from p.schouten/run311_312/
source(paste0(SCRIPTDIR, "/SeqCNV/SeqCNV.R"), chdir=TRUE)
options(stringsAsFactors=TRUE) # needed for loadCovData

opts <- commandArgs(trailingOnly = TRUE);
print(paste0("options: '", paste(opts, collapse="', '"), "'"))
usage <- "Rscript create_NKI_1m.R <length> [mm10/mm9/hg19(=default) [--noplot]]";
doplot <- TRUE

if (length(opts) < 1) {
  print(usage);
} else {
  bw <- opts[1];
  if (length(opts) > 1) {
    ref <- opts[2]
    if ((length(opts) > 2) & (opts[3] == "--noplot")) {
      doplot <- FALSE
    }
  } else {
    ref <- "hg19"
  }
  mapa51 <- paste0(SCRIPTDIR, "/mapa/mapa51-", bw,"K-", ref, ".bed.gz")
  print(paste("Using mapability file", mapa51))
  stopifnot(file.exists(mapa51))

  blacklist <- paste0(SCRIPTDIR, "/mapa/", ref, "-blacklist-nochr.bed")
  print(paste("Using mapability file", blacklist))
  stopifnot(file.exists(blacklist))

  mapa <- read.delim(gzfile(mapa51), header=FALSE, col.names=c("chr","start","end","mapa"))
  black <- read.delim(blacklist, header=F, col.names=c("chr", "start", "end"))

  f <- findCovFiles(pattern=paste0(".*counts-",bw, "000-q37.txt"))
  stopifnot(length(f) > 0)
  gc <- read.delim(paste0("/tmp/gccontent-",bw,"000.txt"))[,c(1,2,3,5)]
  colnames(gc) <- c("chr","start","end","gc")

  data <- loadCovData(f, gc=gc, mapa=mapa, black=black, exclude="MT", datacol=4)
  print("finished loading data")

  colnames(data$cov) <- gsub(paste0("_L00[1-8]-s-counts-", bw, "000-q37.txt"),"",f)
  #fix the names
  sampnames <- paste(sub("(_[CGAT]{6,}).*","\\1", sub("^.*/","", colnames(data$cov)), perl=T), paste0(round(colSums(data$cov) / 1e6,1),"M"), sep="_")
  colnames(data$cov) <- sampnames

  usepoints <-  !(data$anno$black | data$anno$chr %in% c("X", "Y", "MT"))

  ratios <- sapply(1:ncol(data$cov), function(s) {
    if(doplot == TRUE) {
      plots <- paste0("qc", bw, "K/", colnames(data$cov)[s], ".png")
    } else {
      plots <- NULL
    }
    tng(data.frame(count=data$cov[,s], gc=data$anno$gc, mapa=data$anno$mapa), use=usepoints, plot=plots)
  })
  print("finished creating ratios")

  colnames(ratios) <- sampnames
  rd <- list(ratios=ratios[!data$anno$black & !is.na(data$anno$mapa) & data$anno$mapa > .2, ], anno=data$anno[!data$anno$black & !is.na(data$anno$mapa) & data$anno$mapa > .2, ])

  #load bac nki platform
  read.delim(paste0(SCRIPTDIR, "/bac/platformnki_lift", ref, "-nochr.txt"), header=T)->bac
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
  as.data.frame(RangedData(bacr1m)),
  tobac[match(as.character(bac$Order), rownames(tobac)),])
  dir <- getwd()
  write.table(tab, file=paste0("NKI_1M_", gsub("^.*/([^/]+)$", "\\1", getwd()), ".txt"), sep="\t", quote=F, row.names=F)
}

