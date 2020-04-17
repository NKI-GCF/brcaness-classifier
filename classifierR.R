library(readxl)
library(optparse)

option_list = list(
  make_option(c("-i", "--in"), type="character", default=NULL,
              help="input dataset file name", metavar="character"),
  make_option(c("-s", "--snr"), type="character", default=NULL,
              help="sample number", metavar="character"),
  make_option(c("-t", "--type"), type="character", default=NULL,
              help="what type of data (mamma=1/ovarium=2)", metavar="character"),
  make_option(c("-b", "--brca"), type="character", default=NULL,
              help="what BRCA to test for (1 / 2 )", metavar="character"),
  make_option(c("-o", "--out"), type="character", default=NULL,
              help="output file name", metavar="character")

);
opt <- parse_args(OptionParser(option_list=option_list))

brcanum<-as.integer(opt$b)

if(opt$t==1){

  library(nkiBRCA)

  if(brcanum==1){
    newData<-as.data.frame(read_excel(opt$i,col_names = F,skip=1))
    print(b1191(newData[,as.integer(opt$s)]))

  }else if (brcanum==2){
    newData<-as.data.frame(read_excel(opt$i,col_names = F,skip=1, sheet="segments"))
    print(b2704(newData[,as.integer(opt$s)]))
  }

}else if (opt$t==2){
  library(pamr)
  if(brcanum==1){

    load('/app/ov.B1.RDa')
    #print(opt$i)
    newData<-as.data.frame(read_excel(opt$i,col_names = F,skip=1, sheet="segments" ))
    out<-with(ov.B1,pamr.predict(m,as.matrix(newData[,as.integer(opt$s)][c(1:3248)]),delta[sel],"posterior"))

  }else if (brcanum==2){

    load('/app/ov.B2.RDa')
    #print(opt$i)
    newData<-as.data.frame(read_excel(opt$i,col_names = F,skip=1, sheet="segments"))
    out<-with(ov.B2,pamr.predict(m,as.matrix(newData[,as.integer(opt$s)][c(1:3248)]),delta[sel],"posterior"))
  }
  print(out[2]) #this is what we want to have!
}

#print(out)







#notes:
#local({
#r <- getOption("repos");
#r["CRAN"] <- "https://cran.rstudio.com/"
#options(repos=r)
#})
