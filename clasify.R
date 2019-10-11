#'Classify BRCA1 and BRCA2 in ovarian and mama data.
#'
#'takes following input to generate the classification.
#'@param i location of the input data file
#'@param s what column of the input data file to process
#'@param t what type of data, 1 for mama data, 2 for ovarian data
#'@param b what type of BRCA to check 1 for BRCA1, 2 for BRCA2
#'@return returns the BRCA number, <0.5 probably no BRCA, >0.5 probably BRCA
#'@export
#'

classify <- function(i,s,t,b){




  library(readxl)
  brcanum<-as.integer(b)

  if(t==1){

    library(nkiBRCA)

    if(brcanum==1){
      newData<-as.data.frame(read_excel(i,col_names = F,skip=1))
      return(b1191(newData[,as.integer(s)]))

    }else if (brcanum==2){
      newData<-as.data.frame(read_excel(i,col_names = F,skip=1))
      return(b2704(newData[,as.integer(s)]))
    }

  }else if (t==2){
    library(pamr)
    if(brcanum==1){

      load('~/Classifier/ov.B1.RDa')
      #print(opt$i)
      newData<-as.data.frame(read_excel(i,col_names = F,skip=1))
      out<-with(ov.B1,pamr.predict(m,as.matrix(newData[,as.integer(s)][c(1:3248)]),delta[sel],"posterior"))
      return(out[2]) #this is what we want to have!

    }else if (brcanum==2){

      load('~/Classifier/ov.B2.RDa')
      #print(opt$i)
      newData<-as.data.frame(read_excel(i,col_names = F,skip=1))
      out<-with(ov.B2,pamr.predict(m,as.matrix(newData[,as.integer(s)][c(1:3248)]),delta[sel],"posterior"))
      return(out[2]) #this is what we want to have!
    }

  }
}

#print(out)







#notes:
#local({
#r <- getOption("repos");
#r["CRAN"] <- "https://cran.rstudio.com/"
#options(repos=r)
#})
