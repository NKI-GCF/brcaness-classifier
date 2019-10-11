#'Classify BRCA1 and BRCA2 in ovarian and mama data.
#'
#'The object is accessable with the following classify$name classify$type classify$score
#'
#'takes following input to generate the classification.
#'@param i location of the input data file
#'@param s what column of the input data file to process
#'@param t what type of data, 1 for mama data, 2 for ovarian data
#'@param b what type of BRCA to check 1 for BRCA1, 2 for BRCA2
#'@return returns the BRCA number, <0.5 probably no BRCA, >0.5 probably BRCA
#'@export
#'


###Get all arguments###
classify <- function(i,s,t,b){


toReturn<-list(score= NULL,name= "",type="")

###Get all libraries and dependencies###
library(readxl)
library(RODBC)

###Interpertate arguments###
brcanum<-as.integer(b)
SampleName<-colnames(read_excel(i))[s]


###If type is Mama BRCA classification###
if(t==1){

  ###only load the Mama BRCA library now save resources if ovarian would be done###
  library(nkiBRCA)

  ###Check if the BRCA1 or BRCA2 classification should be done###
  if(brcanum==1){
    ###Put the right data in the object###
    newData<-as.data.frame(read_excel(i))
    toReturn$name<-colnames(newData)[s]
    toReturn$type<-"mama"
    toReturn$score<-b1191(newData[1:nrow(newData),as.integer(s)])
    ###Return the object###
    return(toReturn) ###this is what we want to have!###

  }else if (brcanum==2){
    ###Put the right data in the object###
    newData<-as.data.frame(read_excel(i, sheet="segments"))
    toReturn$name<-colnames(newData)[s]
    toReturn$type<-"mama"
    toReturn$score<-b2704(newData[1:nrow(newData),as.integer(s)])
    ###return the object
    return(toReturn) ###this is what we want to have!###
  }

}else if (t==2){
  ###Ovarian BRCA1 classification###
  library(pamr)
  ###load PamR needed for the classification###
  if(brcanum==1){
    ###BRCA1 ovarian classification###
    load('/export/data/apps/cnvseq/CV_02/brcaClassifier/ov.B1.RDa')
    ###do the scoring and set the object right###
    newData<-as.data.frame(read_excel(i,sheet="segments" ))
    toReturn$name<-colnames(newData)[s]
    toReturn$type<-"ovarian"
    out<-with(ov.B1,pamr.predict(m,as.matrix(newData[1:nrow(newData),as.integer(s)][c(1:3248)]),delta[sel],"posterior"))
    toReturn$score<-out[2]
    ###The scoring is done###
    return(toReturn) ###this is what we want to have!###

  }else if (brcanum==2){
    ###Ovarian BRCA2 classification###
    load('/export/data/apps/cnvseq/CV_02/brcaClassifier/ov.B2.RDa')
    #print(opt$i)
    ###Do the scorring and set the object right###
    newData<-as.data.frame(read_excel(i, sheet="segments" ))
    toReturn$name<-colnames(newData)[s]
    toReturn$type<-"ovarian"
    out<-with(ov.B2,pamr.predict(m,as.matrix(newData[1:nrow(newData),as.integer(s)][c(1:3248)]),delta[sel],"posterior"))
    toReturn$score<-out[2]
    ###The scoring is done###
    return(toReturn) ###this is what we want to have!###
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
