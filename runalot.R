install.packages("C:\\Users\\r.moritz\\Documents\\nkiClassifier_0.1.0.zip", repos=NULL)
install.packages("http://ccb.nki.nl/software/nkibrca/nkiBRCA_1.00.tar.gz",repos=NULL)
library(nkiClassifier)
install.packages("nkiBRCA")
library(nkiBRCA)
library(readxl)
library(pamr)
install.packages("pamr")

columnNames<-unlist(as.data.frame(read_excel("C:\\Users\\r.moritz\\Desktop\\Classifier\\AlleRatioData.xlsx", range="A1:ZZ1"))[1,])
BCRA1OV<-c()
BCRA2OV<-c()
BCRA1mama<-c()
BCRA2mama<-c()



for(i in 1:666){
  print(i)
  BCRA1OV[i]<-classify("C:\\Users\\r.moritz\\Desktop\\Classifier\\AlleRatioData.xlsx",i,2,1)
  BCRA2OV[i]<-classify("C:\\Users\\r.moritz\\Desktop\\Classifier\\AlleRatioData.xlsx",i,2,2)
  BCRA1mama[i]<-classify("C:\\Users\\r.moritz\\Desktop\\Classifier\\AlleRatioData.xlsx",i,1,1)
  BCRA2mama[i]<-classify("C:\\Users\\r.moritz\\Desktop\\Classifier\\AlleRatioData.xlsx",i,1,2)
}

totalDF<-as.data.frame(cbind(columnNames[c(1:666)],BCRA1OV[c(1:666)],BCRA2OV[c(1:666)],BCRA1mama[c(1:666)],BCRA2mama[c(1:666)]))
totalDF<-totalDF[,-c(1)]
colnames(totalDF)<-c("BCRA1OV","BCRA2OV","BCRA1mama","BCRA2mama")
write.csv(totalDF,file = "C:\\Users\\r.moritz\\Desktop\\Classifier\\uitslag.csv")
