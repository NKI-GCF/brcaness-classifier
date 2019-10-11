##########################################################################################################
##########################################################################################################
##########################################################################################################
###                                                                                                    ###
###                                                                                                    ###
###        Runs the classifier and gets the scores of all samples in the newest files.                 ###
###        After this it writes the score to the molpa DB                                              ###
###                                                                                                    ###
###                                                                                                    ###
##########################################################################################################
##########################################################################################################
##########################################################################################################

###Load Libraries and scripts needed###
library(readxl)
library(rjson)
library(optparse)

###Arguments to load###
option_list = list(
	make_option(c("-i", "--in"), type="character", default=NULL,help="input file", metavar="character"),
	make_option(c("-c", "--config"), type="character", default=NULL,help="input db config file", metavar="character"),
 	make_option(c("-l", "--inl"), type="character", default=NULL,help="input lib path", metavar="character"),
  	make_option(c("-d", "--db"), action="store_true", default=FALSE,help="insert output to database (default %default)")
);
opt <- parse_args(OptionParser(option_list=option_list))

###Read Arguments###
fileToRun <- opt$`in`
libPath <- opt$`inl`
dbc <- opt$`c`
procedure = ""

## source the libraries
if (dir.exists(libPath)){
  	source(paste0(libPath,"/",'classifierPackage.R'))
  
}

###make DB connections###
if(opt$d){
  	if (file.exists(dbc)){
    		dbConfig <- fromJSON(paste(readLines(dbc), collapse=""))

  	}else{
  	  	quit(status=1)
  	}
}


###old code might be handy for later use###
#p <- '//export//data//CNVSeq_Share//'
#paths <- dir(p, full.names=TRUE)
#fileInf <-file.info(paths)
#newestIndex <- difftime(Sys.time(), fileInf[,"mtime"], units = "days") < 20 #20 days difference!
#newestFiles<-paths[newestIndex]
#for(fileI in newestFiles){
###end old code###

###Start of the script by checking if the directories and files that are needed are in place###
#if(dir.exists(fileI)){
#pathFurther <- dir(paste(fileI,'//CV_01//',sep=''),full.names=TRUE)
#pathOfFile <-paste(pathFurther,'//Data//Intensities//BaseCalls//Alignment//',sep='')
#fileToRun <- paste(pathOfFile,'NKI_1M_CV_01.xls',sep='')

###Check if the file needed is in place###
#if(file.exists(fileToRun) && !file.exists(paste(pathFurther,'//Data//Intensities//BaseCalls//Alignment//Scored',sep='') )){
if (file.exists(fileToRun)){
  	###File found###
  	numberCols<-ncol(read_excel(fileToRun))
  	for(cols in 2:numberCols){

    		###found number of patients###
    		BRCA1<-classify(i=fileToRun,s=cols,t=1,b=1)
    		BRCA2<-classify(i=fileToRun,s=cols,t=1,b=2)

    		###puts the values into vectors and combines them to a Dataframe###
    		BRCA1Vect<-c(BRCA1$name,BRCA1$type,BRCA1$score)
    		BRCA2Vect<-c(BRCA2$name,BRCA2$type,BRCA2$score)
    		fileDataFrame<-data.frame(BRCA1Vect,BRCA2Vect)
    		names(fileDataFrame)<- c("BRCA1 values","BRCA2 values")

    		###the dataframe is writen to the given filename + output.csv###
    		write.table(fileDataFrame,file = paste(strsplit(fileToRun, ".", fixed = TRUE)[[1]][1],".output.csv",sep=""),append = TRUE)

    		#print(BRCA1)
    		#print(BRCA2)
    		TotalNameOfAnalysis<-BRCA1$name

   	 	###Made Classify score###
    		#print(substr(TotalNameOfAnalysis,1,1))
    		if(substr(TotalNameOfAnalysis,1,1)=='M'){

      			###correct Mnumber found###
      			TotalVect<-unlist(strsplit(TotalNameOfAnalysis,'.',T))
      			MNumber<-paste(TotalVect[1],"-",TotalVect[2],sep = '')
      			TestId<-unlist(strsplit(TotalVect[3],'_',T))[1]

      			###Make correct query###
      			#if(opt$d){
      			#  
      			#  dbCon<-paste0(libPath,'DataBaseConnect.py')
      			#  
			#	query<- paste0("UPDATE [MolPaApp].[dbo].[Test] set BRCA1 = ",BRCA1$score,", BRCA2 =",BRCA2$score," WHERE TestID=",TestId)
        #system(paste0("python ",dbCon," -s ",dbConfig$instance," -u ",dbConfig$UserUID," -d ",dbConfig$database," -p ",dbConfig$UserPWD," -q '",query,"'"))			
			#	#insertQuery(query,con)
      			#}
			###Query has been ran###

      			###unzipping needed files (needed for rashmie she can not unzip it###
      			#unzip(paste(pathOfFile,'singlechannel_mapa.2_ratios20000.zip',sep = ''),exdir=pathOfFile)

    		}

    		#uncomment if the scoring is completed.
    		#file.create(paste(pathFurther,"\\Data\\Intensities\\BaseCalls\\Alignment\\Scored",sep=''))

  	}
}
#}


###Complete###
