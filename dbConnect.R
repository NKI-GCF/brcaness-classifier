library(RODBC)

dbConnect <- function (instance,database,UserUID,UserPWD){
	
	connStr <-paste0(
    		"Driver={ODBC Driver 17 for SQL Server};Server=",
    		instance,
    		";Database=",
    		database,
    		";Uid=",
    		UserUID,
    		";Pwd=",
    		UserPWD,
		";"
  	)
	print (connStr)
	con <- odbcDriverConnect(connStr)
	return(con)

}

insertQuery <- function (query,con){
	sqlQuery(con,query)
}

dbConnectionClose <- function(con){
	odbcClose(con)
}
