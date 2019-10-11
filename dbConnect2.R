library(RODBC)

dbConnect <- function (instance,UserUID,UserPWD){
	
	connStr <-paste0(
    		"Driver={ODBC Driver 17 for SQL Server};Server=",
    		instance,
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
