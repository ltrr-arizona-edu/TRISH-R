emssgUNH<- function(emssg,outDir){
# Write an error message to a specified output directory
# D. Meko; last revised 2021-12-27
#
# emmsg (xx?)s message, written to file "error.txt" 
# outDir (1x?)s output directory to which message file is written
#
# Returns emssg, which is just a string of the input message, regurgitated
# Write error file
fileErr<-file(paste(outDir,"error.txt",sep=""))
writeLines(c(emssg), fileErr)
close(fileErr)
return(emssg)
}