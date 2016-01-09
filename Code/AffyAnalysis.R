#Title: AffyAnalysis.R
#Author: SHannon Bell, ILS-INC, sbell@ils-inc.com
#
#
#performs RMA normalization for all arrays in a folder then
#Calcualtes the average probe intensity for each treatment so give both
#Individual animal probe intensity and the treatment average
#Parameters:
#d= the working directory
#od = the output directory for files that are to be collected
#exptinfo= a tab delim file of experimental information
#outdir = the output directory
#chem = the chemical name of interest
#organ = organ name of interest
#runID = run id value


################################
#get the parameters from the environment:
args=commandArgs(TRUE)
print(paste(commandArgs(), collapse = " "))
if(length(args)==0){
  print("No arguments supplied.")
  break
  ##supply default values
  d <-"Z:/TGGATES2/Data/Microarray/acetaminophen.Rat.in_vivo.Liver.Single/celfiles"
  exptinfo<-"Z:/TGGATES2/Files/Open-tggates_AllAttribute_EDT_RatLiver.txt"
  outdir <- 'Z:/TGGATES2/Output/test'
  chem<-'APAP'
  organ<-'Liver'
  runID<-'test'
  Repeat<-"Single"
  
}else{
  print(args)
  args2<-strsplit(args, " ")
  for(i in 1:length(args2)){
    eval(parse(text=args2[[i]]))
  }
}
#Import libraries
library(affy)
#read in experimental paramters for all the microarray experiments
allinfo<-read.table(exptinfo, sep='\t', header=TRUE, na.strings='', quote="\"", stringsAsFactors=FALSE)
#
#Move to the folder containing the data of interest
setwd(d)
print(getwd())

#analysis for affy data with rma normalization
data<-ReadAffy()
eset <- rma(data,verbose=FALSE, background=FALSE )#rma with no bkgd correction
#NOTE: eset from rma returns obj with expression values in log2 space
#create a file name and save the eset object
#fname<-paste("eset", organ, chem, paste(runID, 'txt', sep='.'), sep='-')
#write.exprs(eset, file=fname) # Writes expression values to working directory.
###############################################
#get the average for each treatment level including controls
#convert to dataframe
edata <- as.data.frame(exprs(eset))
#note the columns are the cel file names
#use this to extract array information from the information table
#since there may be samples from the same treatment group that dont have microarrays we will need to bring that in seperatly anyway
#so removing that material now
info<-subset(allinfo, BARCODE %in% as.character(gsub('.CEL','',colnames(edata)) ))[,c('folder','filename','BARCODE','ORGAN_ID','COMPOUND_NAME','COMPOUND.Abbr.','SPECIES','TEST_TYPE','SIN_REP_TYPE','SACRI_PERIOD','DOSE_LEVEL')]
info$Array_ID<-paste(info$BARCODE, '.CEL', sep='')
#Creating a treatment variable that idenitifies individuals getting the same treatment
#remove special characters from treatment as it will become column name
info$treatment<-gsub("[[:space:]]","",as.factor(paste(info$ORGAN_ID, chem, info$DOSE_LEVEL, info$SACRI_PERIOD, info$SIN_REP_TYPE, sep="_")))
#adding in experiment information, this will be several treatments by dose (still seperated by singel or repeat dosage)
# experiment<-gsub("[[:space:],%-/]","",as.factor(paste(info$ORGAN_ID, chem, info$SACRIFICE_PERIOD, info$SINGLE_REPEAT_TYPE, sep="_")))
info$experiment <-gsub("[[:space:]]","",as.factor(paste(info$ORGAN_ID, chem, info$SACRI_PERIOD, info$SIN_REP_TYPE, sep="_")))
info$RunID <- runID
RS<-unique(info$SIN_REP_TYPE)
info$ProcessFileInd<-paste("Ind",organ,chem,RS, paste(runID, 'txt', sep='.'), sep='-')

#write out the efile
write.table(edata, file=file.path(outdir,"ProcArray", unique(info$ProcessFileInd)), sep='\t', row.names=TRUE, col.names=TRUE, quote=FALSE) #mean expression values

###########################################
#This section gets the means for each treatment
treats<-unique(info$treatment)
ArrayAve<-NULL
for(i in 1:length(treats)){
  grp_arrays<-subset(info, treatment == treats[i])$Array_ID
  #updating to catch occurances where samples are missing
  if(length(grp_arrays) <2){
  	grp_ave<-rep(NA, nrow(edata))
  }
  else{
	grp_ave<-apply(edata[,colnames(edata) %in% grp_arrays],1,FUN=mean, na.rm=TRUE)
  }
  ArrayAve<-cbind(ArrayAve, grp_ave)
}
#addind in the identifiers, which correspond to the treatment parameters
colnames(ArrayAve)<-treats
#removal of columns lacking data
ArrayAve<-ArrayAve[,!apply(is.na(ArrayAve),2,all)]
#now write the file for use in integrating datasets
#info2<-subset(info, treatment %in% colnames(ArrayAve))#should be the same dimensions
info$ProcessFileAve<-paste("Mean",organ,chem,RS, paste(runID, 'txt', sep='.'), sep='-')
#write object to table to facilitate input and mergeing
#the pattern of filenames should facilitate the retreival and then the checking
write.table(ArrayAve, file=file.path(outdir, "ProcArray", unique(info$ProcessFileAve)), sep='\t', row.names=TRUE, col.names=TRUE, quote=FALSE) #mean expression values
#and to make a shorter info file that doesnt contain all the clinical chemistry values (which arent helpful here)
########################################

#colnames(info)
# [1] "folder"         "filename"       "BARCODE"        "ORGAN_ID"       "COMPOUND_NAME"  "COMPOUND.Abbr." "SPECIES"       
# [8] "TEST_TYPE"      "SIN_REP_TYPE"   "SACRI_PERIOD"   "DOSE_LEVEL"     "Array_ID"       "treatment"      "experiment"    
# [15] "RunID"          "ProcessFileInd" "ProcessFileAve"
#write the info to a specific file in a subfolder, this will be collected with all the data being processed
fname<-paste("RunInfo", organ,chem,Repeat, paste(runID, 'txt', sep='.'), sep='-')
write.table(info, file=file.path(outdir,"Info", fname), sep='\t', append=TRUE, row.names=FALSE, col.names=FALSE, quote=FALSE)

