#setup to  reorg files for analysis

#
#for desktop-transfer file to destination or build there

baseDir<-"/wikihomedir/sbell/TG-Gates"
data<-read.table(file.path(baseDir,"Data/open_tggates_cel_file_attributeEDT.txt"), sep='\t', header=TRUE, na.strings='', quote="\"")

colnames(data)
# [1] "BARCODE"                   "ARR_DESIGN"               
# [3] "EXP_ID"                    "GROUP_ID"                 
# [5] "INDIVIDUAL_ID"             "ORGAN"                    
# [7] "MATERIAL_ID"               "COMPOUND_NAME"            
# [9] "COMPOUND_ABBREVIATION"     "COMPOUND_NO"              
# [11] "SPECIES"                   "EXP_TEST_TYPE"            
# [13] "SINGLE_REPEAT_TYPE"        "SEX_TYPE"                 
# [15] "STRAIN_TYPE"               "ADMINISTRATION_ROUTE_TYPE"
# [17] "ANIMAL_AGE.week."          "SACRIFICE_PERIOD"         
# [19] "DOSE_LEVEL" 
###
# #note that all the prior work usese "ORGAN_ID" and "TEST_TYPE" 
# #so those 2 names need to be changed so i can ensure that other code doesnt break
# y<-colnames(data)
# y<-sub("ORGAN","ORGAN_ID", y)
# y<-sub("EXP_TEST_TYPE","TEST_TYPE", y)
# colnames(data)<-y
# colnames(data)
# # [1] "BARCODE"                   "ARR_DESIGN"               
# # [3] "EXP_ID"                    "GROUP_ID"                 
# # [5] "INDIVIDUAL_ID"             "ORGAN_ID"                 
# # [7] "MATERIAL_ID"               "COMPOUND_NAME"            
# # [9] "COMPOUND_ABBREVIATION"     "COMPOUND_NO"              
# #[11] "SPECIES"                   "TEST_TYPE"                
# #[13] "SINGLE_REPEAT_TYPE"        "SEX_TYPE"                 
# #[15] "STRAIN_TYPE"               "ADMINISTRATION_ROUTE_TYPE"
# #[17] "ANIMAL_AGE.week."          "SACRIFICE_PERIOD"         
# #[19] "DOSE_LEVEL" 

library(doBy)
#note new file has different colnames so alway dbl check

test<-splitBy(BARCODE ~ SPECIES + EXP_TEST_TYPE + ORGAN + SINGLE_REPEAT_TYPE, data=data)
names(test)
# [1] "Rat|in vivo|Liver|Single"  "Rat|in vivo|Kidney|Single"
# [3] "Rat|in vivo|Liver|Repeat"  "Rat|in vivo|Kidney|Repeat"
# [5] "Rat|in vitro|Liver|NA"     "Human|in vitro|Liver|NA"  


#BREAKING EACH PART INTO THE ORGAN (FIRST FOLDER LEVEL)
#folder names/structure is SPECIES/TEST_TYPE/<ORGAN_ID>/<SINGLE_REPEAT_TYPE/
#fname is cmpd.species.test_type.organid
#or
#cmpd.species.test_type.organid.single.repeat type
#note that these arent needed, 
human<-test$"Human|in vitro|Liver|NA" 
ratCel<-test$"Rat|in vitro|Liver|NA"
kidneyR<-test$"Rat|in vivo|Kidney|Repeat"
kidneyS<-test$"Rat|in vivo|Kidney|Single"
liverR<-test$"Rat|in vivo|Liver|Repeat"
liverS<-test$"Rat|in vivo|Liver|Single"


#ok now for work
data2<-data
#remove odd characters such that chem name can match file name
data2$COMPOUND_NAME<-gsub("[[:space:]:]", '_', data$COMPOUND_NAME)
#remove " ' " and ()for processing...i dont think there are any but to be safe
data2$COMPOUND_NAME<-gsub("['()]", '', data2$COMPOUND_NAME)
#remove odd characters such that organ name can match file name
data2$TEST_TYPE<-gsub("[[:space:]]", '_', data$TEST_TYPE)
#to add on the missing leading 0 
data2$BARCODE <- paste("00", data$BARCODE, sep='')
#adding in the origional compound name for use later on
data2$CHEMICAL<-data$COMPOUND_NAME

#####################
# Folder structure and folder names are odd. Folder names are givn as :
# COMPOUND_NAME.SPECIES.TEST_TYPE.ORGAN_ID for invitro
#and COMPOUND_NAME.SPECIES.TEST_TYPE.ORGAN_ID.SINGLE_REPEAT_TYPE for invivo

spltData<-splitBy(BARCODE ~ SPECIES + TEST_TYPE + ORGAN_ID + COMPOUND_NAME + SINGLE_REPEAT_TYPE, data=data2)
ids<-NULL
filePath<-NULL
folder<-NULL
#fn<-NULL
for(i in 1:length(names(spltData))){
  temp<-spltData[[names(spltData)[i]]]
  spe<-unique(temp$SPECIES)
  tt<-unique(temp$TEST_TYPE)
  org<-unique(temp$ORGAN_ID)
  rep<-unique(temp$SINGLE_REPEAT_TYPE)
  cem<-unique(temp$COMPOUND_NAME)
  if(length(org)!=1 || length(cem)!=1){ 
     cat("error in:", org, " :: ", cem, '\n')
   }
  else if(tt == "in_vitro"){
    idstemp<-as.character(temp$BARCODE)
    ids<-c(ids, idstemp)
    fn<-paste(cem, spe, tt, org, sep=".")
    folder<-c(folder, rep(file.path(spe,tt,fn, "celfiles"), length(idstemp)))
    filePath<-c(filePath, file.path(spe, tt, fn, "celfiles", paste(idstemp, "CEL", sep=".")))
  }
  else if(tt == "in_vivo"){
    idstemp<-as.character(temp$BARCODE)
    ids<-c(ids, idstemp)
    fn<-paste(cem, spe, tt, org,rep, sep=".")
    folder<-c(folder, rep(file.path(spe,tt, org, rep,fn, "celfiles"), length(idstemp)))
    filePath<-c(filePath, file.path(spe, tt, org, rep, fn, "celfiles", paste(idstemp, "CEL", sep=".")))
  }
}

dataNew<-merge(data2, cbind(ids, folder, filePath), by.x='BARCODE', by.y='ids', all.y=TRUE)
#spot check
temp1<-subset(dataNew, COMPOUND_NAME == "phenylanthranilic_acid" & SINGLE_REPEAT_TYPE == "Repeat"& ORGAN_ID == "Liver")
dim(temp1)
#[1] 47 22
temp1[c(1,5,22,47),]


dim(dataNew)
#[1] 21385    22
colnames(dataNew)
# [1] "BARCODE"                   "ARR_DESIGN"               
# [3] "EXP_ID"                    "GROUP_ID"                 
# [5] "INDIVIDUAL_ID"             "ORGAN_ID"                 
# [7] "MATERIAL_ID"               "COMPOUND_NAME"            
# [9] "COMPOUND_ABBREVIATION"     "COMPOUND_NO"              
# [11] "SPECIES"                   "TEST_TYPE"                
# [13] "SINGLE_REPEAT_TYPE"        "SEX_TYPE"                 
# [15] "STRAIN_TYPE"               "ADMINISTRATION_ROUTE_TYPE"
# [17] "ANIMAL_AGE.week."          "SACRIFICE_PERIOD"         
# [19] "DOSE_LEVEL"                "CHEMICAL"                 
# [21] "folder"                    "filePath"                      
fn<-paste("Data/Microarray_Data_ALLLabels_", gsub('-','',Sys.Date()),".txt",sep="")
write.table(dataNew, file =file.path(baseDir, fn) , sep='\t', col.names=TRUE, row.names=FALSE, quote=FALSE)
Rat<-subset(dataNew, SPECIES == "Rat")
fn<-paste("Data/Microarray_Data_RatLabels_", gsub('-','',Sys.Date()),".txt",sep="")
write.table(Rat, file =file.path(baseDir, fn) , sep='\t', col.names=TRUE, row.names=FALSE, quote=FALSE)

############
#do to space issues i need to work with a ss of data.
#I will try just liver (since that is what I will work on first)
#and test out using sym links in my script
LiverData<-subset(Rat, ORGAN_ID == 'Liver' & SINGLE_REPEAT_TYPE == "Single")
dim(LiverData)
fn<-paste("Data/Microarray_Data_RatLabels_LIVER_Single", gsub('-','',Sys.Date()),".txt",sep="")
write.table(LiverData, file =file.path(baseDir, fn) , sep='\t', col.names=TRUE, row.names=FALSE, quote=FALSE)
#make a simple file for the inputs
inputs<-unique(Rat[, c("SPECIES","TEST_TYPE", "SINGLE_REPEAT_TYPE","ORGAN_ID", "COMPOUND_NAME","folder" )])
fn<-paste("Data/Microarray_Data_RatMAPPINGS_", gsub('-','',Sys.Date()),".txt",sep="")
write.table(inputs, file =file.path(baseDir, fn) , sep='\t', col.names=TRUE, row.names=FALSE, quote=FALSE)
#
LiverInputs<-subset(inputs, ORGAN_ID == 'Liver' & SINGLE_REPEAT_TYPE == "Single")
testInputs<-LiverInputs[c(1,10,20,40,50:55),]
fn<-paste("Data/Microarray_Data_RatMAPPINGS_LIVER_Single_", gsub('-','',Sys.Date()),".txt",sep="")
write.table(LiverInputs, file =file.path(baseDir, fn) , sep='\t', col.names=TRUE, row.names=FALSE, quote=FALSE)
#
fn<-paste("Data/Microarray_Data_RatMAPPINGS_TEST_", gsub('-','',Sys.Date()),".txt",sep="")
write.table(testInputs, file =file.path(baseDir, fn) , sep='\t', col.names=TRUE, row.names=FALSE, quote=FALSE)



