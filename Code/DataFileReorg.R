#setup to  reorg files for analysis

#
#for desktop-transfer file to destination or build there
#libraries:
library(doBy)
#set the path
baseDir<-"//IT2/NICEATM_Data/TGGATES2"

#read in the data
#this is afile with all the meta data for the experiments along with the clinical chemistry and the weights
data<-read.table(file.path(baseDir,"/Data/Open-tggates_AllAttribute.tsv"), sep='\t', header=TRUE, na.strings='', quote="\"")

colnames(data)
# [1] "BARCODE"          "ARR_DESIGN"       "EXP_ID"           "GROUP_ID"         "INDIVIDUAL_ID"    "ORGAN_ID"        
# [7] "MATERIAL_ID"      "COMPOUND_NAME"    "COMPOUND.Abbr."   "COMPOUND_NO"      "SPECIES"          "TEST_TYPE"       
# [13] "SIN_REP_TYPE"     "SEX_TYPE"         "STRAIN_TYPE"      "ADM_ROUTE_TYPE"   "ANIMAL_AGE.week." "SACRI_PERIOD"    
# [19] "DOSE"             "DOSE_UNIT"        "DOSE_LEVEL"       "TERMINAL_BW.g."   "LIVER.g."         "KIDNEY_TOTAL.g." 
# [25] "KIDNEY_R.g."      "KIDNEY_L.g."      "RBC.x10_4.ul."    "Hb.g.dL."         "Ht..."            "MCV.fL."         
# [31] "MCH.pg."          "MCHC..."          "Ret..."           "Plat.x10_4.uL."   "WBC.x10_2.uL."    "Neu..."          
# [37] "Eos..."           "Bas..."           "Mono..."          "Lym..."           "PT.s."            "APTT.s."         
# [43] "Fbg.mg.dL."       "ALP.IU.L."        "TC.mg.dL."        "TG.mg.dL."        "PL.mg.dL."        "TBIL.mg.dL."     
# [49] "DBIL.mg.dL."      "GLC.mg.dL."       "BUN.mg.dL."       "CRE.mg.dL."       "Na.meq.L."        "K.meq.L."        
# [55] "Cl.meq.L."        "Ca.mg.dL."        "IP.mg.dL."        "TP.g.dL."         "RALB.g.dL."       "A.G"             
# [61] "AST.IU.L."        "ALT.IU.L."        "LDH.IU.L."        "GTP.IU.L."        "DNA..."           "LDH..."   
###


#note new file has different colnames so alway dbl check

test<-splitBy( ~ SPECIES + TEST_TYPE + ORGAN_ID + SIN_REP_TYPE, data=data)
names(test)
# [1] "Rat|in vivo|Liver|Single"  "Rat|in vivo|Kidney|Single" "Rat|in vivo|NA|Single"     "Rat|in vivo|Liver|Repeat" 
# [5] "Rat|in vivo|Kidney|Repeat" "Rat|in vivo|NA|Repeat"     "Rat|in vitro|Liver|NA"     "Rat|in vitro|NA|NA"       
# [9] "Human|in vitro|Liver|NA"   "Human|in vitro|NA|NA"   


#Need to get the folder names/location for each of the experiments to process the microarry data

#folder structure is SPECIES/TEST_TYPE/<ORGAN_ID>/<SIN_REP_TYPE/
#file name is cmpd.species.test_type.organ_id


#ok now for getting the file paths
#need the descriptors for the file names, do not need all the data as we can merge in the end
data2<-data[,c(1,5:13)]
#remove odd characters such that chem name can match file name
data2$COMPOUND_NAME<-gsub("[[:space:]:]", '_', data$COMPOUND_NAME)
#remove " ' " and ()for processing...i dont think there are any but to be safe
data2$COMPOUND_NAME<-gsub("['()]", '', data2$COMPOUND_NAME)
#Remoiving the space
data2$TEST_TYPE<-gsub("[[:space:]]", '_', data$TEST_TYPE)
#adding in the origional compound name for use later on
data2$CHEMICAL<-data$COMPOUND_NAME

#####################
# Folder structure and folder names are odd. Folder names are givn as :
# COMPOUND_NAME.SPECIES.TEST_TYPE.ORGAN_ID for invitro
#and COMPOUND_NAME.SPECIES.TEST_TYPE.ORGAN_ID.SIN_REP_TYPE for invivo

#need to get out the path for each microarray file. Note that some animals/cell cultures did got get arrays, thus there will be no files
#to account for that the path is set to "NA"
data2$folder<-file.path("Data/Microarray", paste(data2$COMPOUND_NAME,data2$SPECIES, data2$TEST_TYPE, data2$ORGAN_ID, data2$SIN_REP_TYPE, sep='.'), "celfiles")
data2$folder[data2$BARCODE=="No ChipData"]<-NA
data2$filename<-file.path(data2$folder, paste(data2$BARCODE, "CEL", sep='.'))
data2$filename[data2$BARCODE=="No ChipData"]<-NA

#Spot check
temp1<-subset(data2, COMPOUND_NAME == "phenylanthranilic_acid" & SIN_REP_TYPE == "Repeat"& ORGAN_ID == "Liver")
dim(temp1)
#[1] 47 13
temp1[c(1,5,22,47),]
subset(data2, CHEMICAL =="1% cholesterol + 0.25% sodium cholate")[1:2,]

#as everything looks good, going to append the origional data frame with the new identifiers of needed to get out the array data
dataNew<-cbind(data2[,c("folder", "filename")], data)
dim(dataNew)
#[1] 33566    68
colnames(dataNew)
# [1] "folder"           "filename"         "BARCODE"          "ARR_DESIGN"       "EXP_ID"           "GROUP_ID"         "INDIVIDUAL_ID"    "ORGAN_ID"        
# [9] "MATERIAL_ID"      "COMPOUND_NAME"    "COMPOUND.Abbr."   "COMPOUND_NO"      "SPECIES"          "TEST_TYPE"        "SIN_REP_TYPE"     "SEX_TYPE"        
# [17] "STRAIN_TYPE"      "ADM_ROUTE_TYPE"   "ANIMAL_AGE.week." "SACRI_PERIOD"     "DOSE"             "DOSE_UNIT"        "DOSE_LEVEL"       "TERMINAL_BW.g."  
# [25] "LIVER.g."         "KIDNEY_TOTAL.g."  "KIDNEY_R.g."      "KIDNEY_L.g."      "RBC.x10_4.ul."    "Hb.g.dL."         "Ht..."            "MCV.fL."         
# [33] "MCH.pg."          "MCHC..."          "Ret..."           "Plat.x10_4.uL."   "WBC.x10_2.uL."    "Neu..."           "Eos..."           "Bas..."          
# [41] "Mono..."          "Lym..."           "PT.s."            "APTT.s."          "Fbg.mg.dL."       "ALP.IU.L."        "TC.mg.dL."        "TG.mg.dL."       
# [49] "PL.mg.dL."        "TBIL.mg.dL."      "DBIL.mg.dL."      "GLC.mg.dL."       "BUN.mg.dL."       "CRE.mg.dL."       "Na.meq.L."        "K.meq.L."        
# [57] "Cl.meq.L."        "Ca.mg.dL."        "IP.mg.dL."        "TP.g.dL."         "RALB.g.dL."       "A.G"              "AST.IU.L."        "ALT.IU.L."       
# [65] "LDH.IU.L."        "GTP.IU.L."        "DNA..."           "LDH..."                             


write.table(dataNew, file =file.path(baseDir, "Data/Open-tggates_AllAttribute_EDT.txt") , sep='\t', col.names=TRUE, row.names=FALSE, quote=FALSE)

############
#getting out just the rat liver data
LiverData<-subset(dataNew, SPECIES=="RAT" & ORGAN_ID == 'Liver')
dim(LiverData)
write.table(LiverData, file =file.path(baseDir, "Files/Open-tggates_AllAttribute_EDT_RatLiver.txt") , sep='\t', col.names=TRUE, row.names=FALSE, quote=FALSE)


#make a smaller file for testing
#
LiverInputs<-subset(LiverData, ORGAN_ID == 'Liver' & SIN_REP_TYPE == "Single")[c(1,10,20,40,50:55),]
write.table(LiverInputs, file =file.path(baseDir, "Files/Open-tggates_AllAttribute_EDT_RatLiver_Test.txt") , sep='\t', col.names=TRUE, row.names=FALSE, quote=FALSE)




