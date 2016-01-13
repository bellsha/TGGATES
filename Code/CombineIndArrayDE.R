#this scrip takes all the output files from AffyAnalysis.R
#and calculate the summary for community analysuis
#process: 1) standardize the controls, 2) scale each experiment by the controls
#3)remove unresponsive probes

#get the parameters from the environment:
#Parameters:
#od = the output directory for files that are to be collected
#exptinfo= a tab delim file of experimental information
#outdir = the output directory
#organ = organ name of interest
#runID = run id value
#t1 = the quantile/percential value for threshold based on expression intensity
#t2 = the quantile/percential value for threshold based on max rank change
#k = the number of communities
#note fc= fc cutoff, so fc2 is converted using log(2,2) to become 1 in the code

args=(commandArgs(TRUE))
if(length(args)==0){
  print("No arguments supplied.")
  break
  ##supply default values
  
  d <-"Z:/TGGATES2/Data/Microarray/acetaminophen.Rat.in_vivo.Liver.Single/celfiles"
  exptinfo<-"Z:/TGGATES2/Files/Open-tggates_AllAttribute_EDT_RatLiver.txt"
  basedir="Z:/TGGATES2"
  outdir <- 'Z:/TGGATES2/Output/RatLiv20160108'
  #outdir <- 'Z:/TGGATES2/Output/RatLiv20160106'
  chem<-'APAP'
  organ<-'Liver'
  runID<-"RatLiv20160108"
  #runID<-"RatLiv20160106"
  Repeat<-"MULTI"
  PathFile="Z:/TGGATES2/Files/Rat2MouseReactomeEDIT201312.txt"
  AnnFile="Z:/TGGATES2/Files/rat2302.probe.entrez.go_20141029.txt"
  t1<-0.75
#   t2<-0.75
#   k<-50
  fc<- 2
}else{
  print(args)
  for(i in 1:length(args)){
    eval(parse(text=args[[i]]))
  }
}
#libraries
#library(gdata)
library(doBy)
# library(gplots)
# library(RColorBrewer)
# mypalette <- rev(brewer.pal(9,"RdBu"))
###################################
#read in the file with the meta data for the experiments, and file names
#to make it simple as these parameters are required for AffyAnalysis to create the file name
fname<-paste("RunInfo", organ, paste(runID, 'txt', sep='.'), sep='-')
#read in the file info
info<-read.table(file.path(outdir, fname), sep='\t', header=FALSE, na.strings='', quote="\"")
colnames(info)<-c("folder","filename","BARCODE","ORGAN_ID","COMPOUND_NAME","COMPOUND.Abbr.","SPECIES","TEST_TYPE","SIN_REP_TYPE",
                  "SACRI_PERIOD","DOSE_LEVEL","Array_ID","treatment","experiment","RunID","ProcessFileInd","ProcessFileAve")
#adding in factorization of some levels
info$DOSE_LEVEL<-factor(info$DOSE_LEVEL, levels=c("Control","Low","Middle","High"), ordered=TRUE)
info$SACRI_PERIOD<-factor(info$SACRI_PERIOD, levels=c("3 hr", "6 hr", "9 hr", "24 hr", "4 day","8 day", "15 day", "29 day"), ordered=TRUE)
#make a variable for the column names on the ind processed data
info$IndColNam<-paste('X',info$Array_ID,sep='')

#output for verification:
info[1,]
#now get in all the output files from ArryAnalysis and assemble into one dataframe 
#note this is for the INDIVIDUAL treatment arrays, not the treatment averages
#due to the number of individual arrays breaking it up by treatment time
#Single dosage
filesSin<-as.character(unique(subset(info, SIN_REP_TYPE =="Single")$ProcessFileInd))
for(i in 1:length(filesSin)){
  if(i == 1){
    dataTempSin<-read.table(file.path(outdir,"ProcArray",filesSin[i]), sep='\t', header=TRUE, na.strings='', quote="\"")
  }
  else{
    inFile<-read.table(file.path(outdir,"ProcArray",filesSin[i]), sep='\t', header=TRUE, na.strings='', quote="\"")
    dataTempSin<-cbind(dataTempSin, inFile)
  }  
}
#changing class to facilitate the downstream steps
#dataTempSin<-as.matrix(dataTempSin)

# #repeat dosage
# filesRep<-as.character(unique(subset(info, SIN_REP_TYPE =="Repeat")$ProcessFileInd))
# for(i in 1:length(filesRep)){
#   if(i == 1){
#     dataTempRep<-read.table(file.path(outdir,"ProcArray",filesRep[i]), sep='\t', header=TRUE, na.strings='', quote="\"")
#   }
#   else{
#     inFile<-read.table(file.path(outdir,"ProcArray",filesRep[i]), sep='\t', header=TRUE, na.strings='', quote="\"")
#     dataTempRep<-cbind(dataTempRep, inFile)
#   }  
# }
# #changing class to facilitate the downstream steps
# dataTempRep<-as.matrix(dataTempRep)



#extract the controls and so that we use only data with controls
#NOTE that any experiments that do NOT have a matched control will be excluded!!!!!
ctrlInfo<-subset(info, DOSE_LEVEL == 'Control')
#ctrls<-ctrlInfo$treatment
treatInfo<-subset(info, DOSE_LEVEL != 'Control')
#treats<-treatInfo$treatment
validExpts<-intersect(treatInfo$experiment, ctrlInfo$experiment)
infoNew<-subset(info, experiment %in% validExpts & SIN_REP_TYPE =="Single")
#reorder data according to experimental conditions etc
dataTempSinx<-dataTempSin[,colnames(dataTempSin) %in% infoNew$IndColNam]
ctrlDataSin<-dataTempSin[, colnames(dataTempSin) %in% ctrlInfo$IndColNam]
###########
#write out the excluded experimental units
noCtrl<-setdiff(treatInfo$experiment, ctrlInfo$experiment)
if (length(noCtrl) >0){
  noCtrlData<-subset(info, experiment %in% noCtrl)
  #save output
  fname<-paste("NoControlsAVE", runID, organ, paste(gsub('-','',Sys.Date()), 'txt', sep='.'), sep='-')
  write.table(noCtrlData, file=file.path(outdir, fname), sep='\t', col.names=TRUE, quote=FALSE)
}
rm(ctrlInfo, treatInfo, noCtrlData,validExpts, inFile )
gc()
######################
#now to scale the data. Assumption is that the probe intensity distribtuion from
#control arrays is the same. using the mean value to scale
#figure out the replication factor for scaling- because individual not treatment summary do array ID
exptReps<-summaryBy(IndColNam~experiment, data=infoNew, FUN=length, keep.names=TRUE)
#summaryBy I believe does a sort, but just to be sure
exptReps<-exptReps[order(exptReps$experiment),]
#get order by EXPERIMENT so that the timing of the experiments are better ordered
#because the treatments are the column names
ctrlInfo<-subset(infoNew, DOSE_LEVEL == 'Control' & SIN_REP_TYPE =='Single')
#order the columns of the output dataframe based based on the rows of the info df
ctrlData<-dataTempSin[,as.character(ctrlInfo[do.call(order, subset(ctrlInfo, select=experiment)),]$IndColNam)]
orderAssays<-dataTempSin[,as.character(infoNew[do.call(order, subset(infoNew, select=experiment)),]$IndColNam)]

#now to get the median of the control and scaling factor
#even though we are looking at individual data, we need to develop a scaling factor based on a representative for
#each experimental unit (defined as chemical x time as the same controls are used across dose levels)
#get experimental median
ctrlExptMed<-NULL
experiments<-unique(ctrlInfo$experiment)
for(i in 1:length(experiments)){
  temp<-subset(ctrlInfo, experiment ==experiments[i])$IndColNam
  ctrlExptMed<-cbind(ctrlExptMed, apply(ctrlData[,temp],1,median ))
}


#this gets each array median
ctrlMed<-apply(as.matrix(ctrlExptMed),2,median)
#this is the median of all individuals
ctrlGMed<-median(apply(ctrlExptMed,2,as.numeric))
sf<-ctrlGMed/ctrlMed
fullsf<-rep(sf, exptReps$IndColNam)
#and scale the data
scaleArray<-sweep(orderAssays,2,FUN='*',fullsf)
dim(scaleArray)
dim(orderAssays)
summary(c(as.matrix(scaleArray)))
summary(c(as.matrix(orderAssays)))
scaleCtrl<-scaleArray[,as.character(ctrlInfo$IndColNam)]
summary(c(as.matrix(scaleCtrl)))
#scaleCtrl<-scaleCtrl[,as.character(ctrlInfo[do.call(order, subset(ctrlInfo, select=c(CHEMICAL, DOSE_LEVEL, SACRIFICE_PERIOD))),]$treatment)]
#scaleArray<-scaleArray[,as.character(infoNew[do.call(order, subset(infoNew, select=c(CHEMICAL, DOSE_LEVEL, SACRIFICE_PERIOD))),]$treatment)]

#lets write out this massive data object, since it is big using rds
fname<-paste("StdExprDataSINGLE",organ, paste(runID, 'rds', sep='.'), sep='-')
saveRDS(scaleArray, file=file.path(outdir, fname))

###################
# #plot
# #write plot to file
# pnm<-paste("DensityPlotv2",runID, paste(gsub('-','',Sys.Date()), 'png', sep='.'), sep='-')
# png(file=file.path(outdir, pnm), width=7, height=7, units="in",pointsize=16, bg='white', res=150)
# plot(density(as.matrix(orderAssays)), main='density of log2 intensities')
# lines(density(scaleArray), col=2)
# lines(density(scaleCtrl), col=3)
# legend("topright",legend=c('All data', 'Scaled data', 'Scaled controls'), col=c(1,2,3), lty=1)
# dev.off()
##################################
#Data thresholding:
#remove the probes that are consistantly low expressors
#implying that they are non-responsive to any treatment
#then remove probes that are not changing over all experiments, incd control

#first threshold
#Does a probe ever have an intensity > t1 quantile intensity of all data
#asked a different way, is the probe always in the bottom t1 quantile for expression
quantT1<-scaleArray - quantile(as.matrix(scaleArray),t1)
quantT1[quantT1 <=0]<-NA
temp<-quantT1[!apply(is.na(quantT1),1,all),]
scaleArrayT1<-scaleArray[row.names(temp),]
dim(scaleArrayT1)
summary(c(as.matrix(scaleArrayT1)))
#####
#remove probes that are nonresponsive to treatment (differentially expressed)
#get order by EXPERIMENT so that the timing of the experiments are better ordered
#because the treatments are the column names
#note the order should have been preserved, this is just insurance :)
ctrlData<-scaleArrayT1[,as.character(ctrlInfo[do.call(order, subset(ctrlInfo, select=(experiment))),]$IndColNam)]
orderAssays<-scaleArrayT1[,as.character(infoNew[do.call(order, subset(infoNew, select=(experiment))),]$IndColNam)]
#And getting full control array:
#repeat each scaling factor according to the number of individuals in an experiment
#get experimental median
ScaleCtrlExptMed<-NULL
experiments<-unique(ctrlInfo$experiment)
for(i in 1:length(experiments)){
  temp<-subset(ctrlInfo, experiment ==experiments[i])$IndColNam
  ScaleCtrlExptMed<-cbind(ScaleCtrlExptMed, apply(ctrlData[,temp],1,median ))
}
#now replicate median value of each set of controls by the number of arrays that use that control
tmp<-rep(1:ncol(ScaleCtrlExptMed), times=exptReps$IndColNam)

#now getting the corresponding columns (columns are repeated)
fullCtrl<-ScaleCtrlExptMed[,tmp]
######################
#and calculate the DE#
######################
#in log2 space so subtract
deArray<-orderAssays - fullCtrl
#remove those that are not DE based on (log2 fc <1.5)
det2<-deArray
det2[abs(deArray) < log(fc,2)] <-NA
temp<-det2[!apply(is.na(det2),1,all),]
scaleArrayT1FC<-orderAssays[row.names(temp),]
scaleArrayT1FC<-as.matrix(scaleArrayT1FC[,as.character(infoNew[do.call(order, subset(infoNew, select=treatment)),]$IndColNam)])
#
deArray<-deArray[,as.character(infoNew[do.call(order, subset(infoNew, select=treatment)),]$IndColNam)]
deData<-as.matrix(deArray[row.names(temp),colnames(deArray) %in% subset(infoNew, DOSE_LEVEL != 'Control')$IndColNam])
#save output
fname<-paste("ProcessedDatav2Single", paste("T1-", t1, sep=''),paste("FC-", fc, sep=''),organ, paste(runID, 'txt', sep='.'), sep='-')
write.table(scaleArrayT1FC, file=file.path(outdir, fname), sep='\t', row.names=TRUE, col.names=TRUE, quote=FALSE) # Writes expression values to text file in working directory.
#
fname<-paste("DiffExprDatav2Single", paste("T1-", t1, sep=''),paste("FC-", fc, sep=''),organ, paste(runID, 'txt', sep='.'), sep='-')
write.table(deData, file=file.path(outdir, fname), sep='\t', row.names=TRUE, col.names=TRUE, quote=FALSE) # Writes expression values to text file in working directory.
#generates the files needed for downstream analusis (updated 01/2015)

#descritize the data to decrease the search space
deDesc<-deData
#using a cutoff of 2 fold change
deDesc[deDesc <= -1]<--1
deDesc[deDesc > -1]<-0
deDesc[deData >= 1]<-1

fname<-paste("Desc2FCSingle",paste("T1-", t1, sep=''),paste("FC-", fc, sep=''),organ, paste(runID, 'txt', sep='.'), sep='-')
write.table(deDesc, file=file.path(outdir, fname), sep='\t', col.names=TRUE, row.names=TRUE,quote=FALSE )

# #
# ###########
# #plotting
# pnm<-paste("DataThresholdPlotv2",paste("T1-", t1, sep=''),paste("FC-", fc, sep=''),runID, paste(gsub('-','',Sys.Date()), 'png', sep='.'), sep='-')
# png(file=file.path(outdir, pnm), width=7, height=7, units="in",pointsize=16, bg='white', res=150)
# plot(density(scaleArray), main="Data Processing", col="black", xlab='log2 intensity')
# lines(density(scaleArrayT1), col="red")
# lines(density(scaleArrayT1FC),col='blue')
# legend("topright", legend=c('Inital data', paste("Post T1=",t1), paste("Post FC=",fc)), col=c('black','red','blue'), lty=1)
# dev.off()
# 
# ########
# mapkey<-as.data.frame(cbind(key=c(1:length(colnames(deData))), label=colnames(deData)))
# pnmk<-paste("heatmapKey",t1, fc,runID, paste(gsub('-','',Sys.Date()), 'txt', sep='.'), sep='-')
# write.table(mapkey, file=file.path(outdir, pnmk), sep='\t', col.names=TRUE, quote=FALSE)
# 
# dcor<-cor(t(deData), method='pearson')
# pctree<-hclust(as.dist(1-dcor), method='complete')
# pnm<-paste("DiffExprplot",runID, paste(gsub('-','',Sys.Date()), 'png', sep='.'), sep='-')
# png(file=file.path(outdir, pnm), width=8, height=8, units="in",pointsize=16, bg='white', res=150)
# heatmap.2(deData, dendrogram='row', Rowv=as.dendrogram(pctree), trace='none', col=rev(brewer.pal(9,"RdBu")), 
#           symbreaks=TRUE,labRow='', labCol=mapkey$key, main=bquote("DiffExpr,"~T1==.(t1)~FC==.(fc)))
# dev.off()
# AssayTree<-hclust(dist(t(scaleArrayT1FC), method="euclidean"), method="ward")
# ##################
# #plotting intensity heatmap
# Acor<-cor(t(scaleArrayT1FC), method='pearson')
# ProbeTree<-hclust(as.dist(1-Acor), method='complete')
# pnm<-paste("IntensityHeatmap",paste("T1-", t1, sep=''),paste("FC-", fc, sep=''),runID, paste(gsub('-','',Sys.Date()), 'png', sep='.'), sep='-')
# png(file=file.path(outdir, pnm), width=12, height=8, pointsize=14,units='in', bg='white', res=150)
# heatmap.2(as.matrix(scaleArrayT1FC), Rowv=as.dendrogram(ProbeTree),Colv='', col=mypalette[5:9], dendrogram='row', trace='none', labRow='',labCol='', main=bquote("Scaled Intensities,"~T1==.(t1)~FC==.(fc)))
# dev.off()




