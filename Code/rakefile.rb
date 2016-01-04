#this worflow takes input from  a 'parameters' file
# and startes a set number of R process
#modifying for use in the Affy analysis

require "yaml"
require "fileutils"
#require "forkmanager"
require "fork_manager"
require "timeout"

#get inputs
#these need to be outside the task inorder for the file tasks to execute correctly
@parameters = ENV["parameters"] || "Z:/TGGATES2/Files/InputParametersTest-20160104.yml"
@OrganType = ENV["OrganType"] || "Liver"
@Cell = ENV["TestType"] || "in vivo"
@Repeat = ENV["RepeatType"] || "MULTI" #note MULTI refers to both single and repeat dosing stragaies
@runID = ENV["runID"] || "LiverRun" #note should include "single" "rep" or "invitro" in name or similar
@t1 = ENV["t1"] || 0.75
@t2 = ENV["t2"] || 0.75
#@k = ENV["k"] || 50
@fc = ENV["fc"] || 2

@config = Hash[YAML.load_file(@parameters).map{|(k,v)| [k.to_sym,v]}]
od = @config[:runParameters]["basedir"]
@outDir = "#{@config[:runParameters]["basedir"]}/Output/#{@runID}"
@logDir = "#{od}/Log/#{@runID}_#{@Repeat}"
@AffyProcDir = "#{@outDir}/ProcArray"
@AffyProcFiles = "#{@outDir}/AffyProcessedFiles-#{@OrganType}-#{@runID}.txt"
#@combfile = "#{@outDir}/ProcessedData-T1-#{@t1}-T2-#{@t2}-#{@OrganType}-#{@runID}.txt"
#@defile = "#{@outDir}/DifExpData-T1-#{@t1}-T2-#{@t2}-#{@OrganType}-#{@runID}.txt"
#@AdjMat = "#{@outDir}/AdjMatrix-T1-#{@t1}-T2-#{@t2}-#{@OrganType}-#{@runID}.RData"
@combfile = "#{@outDir}/ProcessedDatav2-T1-#{@t1}-FC-#{@fc}-#{@OrganType}-#{@runID}.txt"
@defile = "#{@outDir}/DifExpDatav2-T1-#{@t1}-FC-#{@fc}-#{@OrganType}-#{@runID}.txt"
@pathenrich = "#{@outDir}/AnnoSummaryResults-T1-#{@t1}-FC-#{@fc}-#{@OrganType}-#{@runID}.txt"
@AdjMat = "#{od}/Output/AdjMatrix-T1-#{@t1}-T2-#{@t2}-#{@OrganType}-#{@runID}.txt"


@AnnFile = ENV["AnnFile"] || "#{od}/Files/rat2302.probe.entrez.go_20141029.txt"
@PathFile = ENV["PathFile"] || "#{od}/Files/Rat2MouseReactomeEDIT201312.txt"

@DefCom = "#{@outDir}/CommunityAssignment-T1-#{@t1}-T2-#{@t2}-k-#{@k}-#{@OrganType}-#{@runID}.txt"
@GO_db = "#{od}/GO/GO_rat2302"
#@GO_db = "#{@outDir}/testGO"
#@go_df = "#{@outDir}/GO_rat2302/go_df"
#@go_df = "#{@GO_db}/go_df"
@Enrich = "#{@outDir}/AdjNetComGrp-T1-#{@t1}-T2-#{@t2}-k-#{@k}-#{@OrganType}-#{@runID}.txt"
@EucComm = "#{@outDir}/EUCEnruch-T1-#{@t1}-T2-#{@t2}-k-#{@k}-#{@OrganType}-#{@runID}.txt"
@DeCom = "#{@outDir}/DeCommunities-T1-#{@t1}-T2-#{@t2}-k-#{@k}-#{@OrganType}-#{@runID}.RData"
############################################
##rake parameters="../Files/TestInputParameters-20130711.yml" Go_database --trace
##make the go files and the database
#directory "#{@GO_db}"
#file "#{@go_df}" => "#{@GO_db}" do
#  system("R CMD BATCH --no-save --no-restore '--args GODir=#{@GO_db.inspect} annFile=#{@AnnFile.inspect}' #GoMapping.R #{@logDir}/RGOMAP-#{Time.now.strftime('%Y_%m_%d')}.out")
#end
#desc "Creates the GO database for use in enrichment mapping"
#task :Go_database => "#{@go_df}" do
#  puts "Creating GO reference databases. See folder in #{@GO_db} and log in #{@logDir}/RGOMAP-#{Time.now.strftime('%Y_%m_%d')}.out"
#end
####################################################
#to run:
# rake parameters="Z:/TGGATES2/Files/InputParametersRatLiver-20160104.yml" OrganType="Liver" runID="testLivSviv"  get_input --dry-run
#desc "Prep parameters"
task :get_input do |t, args|
  #getting the parameters in a way that can be used by R scripts
  #so extract values from the hash table and return them as 'key=value'
  #here I am making a string that includes all the parameters needed to execute the R scripts
  #which are collected in an array that is accessible from anywhere in the rakefile 
  #note that I need: d=directory, exptinfo= experimental info, chem=chemical, organ=organ/cell type
  #get out the run parameters into a string
  rpStr = ''
  @config[:runParameters].each do |k,v|
    rpStr.concat("#{k}=#{v.inspect} ")
  end
  @pAry = []
  if @Repeat !='MULTI'
    folders = @config[:FolderbyExpt][@OrganType][@Cell][@Repeat]
  else
    folders = @config[:FolderbyExpt][@OrganType][@Cell]['Single']
    tmp = @config[:FolderbyExpt][@OrganType][@Cell]['Repeat']
    folders = folders.merge(tmp)
  end
  kFol = folders.keys
  nFol = kFol.length
  for j in 0..nFol-1 do
    pStr = ''
    folders[kFol[j]].each do |k,v|
      #regular expression check to add quotes to character
      #this is needed for R
      pStr.concat("#{k}=#{v.inspect} ")
    end
    #note that there needs to be a space b/w variables!!!
#    @pAry[j] = pStr + rpStr + "t1=#{@t1} "+"t2=#{@t2} "+"fc=#{@fc} "+"k=#{@k} "+"runID=#{@runID.inspect} "+"outdir=#{@outDir.inspect} "+"Repeat=#{@Repeat.inspect} "+"AnnFile=#{@AnnFile.inspect} "+"PathFile=#{@PathFile.inspect}"
	@pAry[j] = pStr + rpStr + "t1=#{@t1} "+"t2=#{@t2} "+"fc=#{@fc} "+"runID=#{@runID.inspect} "+"outdir=#{@outDir.inspect} "+"Repeat=#{@Repeat.inspect} "+"AnnFile=#{@AnnFile.inspect} "+"PathFile=#{@PathFile.inspect}"
  end
#   puts "Input Complete!!"
   puts "Used parameters located #{@parameters}"
   puts "For tissue type #{@OrganType}"
end

###########################################################

########
#this file task will generate the input for the R_Affy, but still allow R_Affy to be called
#note that this will always generate a new file (or rewrite an existing one), 
#but that file will have the files uesed for that ru    @pAry[j] = pStr + rpStr + "t1=#{@t1} "+"t2=#{@t2} "+"k=#{@k} "+"runID=#{@runID.inspect} "+"outdir=#{@outDir.inspect} "+"Repeat=#{@Repeat.inspect} "nID
desc "RMA normalization by folder"
file "#{@AffyProcFiles}" => ["#{@AffyProcDir}", "#{@outDir}", "#{@logDir}"] do
#file "#{@AffyProcFiles}" =>[:get_input] do |t, args|
  #write a file now for file task
  od=@outDir
  chemList = []
  if @Repeat !='MULTI'
    folders = @config[:FolderbyExpt][@OrganType][@Cell][@Repeat]
  else
    folders = @config[:FolderbyExpt][@OrganType][@Cell]['Single']
    tmp = @config[:FolderbyExpt][@OrganType][@Cell]['Repeat']
    folders = folders.merge(tmp)
  end
#  temp = @config[:FolderbyExpt][@OrganType][@Cell][@Repeat]
  tk = folders.keys
  for i in 0..tk.length-1 do
    cpd = folders[tk[i]]["chem"]
    chemList << "ExpAVE-" + "#{@OrganType}-" "#{cpd}" + "-#{@runID}.txt"
  end
  test = Dir["#{od}/ProcArray/ExpAVE-#{@OrganType}*-#{@runID}.txt"].collect{|f| f.gsub("#{od}/ProcArray",'')}

  if chemList.all?{|c| test.include?(c)}
    puts "Data are already pre-processed"
  else 
    #add in functionality to do those not on the list???
    numProc = ENV["numProc"].to_i || 4
    fm = ForkManager.new(numProc)
    stats = fm.manage do
      n = @pAry.length
      for i in 0..n-1 do
       # puts "Starting iteration #{i}"
        fm.fork do
          ccem=@pAry[i].match(/chem=\"(\S+)\"/)[1]        
          system("R CMD BATCH --no-save --no-restore '--args #{@pAry[i]}' AffyAnalysis.R #{@logDir}/RAffy-#{ccem}-#{@runID}.out")
        # puts " CHILD PID #{Process.pid} for #{ccem} done!"
        end
      end
    end
  end  
  filesEnd = Dir["#{od}/ProcArray/ExpAVE-#{@OrganType}*-#{@runID}.txt"]
  File.open("#{@AffyProcFiles}", 'w') {|f| f.puts(filesEnd) }  
end
###
directory "#{@AffyProcDir}"
directory "#{@outDir}"
directory "#{@logDir}"
##############################
#to run: (use -n for a dryrun first)
# rake parameters="Z:/TGGATES2/Files/InputParametersRatLiver-20160104.yml" OrganType="Liver" runID="RatLiv20160104" numProc=3  R_Affy --dry-run

desc "Perform RMA normilzation and treatment-level summerization for each chemical by calling a file task"
task :R_Affy =>[:get_input, "#{@AffyProcFiles}"] do
   puts "Processing files complete"
end 
########################################################
##UPDATE!!!!
########################################################
#to run
#rake parameters="../Files/InputParameters-20141030.yml" OrganType="Liver" runID="NetRun1410" numProc=15 t1=0.75 fc=2 R_Combine --trace

desc "Combine expression dataframes from the individual experiments"
file "#{@combfile}" =>["#{@AffyProcFiles}"] do
  n = @pAry.length
  #note that the extra parameters in @pAry wont matter here
  #they are chem, d, exptinfo. chem and d are the only ones that change w/iterations
  system("R CMD BATCH --no-save --no-restore '--args #{@pAry[1]}' CombineArrayDE.R #{@logDir}/RComb-#{@runID}.out")
  puts "Data Aggregation Finished :) "
  puts "see log file #{@logDir}/RComb-#{@runID}.out"
  puts "#{@combfile}"
end
#############################3
desc "Prepare the the Combined expression dataframe"
task :R_Combine  =>[:get_input,"#{@combfile}"] do
  puts "Combination dataframe ready"
  puts "#{@combfile}"
end
##########################################
#################################################33
#########Updated workflow
##############################################
#to run
#rake parameters="../Files/InputParameters-20141030.yml" OrganType="Liver" runID="NetRun1410" numProc=20 t1=0 fc=2 R_PathEnrich --trace

desc "Calculate the pathway enrichment for various clustering approaches"
file "#{@pathenrich}" =>["#{@combfile}"] do
  n = @pAry.length
  #note that the extra parameters in @pAry wont matter here
  #they are chem, d, exptinfo. chem and d are the only ones that change w/iterations
  system("R CMD BATCH --no-save --no-restore '--args #{@pAry[1]}' EnrichTestPathAuto2.R #{@logDir}/RpEnrich-#{@runID}.out")
  puts "Enrichment Calcs Finished :) "
  puts "see log file #{@logDir}/RpEnrich-#{@runID}.out"
  puts "#{@pathenrich}"
end
#############################################

desc "Calculate the pathway enrichment"
task :R_PathEnrich  =>[:get_input,"#{@pathenrich}"] do
  puts "Enrichment object and summery ready"
  puts "#{@pathenrich}"
end


#######################################################################
#rake parameters="../Files/TestInputParameters-20130718.yml" OrganType="Liver" runID="testLiv" RepeatType="Single" numProc=4 t1=0.75 t2=0.75 R_DiffExpr --trace
desc "Calculate the differential expression (log2) of experiments"
file "#{@defile}" =>["#{@combfile}"] do
  n = @pAry.length
  #note that the extra parameters in @pAry wont matter here
  #they are chem, d, exptinfo. chem and d are the only ones that change w/iterations
  system("R CMD BATCH --no-save --no-restore '--args #{@pAry[1]}' DiffExp.R #{@logDir}/RDiffE-#{@t1}-#{@t2}-#{@runID}.out")
  puts "Differential Expression DF generated! "
  puts "see log file #{@logDir}/RDiffE-#{@runID}.out"
end
#############################################
desc "Calculate the DE of the combined arrays"
task :R_DiffExpr  =>[:get_input,"#{@defile}"] do
  puts "DiffExp dataframe ready"
end

#############################################
#######################################################################
#rake parameters="../Files/TestInputParameters-20130718.yml" OrganType="Liver" runID="testLiv" RepeatType="Single" numProc=4 t1=0.75 t2=0.75 k=50 R_AdjMat --trace
desc "Calculate the adjaceny matrices for combined data (log2 intensity)"
file "#{@AdjMat}" =>["#{@combfile}"] do
  n = @pAry.length
  #note that the extra parameters in @pAry wont matter here
  system("R CMD BATCH --no-save --no-restore '--args #{@pAry[1]}' CalcAdj.r #{@logDir}/CalcAdj-#{@t1}-#{@t2}-#{@runID}.out")
  puts "Adjacency matrices generated! "
  puts "see log file #{@logDir}/CalcAdj-#{@runID}.out"
end
#############################################
desc "Calculate the adjancy matrices for the combined arrays"
task :R_AdjMat  =>[:get_input,"#{@AdjMat}"] do
  puts "Adjacency matrices are ready"
end
###############################################################################################
#######################################################################
#rake parameters="../Files/TestInputParameters-20130718.yml" OrganType="Liver" runID="testLiv" RepeatType="Single" numProc=4 t1=0.75 t2=0.75 k=50 R_DefCom --trace

#AnnFile="/datadrive/GitProjects/DrugMatrix/Files/rat2303.probe.entrez.go130513.txt"
#rake parameters="../Files/TestInputParameters-20130711.yml" OrganType="LIVER" runID="testRun" t1=0.75 t2=0.75 k=50 AnnFile="/datadrive/GitProjects/DrugMatrix/Files/rat2303.probe.entrez.go130513.txt" R_DefCom --trace
desc "Probe mapping and community definition"
file "#{@DefCom}" =>["#{@AdjMat}"] do
  #note that the extra parameters in @pAry wont matter here
  #they are chem, d, exptinfo. chem and d are the only ones that change w/iterations
  args2 = "#{@pAry[1]} AnnFile=#{@AnnFile.inspect} "
  system("R CMD BATCH --no-save --no-restore '--args #{args2}' DefCommunity.r #{@logDir}/DefCom-#{@t1}-#{@t2}-#{@k}-#{@runID}.out")
  puts "Communities defined and probe ID's mapped using #{@AnnFile}"
  puts "see log file #{@logDir}/DefCom-#{@runID}.out"
end
#############################################
desc "Calculate the communities from adjancy matrices and map the probes to ID's"
task :R_DefCom  =>[:get_input,"#{@DefCom}"] do
  puts "Communities and mapping ready"
end
###############################################################################################
#######################################################################
#rake parameters="../Files/TestInputParameters-20130718.yml" OrganType="Liver" runID="testLiv" RepeatType="Single" numProc=4 t1=0.75 t2=0.75 k=50 R_Enrich --trace
desc "GO enrichment"
file "#{@Enrich}" =>["#{@DefCom}", "#{@go_df}"] do
  #note that the extra parameters in @pAry wont matter here
  #they are chem, d, exptinfo. chem and d are the only ones that change w/iterations
  system("R CMD BATCH --no-save --no-restore '--args #{@pAry[1]} GODir=#{@GO_db.inspect}' GoEnrich.R #{@logDir}/GoEnrich-#{@t1}-#{@t2}-#{@k}-#{@runID}.out")
  puts "Communities defined and probe ID's mapped! "
  puts "see log file #{@logDir}/GoEnrich-#{@runID}.out"
end
#############################################
desc "Calculates the community enrichment (based on GO)"
task :R_Enrich => [:get_input, "#{@Enrich}"] do 
  #note that the extra parameters in @pAry wont matter here
  #they are chem, d, exptinfo. chem and d are the only ones that change w/iterations
  #all parameters are in the script
  puts "Community Analysis Finished :) "
end
############################################################################################
#######################################################################
#rake parameters="../Files/TestInputParameters-20130718.yml" OrganType="Liver" runID="testLiv" RepeatType="Single" numProc=4 t1=0.75 t2=0.75 k=50 R_Enrich --trace
desc "Euc Communities"
file "#{@EucComm}" =>["#{@AdjMat}", "#{@go_df}"] do
  #note that the extra parameters in @pAry wont matter here
  #they are chem, d, exptinfo. chem and d are the only ones that change w/iterations
  system("R CMD BATCH --no-save --no-restore '--args #{@pAry[1]} GODir=#{@GO_db.inspect}' EucComm.R #{@logDir}/EucComm-#{@t1}-#{@t2}-#{@k}-#{@runID}.out")
  puts "Communities defined and probe ID's mapped! "
  puts "see log file #{@logDir}/EucComm-#{@runID}.out"
end
#############################################
desc "Calculates the community enrichment (based on GO)"
task :R_EucComm => [:get_input, "#{@EucComm}"] do 
  #note that the extra parameters in @pAry wont matter here
  #they are chem, d, exptinfo. chem and d are the only ones that change w/iterations
  #all parameters are in the script
  puts "AdjMat through Go analysis via Eucledian distance Finished :) "
end
############################################################################################

##################################
desc "DE Communities"
file "#{@DeCom}" =>["#{@Enrich}","#{@defile}"] do
  system("R CMD BATCH --no-save --no-restore '--args #{@pAry[1]}' DEcomm.R #{@logDir}/DEComm-#{@runID}.out")
  puts "DE by Communities defined"
  puts "see log file #{@logDir}/DEComm-#{@runID}.out"
end
##################
desc "Calculates the average DE for each community"
task :R_DeCom => [:get_input, "#{@DeCom}"] do 
  puts "DE by Community Finished :) "
end
#########################################################3
#rake parameters="../Files/InputParameters-20130718.yml" OrganType="Liver" runID="LivS" RepeatType="Single" numProc=6 t1=0.75 t2=0.75 k=50 DEGOworkFlow --trace
desc "calculates the DE values as well as GO enrichment on intensities"
multitask :DEGOworkFlow => [:R_AdjMat, :R_DeCom] do
   puts "parallel execution of DE and adjmat on intensity"
end

###############################################################
#This is a little chunk to run the NetworkBreakTest.R script
desc "run network breakdown testing script"
task :NetworkBreakTest do
    fm = ForkManager.new(4)
    runs = ['chem','chemS','chemR','mech']
    stats = fm.manage do
      for i in 0..3 do
       # puts "Starting iteration #{i}"
        fm.fork do
          system("R CMD BATCH --no-save --no-restore '--args type=#{runs[i].inspect}' NetworkBreakTest.R NetworkBreak-#{runs[i].inspect}.out")
        end
      end
    end
end
    
    
    
    
    
