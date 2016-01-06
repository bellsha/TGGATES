#this take the information from an experiment file
#and creates a list of folders to iterate through
require 'yaml'
####################
#set path, to be changed as needed
#pth = "Z:/IT2/NICEATM_Data/TGGATES2"
#local
pth = "Z:/TGGATES2"
#change this as needed
#test
#infilename="/Files/Open-tggates_AllAttribute_EDT_RatLiver_Test.txt"
#outfilename="/Files/InputParametersTest"
#all the microarray
#infilename = "/Files/Open-tggates_AllAttribute_EDT.txt"
#rat liver data (does include in vivo)
infilename = "/Files/Open-tggates_AllAttribute_EDT_RatLiver.txt"
outfilename="/Files/InputParametersRatLiver"

infile = File::open("#{pth}#{infilename}", 'r')


header = infile.gets
hdr_cols = header.split("\t").collect {|h| h.strip}

#this section pulls out all the relevent information and puts it into a hash so it can be easily submitted as
#run parameters for R scripts
file_idx = hdr_cols.index('folder')
organ_idx = hdr_cols.index('ORGAN_ID') #Liver or Kidney
sr_idx = hdr_cols.index('SIN_REP_TYPE') #Single, Repeat, or NA
tt_idx = hdr_cols.index('TEST_TYPE') #in_vivo or in_vitro
cmpd_idx = hdr_cols.index('COMPOUND.Abbr.')
organ = []
folder = []
chemical = []
wkdir = []
sr = []
tt = []
infile.each do |line|
#note the scrub us to take care of some utf-8 issues
  cols = line.scrub.split("\t")
  organ << cols[organ_idx].chomp
  folder << cols[file_idx].chomp 
  chemical << cols[cmpd_idx].chomp
  sr << cols[sr_idx].chomp
  tt << cols[tt_idx].chomp
  wkdir << "#{pth}" + "/#{cols[file_idx].chomp}"
end
#files = folder

inputs = Hash.new{|h,k| h[k]=Hash.new(&h.default_proc) }
folder.each_with_index do |f,i|
   #this creates a nested hash where the top level is the organ inputs[organ]
   #the second level is the folder inputs[organ][folder]
   #and from that you get the d/chem/organ bit
   inputs[organ[i]][tt[i]][sr[i]][folder[i]] = {"d" => wkdir[i], "chem" => chemical[i], "organ" =>organ[i]}
  
end
puts "Organ keys in input file:"
puts inputs.keys

###############################
#experiment specific...these are the ones you SHOULD change
runP ={"exptinfo" => "#{pth}#{infilename}", 
  "basedir" => "#{pth}"}
moreInputs = {"runParameters" =>runP}
rInputs = moreInputs.merge({"FolderbyExpt" => inputs})
#######################################

#save to the new file
#write the test
#fname = "#{pth}/Files/InputParametersTest-#{Time.now.strftime('%Y%m%d')}.yml"
#write the actual
fname = "#{pth}#{outfilename}-#{Time.now.strftime('%Y%m%d')}.yml"
File.open(fname,"w"){|f| f.write(rInputs.to_yaml)}
#note that the working directory,chem, and organ are obtanabile from the 'folderByExpt' hash and it may be better to take it from there
puts "Parameter file #{fname} is complete."

