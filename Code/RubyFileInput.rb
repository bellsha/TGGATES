#this take the information from an experiment file
#and creates a list of folders to iterate through
require 'yaml'
####################
#set path, to be changed as needed
pth = "/wikihomedir/sbell/TG-Gates"

#input file containes the mappings of cel files to the appropriate folder 
infile = File::open("#{pth}/Data/Microarray_Data_RatMAPPINGS_20141030.txt", 'r')
#infile = File::open("#{pth}/Data/Microarray_Data_RatMAPPINGS_LIVER_Single_20130718.txt", 'r')
#infile = File::open("#{pth}/Data/Microarray_Data_RatMAPPINGS_TEST_20130718.txt", 'r')
header = infile.gets
hdr_cols = header.split("\t").collect {|h| h.strip}

#this section pulls out all the relevent information and puts it into a hash so it can be easily submitted as
#run parameters for R scripts
file_idx = hdr_cols.index('folder')
organ_idx = hdr_cols.index('ORGAN_ID') #Liver or Kidney
sr_idx = hdr_cols.index('SINGLE_REPEAT_TYPE') #Single, Repeat, or NA
tt_idx = hdr_cols.index('TEST_TYPE') #in_vivo or in_vitro
cmpd_idx = hdr_cols.index('COMPOUND_NAME')
organ = []
folder = []
chemical = []
wkdir = []
sr = []
tt = []
infile.each do |line|
  cols = line.split("\t")
  organ << cols[organ_idx].chomp
  folder << cols[file_idx].chomp 
  chemical << cols[cmpd_idx].chomp
  sr << cols[sr_idx].chomp
  tt << cols[tt_idx].chomp
  wkdir << "#{pth}/Data/" + "#{cols[file_idx].chomp}"
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
runP ={"exptinfo" => "#{pth}/Data/Microarray_Data_RatLabels_20141030.txt", 
  "basedir" => "#{pth}"}
moreInputs = {"runParameters" =>runP}
rInputs = moreInputs.merge({"FolderbyExpt" => inputs})
#######################################

#allInputs = {"Rparam" => rInputs}
#save to the new file
fname = "#{pth}/Files/InputParameters-#{Time.now.strftime('%Y%m%d')}.yml"

#File.open(fname,"w"){|f| f.write(allInputs.to_yaml)}
File.open(fname,"w"){|f| f.write(rInputs.to_yaml)}
#note that the working directory,chem, and organ are obtanabile from the 'folderByExpt' hash and it may be better to take it from there
puts "Parameter file #{fname} is complete."

