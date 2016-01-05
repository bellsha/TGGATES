This set of scripts carries out the workflow for processing the TGGATEs data. 

## Data
Data was downloaded from the TGGATEs ftp site in 2013
ftp://ftp.biosciencedbc.jp/archive/open-tggates/LATEST/

## Steps
### Step1 making folders and generating the data summary file
The first step is to generate the Data Summary File using "DataFileReorg.R". 
THe output file is:Open-tggates_AllAttribute_EDT_RatLiver.txt for the rat liver data
This output contains all the experimental paramters for each animal/tissue, the clinical chemistry data, and the file path for the array data, if applicable

Before you start: DataFileReorg.R assumes the working directory has 2 folders: Data and Files. The microarray data is assumed to be unzipped in a folder labeled Microarray (Data/Microarray)
You may choose to unzip the microarray data before or after, but you will need to know where the data will be and what data you are unzipping
For the workflow described here as of 12/2015 the only microarray data being used are the Rat Liver data
To unzip (doing both Single and Repeat): 

sbell@CG34D42 //IT2/NICEATM_Data/TGGATES/ftp.biosciencedbc.jp/archive/open-tggates/LATEST/Rat/in_vivo/Liver/Single
$ unzip \*.zip -d //IT2/NICEATM_Data/TGGATES2/Data/Microarray/

sbell@CG34D42 //IT2/NICEATM_Data/TGGATES/ftp.biosciencedbc.jp/archive/open-tggates/LATEST/Rat/in_vivo/Liver/Repeat
$ unzip \*.zip -d //IT2/NICEATM_Data/TGGATES2/Data/Microarray/

### Step2 generate the YAML file
a yaml file is used in performing the batch processing on all of the array data.
A ruby script, RubyFileInput.rb, takes the output file from step 1 and creates a hash describing where data are located
for each of the "treatments"
the output file: fname = "#{outfilename}-#{Time.now.strftime('%Y%m%d')}.yml" is used in step 3

### Step3 rake task to process the  affy data
this is run by the rake file and parameters from the command line.
Calls AffyAnalysis.R and CombineArrayDE.R