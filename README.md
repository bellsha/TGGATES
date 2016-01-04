This set of scripts carries out the workflow for processing the TGGATEs data. 

## Data
Data was downloaded from the TGGATEs ftp site in 2013
ftp://ftp.biosciencedbc.jp/archive/open-tggates/LATEST/

## Steps
The first step is to generate the Data Summary File using "DataFileReorg.R". 
Before you start: DataFileReorg.R assumes the working directory has 2 folders: Data and Files. The microarray data is assumed to be unzipped in a folder labeled Microarray (Data/Microarray)
You may choose to unzip the microarray data before or after, but you will need to know where the data will be and what data you are unzipping
For the workflow described here as of 12/2015 the only microarray data being used are the Rat Liver data
To unzip (doing both Single and Repeat): 

sbell@CG34D42 //IT2/NICEATM_Data/TGGATES/ftp.biosciencedbc.jp/archive/open-tggates/LATEST/Rat/in_vivo/Liver/Single
$ unzip \*.zip -d //IT2/NICEATM_Data/TGGATES2/Data/Microarray/

sbell@CG34D42 //IT2/NICEATM_Data/TGGATES/ftp.biosciencedbc.jp/archive/open-tggates/LATEST/Rat/in_vivo/Liver/Repeat
$ unzip \*.zip -d //IT2/NICEATM_Data/TGGATES2/Data/Microarray/
