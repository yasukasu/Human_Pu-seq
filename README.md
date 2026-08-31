# Human polymerase usage sequencing (Pu-seq) R-scripts & demo data  v1.01

# Initial setting
\- Rcpp package has to be installed in advance to R environment.  
\- The R working directory have to be set to the directory which contain 'README.md'(./).  
\- './sub' and './rcpp' contains R and rcpp source codes which are required to be implemented with script above. For this reason, the directory structure within 'script' has to be maintained.  
\- wigToBigWig have to be executable (located in $PATH)   

## Provided R-script
#### ./script/bincount-csv_to_pol-usage-wig.R
\- Generating polymerase usage data in the wig format from genome-wide bin count data, which is produced by Perl script: pe-sam-to-bincount.pl (available at the GitHub site: https://github.com/yasukasu/sam-to-bincount).  
\- In the default setting, files in ./data/count.csv will be inputted and outputted files will appear in the working directory.  

#### ./script/pol-usage-wig_to_ini-index-wig.R
\- Generating initiation index data in the wig format from polymerase usage data.  
\- Inputting files (.wig or .bw) are required to be written in the code. In the default setting, polymerase usage data in ./data/wig will be inputted and outputted files will appear in the working directory.  

#### ./script/pol-usage-wig_to_fork-index-wig.R
\- Generating fork index data in the wig format from polymerase usage data.  
\- Inputting files (.wig or .bw) are required to be written in the code. In the default setting, polymerase usage data in ./data/wig will be inputted and outputted files will appear in the working directory.

#### ./script/pol-usage-wig_to_coupling-index-wig.R
\- Generating coupling index data in the wig format from polymerase usage data.  
\- Inputting files (.wig or .bw) are required to be written in the code. In the default setting, polymerase usage data in ./data/wig will be inputted and outputted files will appear in the working directory.

#### ./script/bincount-csv_to_pol-usage-wig.R
\- Generating coupling index data in the wig format from genome-wide bin count data.  
\- In the default setting, files in ./data/count.csv will be inputted and outputted files will appear in the working directory.   

#### ./script/pol-usage-wig_to_strand-bias-wig.R
\- Generating strand bias data in the wig format from polymerase usage data.
\- Input files (.wig or .bw) must be specified in the script. By default, polymerase-usage data in ./data/wig are used as input, and output files are written to the working directory.

#### ./script/add-RT_to_Pu-seq-intiation-zone-peaks.R
\- Assigns a replication-timing value to each Pu-seq initiation-site peak using the nearest replication-timing measurement on the same chromosome. Replication-timing data were derived from Zhao et al. (2020).
\- By default, replication-timing data in ./data/wig and initiation-site peak data in ./data/bed are used as input. Output files are written to the working directory.

#### ./script/convert-RNA-seq-FPKM-bed_to_wig.R
\- Generates a genome-wide gene-expression track in wig format from a BED file containing gene-level FPKM values.
\- By default, files in ./data/bed are used as input, and output files are written to the working directory.


## Provided wig format datasets
#### Genome-wide bin counts datasets
./data/count.csv/pol-e-m630f.f-w1000.count.csv  
./data/count.csv/pol-e-m630f.r-w1000.count.csv  
./data/count.csv/pol-e-plus.f-w1000.count.csv  
./data/count.csv/pol-e-plus.r-w1000.count.csv  
./data/count.csv/pol-a-y865f.rep1.f-w1000.count.csv  
./data/count.csv/pol-a-y865f.rep1.r-w1000.count.csv  
./data/count.csv/pol-a-plus.f-w1000.count.csv  
./data/count.csv/pol-a-plus.r-w1000.count.csv    

#### Polymerase usage data
./data/wig/pol-e-m630f_pol-usage_watson_w1000_MA30.bw  
./data/wig/pol-e-m630f_pol-usage_crick_w1000_MA30.bw  
./data/wig/pol-a-y865f_pol-usage_watson_w1000_MA30.bw  
./data/wig/pol-a-y865f_pol-usage_crick_w1000_MA30.bw  
./data/wig/pol-h-f18a_pol-usage_watson_w1000_MA30.bw
./data/wig/pol-h-f18a_pol-usage_crick_w1000_MA30.bw

#### Initiation index
./data/wig/ini-index.pol-e-a.w1000.MA30-15.rep1.wig  
./data/bed/ini-index.peak.pol-e-a.w1000.MA30-15.rep1.bedgraph

#### Fork index
./data/wig/Pol-e-a.fork-index.leftward.norm_z.w1000.MA30-15.bw  
./data/wig/Pol-e-a.fork-index.rightward.norm_z.w1000.MA30-15.bw  

#### Coupling index
./data/wig/Pol-e-a.coupling-index.rightward.w1000.MA30.bw  
./data/wig/Pol-e-a.coupling-index.leftward.w1000.MA30.bw  

#### Replication fork directionality (RFD)
./data/wig/RFD.pol-e-a.w1000.MA30.wig  

#### Replication timing (RT)
./data/wig/hct116.mean-RT.wig

#### Gene expression (FPKM)
./data/bed/hct116.FPKM.TSS-TTS.gold.bed
