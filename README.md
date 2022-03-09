# Human polymerase usage sequencing (Pu-seq) R-scripts & demo data  v1.0

## Provided R-script
#### ./script/bincount-csv_to_pol-usage-wig.R
\- Generating polymerase usage data in the wig format from genome-wide bin count data.  
\- Inputting files (.count.ccv) are required to be written in the code.  

#### ./script/pol-usage-wig_to_ini-index-wig.R
\- Generating initiation index data in the wig format from polymerase usage data.  
\- Inputting files (.wig) are required to be written in the code.  

#### ./script/pol-usage-wig_to_fork-index-wig.R
\- Generating fork index data in the wig format from polymerase usage data.  
\- Inputting files (.wig) are required to be written in the code.  

#### ./script/pol-usage-wig_to_coupling-index-wig.R
\- Generating coupling index data in the wig format from polymerase usage data.  
\- Inputting files (.wig) are required to be written in the code.  

#### ./script/bincount-csv_to_pol-usage-wig.R
\- Generating coupling index data in the wig format from genome-wide bin count data.  
\- Inputting files (.count.csv) are required to be written in the cod  

'./sub' and './rcpp' contains R and rcpp source codes which are required to be implemented with script above.  
Inputting files are required to be written in the codes.  

## Provided wig format datasets
####Genome-wide bin counts datasets
./data/count.csv/pol-e-m630f.f-w1000.count.csv  
./data/count.csv/pol-e-m630f.r-w1000.count.csv  
./data/count.csv/pol-e-plus.f-w1000.count.csv  
./data/count.csv/pol-e-plus.r-w1000.count.csv  
./data/count.csv/pol-a-y865f.rep1.f-w1000.count.csv  
./data/count.csv/pol-a-y865f.rep1.r-w1000.count.csv  
./data/count.csv/pol-a-plus.f-w1000.count.csv  
./data/count.csv/pol-a-plus.r-w1000.count.csv  

#### Polymerase usage data
./data/wig/pol-a-usage.crick.w1000.MA30.rep1.wig  
./data/wig/pol-a-usage.watson.w1000.MA30.rep1.wig  
./data/wig/pol-e-usage.crick.w1000.MA30.rep1.wig  
./data/wig/pol-e-usage.watson.w1000.MA30.rep1.wig  

#### Initiation index
./data/wig/ini-index.pol-e-a.w1000.MA30-15.rep1.wig  

#### Fork index
./data/wig/fork-index.leftward.w1000.MA30-15.rep1.wig  
./data/wig/fork-index.rightward.w1000.MA30-15.rep1.wig  

#### Coupling index
./data/wig/couplig-index.leftward.w1000.MA30.rep1.wig  
./data/wig/couplig-index.rightward.w1000.MA30.rep1.wig  

#### Replication fork directionality (RFD) 
./data/wig/RFD.pol-e-a.w1000.MA30.wig  
