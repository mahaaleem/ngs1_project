# 1. Data Download
#1.1. Download the RNA-Seq fragements from https://www.ncbi.nlm.nih.gov/sra/?term=SRR8797509

~/Documents/NGS1_Assign./ngs1_project/Raw_Data$ wget -c ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR879/SRR8797509/SRR8797509.sra

#1.2. Convert 5M Human RNA-Seq fragements from SRA to Fastq using fastq-dump and split R1 & R2

~/Documents/NGS1_Assign./ngs1_project/Raw_Data$ fastq-dump --split-files -X 5000000 SRR8797509

#2. Prepare the data
 
