# 1. Data Download
#1.1. Download the RNA-Seq fragements from https://www.ncbi.nlm.nih.gov/sra/?term=SRR8797509

~/Documents/NGS1_Assign./ngs1_project/Raw_Data$ wget -c ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR879/SRR8797509/SRR8797509.sra

#1.2. Convert from SRA to Fastq using fastq-dump