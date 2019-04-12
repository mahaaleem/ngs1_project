# 1. Data Download
#1.1. Download the RNA-Seq fragements from https://www.ncbi.nlm.nih.gov/sra/?term=SRR8797509

~/Documents/NGS1_Assign./ngs1_project/Raw_Data$ wget -c ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR879/SRR8797509/SRR8797509.sra

#1.2. Convert 5M Human RNA-Seq fragements from SRA to Fastq using fastq-dump and split R1 & R2

~/Documents/NGS1_Assign./ngs1_project/Raw_Data$ fastq-dump --split-files -X 5000000 SRR8797509

#2. Prepare the data
~/Documents/NGS1_Assign./ngs1_project/Raw_Data$ seqkit split2 -1 SRR8797509_1.fastq -2 SRR8797509_2.fastq -p 5 -O notshuffled_splitted -f

~/Documents/NGS1_Assign./ngs1_project/Raw_Data$ seqkit shuffle SRR8797509_1.fastq > SRR8797509_1shuffled.fastq

~/Documents/NGS1_Assign./ngs1_project/Raw_Data$ seqkit shuffle SRR8797509_2.fastq > SRR8797509_2shuffled.fastq

~/Documents/NGS1_Assign./ngs1_project/Raw_Data$ seqkit split2 -1 SRR8797509_1shuffled.fastq -2 SRR8797509_2shuffled.fastq -p 5 -O shuffled_splitted -f 

#3. FASTQ Quality Control
~/Documents/NGS1_Assign./ngs1_project/Raw_Data$ mkdir all_S1 && cd notshuffled_splitted

~/Documents/NGS1_Assign./ngs1_project/Raw_Data/notshuffled_splitted$ cp {SRR8797509_1.part_001.fastq,SRR8797509_2.part_001.fastq} ../all_S1

~/Documents/NGS1_Assign./ngs1_project/Raw_Data/notshuffled_splitted$ cd ..

~/Documents/NGS1_Assign./ngs1_project/Raw_Data$ cd shuffled_splitted

 cp {SRR8797509_1shuffled.part_001.fastq,SRR8797509_2shuffled.part_001.fastq} ../all_S1

~/Documents/NGS1_Assign./ngs1_project/Raw_Data$ cd all_S1
