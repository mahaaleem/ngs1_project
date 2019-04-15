# 1. Data Download
#1.1. Download the RNA-Seq fragements from https://www.ncbi.nlm.nih.gov/sra/?term=SRR8797509

mkdir Raw_Data && cd Raw_Data
wget -c ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR879/SRR8797509/SRR8797509.sra

#1.2. Convert 5M Human RNA-Seq fragements from SRA to Fastq using fastq-dump and split R1 & R2
source activate ngs1
fastq-dump --split-files -X 5000000 SRR8797509

#2. Prepare the data

seqkit split2 -1 SRR8797509_1.fastq -2 SRR8797509_2.fastq -p 5 -O notshuffled_splitted -f

seqkit shuffle SRR8797509_1.fastq > SRR8797509_1shuffled.fastq

seqkit shuffle SRR8797509_2.fastq > SRR8797509_2shuffled.fastq

seqkit split2 -1 SRR8797509_1shuffled.fastq -2 SRR8797509_2shuffled.fastq -p 5 -O shuffled_splitted -f

#3. FASTQ Quality Control

mkdir all_S1 && cd notshuffled_splitted

cp {SRR8797509_1.part_001.fastq,SRR8797509_2.part_001.fastq} ../all_S1

cd ..

cd shuffled_splitted

cp {SRR8797509_1shuffled.part_001.fastq,SRR8797509_2shuffled.part_001.fastq} ../all_S1

cd all_S1

for f in  *.fastq  ;

do

fastqc -t 1 -f fastq -noextract $f;

done

~/Documents/NGS1_Assign./ngs1_project/Raw_Data/all_S1$ multiqc -z -o . .

#4. Trimming
#4.1. Mild Trimming for unshufflled samples

mkdir notshuffled_Trimmed

for ((i=1;i<=5;i++)) ; 

do

f1="notshuffled_splitted/SRR8797509_1.part_00$i.fastq"
f2="notshuffled_splitted/SRR8797509_2.part_00$i.fastq"

newf1="notshuffled_Trimmed/SRR8797509_1.part_00$i.pe.trim.fq"
newf2="notshuffled_Trimmed/SRR8797509_2.part_00$i.pe.trim.fq"

newf1U="notshuffled_Trimmed/SRR8797509_1.part_00$i.se.trim.fq"
newf2U="notshuffled_Trimmed/SRR8797509_2.part_00$i.se.trim.fq"

adap="/home/maha/miniconda3/envs/ngs1/share/trimmomatic-0.38-1/adapters"

trimmomatic PE -threads 1 -phred33 -trimlog trimLogFile -summary statsSummaryFile  $f1 $f2 $newf1 $newf1U $newf2 $newf2U \
ILLUMINACLIP:$adap/TruSeq3-PE.fa:2:30:10:1 SLIDINGWINDOW:4:15 MINLEN:36

done

#4.2. Aggressive Trimming for shufflled samples

mkdir shuffled_A_Trimmed

for ((i=1;i<=5;i++)) ; 

do

f1="shuffled_splitted/SRR8797509_1shuffled.part_00$i.fastq"
f2="shuffled_splitted/SRR8797509_2shuffled.part_00$i.fastq"

newf1="shuffled_A_Trimmed/SRR8797509_1shuffled.part_00$i.pe.trim.fq"
newf2="shuffled_A_Trimmed/SRR8797509_2shuffled.part_00$i.pe.trim.fq"

newf1U="shuffled_A_Trimmed/SRR8797509_1shuffled.part_00$i.se.trim.fq"
newf2U="shuffled_A_Trimmed/SRR8797509_2shuffled.part_00$i.se.trim.fq"

adap="/home/maha/miniconda3/envs/ngs1/share/trimmomatic-0.38-1/adapters"

trimmomatic PE -threads 1 -phred33 -trimlog trimLogFile -summary statsSummaryFile  $f1 $f2 $newf1 $newf1U $newf2 $newf2U \
ILLUMINACLIP:$adap/TruSeq3-PE.fa:2:30:10:1 SLIDINGWINDOW:4:25 MINLEN:36

done

#5- Alignment [Align all the samples (1:5) using BWA and Hisat against the human reference file. (BWA for unshuffled and HISAT for huffled)]

#5.1. download gencode.v29.pc_transcripts & gencode.v29.annotation

wget -c ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_29/gencode.v29.pc_transcripts.fa.gz

gunzip gencode.v29.pc_transcripts.fa.gz

wget -c ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_29/gencode.v29.annotation.gtf.gz
gunzip gencode.v29.annotation.gtf.gz

#5.2. Select the transcripts of Chr22

READS=$(grep "^chr22" gencode.v29.annotation.gtf | awk -F'\t' '{print $9}' | awk -F';' '{print $1}' | awk -F' ' '{print $2}' | awk -F'"' '{print $2}' | sort | uniq)

for value in $READS ;
    
    do 
        echo "Processing: $value"
        seqkit grep -r -p ${value} gencode.v29.pc_transcripts.fa | awk -F'|' '{print $1}' >> gencode.v29.pc_transcripts.chr22.simplified.fa
    
    done
    
 #5.3. Indexing gencode.v29.pc_transcripts.chr22.simplified.fa

mkdir -p bwa_align/bwaIndex && cd bwa_align/bwaIndex

ln -s  ln -s /home/maha/Documents/NGS1_Assign./ngs1_project/Raw_Data/gencode.v29.pc_transcripts.chr22.simplified.fa .
 . 

bwa index -a bwtsw gencode.v29.pc_transcripts.chr22.simplified.fa

#5.4. sequence alignment for mild trimmed unshuffled samples to human chr 22 by using bwa

~/Documents/NGS1_Assign./ngs1_project/Raw_Data/bwa_align/Index4bwa$ cd ..
~/Documents/NGS1_Assign./ngs1_project/Raw_Data/bwa_align$ for ((i=1;i<=5;i++)) ; 
 do
R1="/home/maha/Documents/NGS1_Assign./ngs1_project/Raw_Data/notshuffled_Trimmed/SRR8797509_1.part_00$i.pe.trim.fq"
R2="/home/maha/Documents/NGS1_Assign./ngs1_project/Raw_Data/notshuffled_Trimmed/SRR8797509_2.part_00$i.pe.trim.fq"
/usr/bin/time -v bwa mem bwaIndex/gencode.v29.pc_transcripts.chr22.simplified.fa $R1 $R2 > SRR8797509_part_00$i.sam
done



