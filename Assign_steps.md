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

#5.5. sequence alignment for aggressive trimmed shuffled samples to human chr 22 by using hisat 
#5.5.1. Data dwnload

~/Documents/NGS1_Assign./ngs1_project/Raw_Data/bwa_align/Index4bwa$ cd ../..

~/Documents/NGS1_Assign./ngs1_project/Raw_Data$ wget http://genomedata.org/rnaseq-tutorial/fasta/GRCh38/chr22_with_ERCC92.fa

~/Documents/NGS1_Assign./ngs1_project/Raw_Data$ wget http://genomedata.org/rnaseq-tutorial/annotations/GRCh38/chr22_with_ERCC92.gtf

#5.5.2. Indexing

~/Documents/NGS1_Assign./ngs1_project/Raw_Data$ mkdir -p hisat_align/hisatIndex && cd hisat_align/hisatIndex

~/Documents/NGS1_Assign./ngs1_project/Raw_Data/hisat_align/hisatIndex$ ln -s /home/maha/Documents/NGS1_Assign./ngs1_project/Raw_Data/chr22_with_ERCC92.fa .

~/Documents/NGS1_Assign./ngs1_project/Raw_Data/hisat_align/hisatIndex$ hisat2_extract_splice_sites.py /home/maha/Documents/NGS1_Assign./ngs1_project/Raw_Data/chr22_with_ERCC92.gtf > splicesites.tsv

~/Documents/NGS1_Assign./ngs1_project/Raw_Data/hisat_align/hisatIndex$ hisat2_extract_exons.py /home/maha/Documents/NGS1_Assign./ngs1_project/Raw_Data/chr22_with_ERCC92.gtf > exons.tsv

~/Documents/NGS1_Assign./ngs1_project/Raw_Data/hisat_align/hisatIndex$ hisat2-build -p 1 --ss splicesites.tsv --exon exons.tsv chr22_with_ERCC92.fa chr22_with_ERCC92

#5.5.3. Aligning aggressive trimmed shuffled saples to human ch22

~/Documents/NGS1_Assign./ngs1_project/Raw_Data/hisat_align/hisatIndex$ cd ..

~/Documents/NGS1_Assign./ngs1_project/Raw_Data/hisat_align/ for ((i=1;i<=5;i++)) ;

do
 
R1="/home/maha/Documents/NGS1_Assign./ngs1_project/Raw_Data/shuffled_A_Trimmed/SRR8797509_1shuffled.part_00$i.pe.trim.fq"
R2="/home/maha/Documents/NGS1_Assign./ngs1_project/Raw_Data/shuffled_A_Trimmed/SRR8797509_2shuffled.part_00$i.pe.trim.fq"
hisat2 -p 1 -x hisatIndex/chr22_with_ERCC92 --dta --rna-strandness RF -1 $R1 -2 $R2 -S shuff_SRR8797509_part_00$i.sam

done

#6. Assembly
#6.1. Assembly aligned unshuffled samples
~/Documents/NGS1_Assign./ngs1_project/Raw_Data/hisat_align/ cd ..
~/Documents/NGS1_Assign./ngs1_project/Raw_Data$ cd bwa_align

#6.1.1. Convert sam files into bam files

~/Documents/NGS1_Assign./ngs1_project/Raw_Data/bwa_align$ samtools view -bS SRR8797509_part_001.sam > SRR8797509_part_001.bam
~/Documents/NGS1_Assign./ngs1_project/Raw_Data/bwa_align$ samtools view -bS SRR8797509_part_002.sam > SRR8797509_part_002.bam
~/Documents/NGS1_Assign./ngs1_project/Raw_Data/bwa_align$ samtools view -bS SRR8797509_part_003.sam > SRR8797509_part_003.bam
~/Documents/NGS1_Assign./ngs1_project/Raw_Data/bwa_align$ samtools view -bS SRR8797509_part_004.sam > SRR8797509_part_004.bam
~/Documents/NGS1_Assign./ngs1_project/Raw_Data/bwa_align$ samtools view -bS SRR8797509_part_005.sam > SRR8797509_part_005.bam

#6.1.2. convert the BAM files to a sorted BAM files
~/Documents/NGS1_Assign./ngs1_project/Raw_Data/bwa_align$ samtools sort SRR8797509_part_001.bam -o SRR8797509_part_001.sorted.bam
~/Documents/NGS1_Assign./ngs1_project/Raw_Data/bwa_align$ samtools sort SRR8797509_part_002.bam -o SRR8797509_part_002.sorted.bam
~/Documents/NGS1_Assign./ngs1_project/Raw_Data/bwa_align$ samtools sort SRR8797509_part_003.bam -o SRR8797509_part_003.sorted.bam
~/Documents/NGS1_Assign./ngs1_project/Raw_Data/bwa_align$ samtools sort SRR8797509_part_004.bam -o SRR8797509_part_004.sorted.bam
~/Documents/NGS1_Assign./ngs1_project/Raw_Data/bwa_align$ samtools sort SRR8797509_part_005.bam -o SRR8797509_part_005.sorted.bam

#6.1.3 # Export statistical report for each sample indvidually

for ((i=1;i<=5;i++)) ;
    
do
  samtools flagstat SRR8797509_part_00$i.sorted.bam > useful_stat_$i.txt;

done

#6.1.4. Assembly without known annotations

for ((i=1;i<=5;i++)) ;
do

stringtie SRR8797509_part_00$i.sorted.bam --rf -l ref_free_$i -o ref_free_$i.gtf

done

#6.1.5. Assembly with known annotations
for ((i=1;i<=5;i++));
do
stringtie SRR8797509_part_00$i.sorted.bam --rf -l ref_sup_$i -G /home/maha/Documents/NGS1_Assign./ngs1_project/Raw_Data/chr22_with_ERCC92.gtf -o ref_sup_$i.gtf 
done

# 6.2.Assembly for aligned shuffled data

~/Documents/NGS1_Assign./ngs1_project/Raw_Data/bwa_align$ cd ..

~/Documents/NGS1_Assign./ngs1_project/Raw_Data$ cd hisat_align

# 6.2.1. convert SAM files into BAM files

~/Documents/NGS1_Assign./ngs1_project/Raw_Data/hisat_align$ samtools view -bS shuff_SRR8797509_part_001.sam > shuff_SRR8797509_part_001.bam
~/Documents/NGS1_Assign./ngs1_project/Raw_Data/hisat_align$ samtools view -bS shuff_SRR8797509_part_002.sam > shuff_SRR8797509_part_002.bam
~/Documents/NGS1_Assign./ngs1_project/Raw_Data/hisat_align$ samtools view -bS shuff_SRR8797509_part_003.sam > shuff_SRR8797509_part_003.bam
~/Documents/NGS1_Assign./ngs1_project/Raw_Data/hisat_align$ samtools view -bS shuff_SRR8797509_part_004.sam > shuff_SRR8797509_part_004.bam
~/Documents/NGS1_Assign./ngs1_project/Raw_Data/hisat_align$ samtools view -bS shuff_SRR8797509_part_005.sam > shuff_SRR8797509_part_005.bam

# 6.2.2. convert the BAM files into sorted BAM files 

~/Documents/NGS1_Assign./ngs1_project/Raw_Data/hisat_align$ samtools sort shuff_SRR8797509_part_001.bam -o shuff_SRR8797509_part_001.sorted.bam
~/Documents/NGS1_Assign./ngs1_project/Raw_Data/hisat_align$ samtools sort shuff_SRR8797509_part_002.bam -o shuff_SRR8797509_part_002.sorted.bam
~/Documents/NGS1_Assign./ngs1_project/Raw_Data/hisat_align$ samtools sort shuff_SRR8797509_part_003.bam -o shuff_SRR8797509_part_003.sorted.bam
~/Documents/NGS1_Assign./ngs1_project/Raw_Data/hisat_align$ samtools sort shuff_SRR8797509_part_004.bam -o shuff_SRR8797509_part_004.sorted.bam
~/Documents/NGS1_Assign./ngs1_project/Raw_Data/hisat_align$ samtools sort shuff_SRR8797509_part_005.bam -o shuff_SRR8797509_part_005.sorted.bam

#6.2.3. Export statistical report for each sample indvidually
for ((i=1;i<=5;i++)) ;
do
samtools flagstat shuff_SRR8797509_part_00$i.sorted.bam > shuff_useful_stat_$i.txt;
done

#6.2.3. Assembly without known annotations
for ((i=1;i<=5;i++));
do

stringtie shuff_SRR8797509_part_00$i.sorted.bam --rf -l ref_free_$i -o ref_free_$i.gtf

done

#6.2.4. Assembly with known annotations
for ((i=1;i<=5;i++)) ;
do
stringtie shuff_SRR8797509_part_00$i.sorted.bam --rf -l ref_sup_$i -G /home/maha/Documents/NGS1_Assign./ngs1_project/Raw_Data/chr22_with_ERCC92.gtf -o ref_sup_$i.gtf 
done

#7. Using GTF-Compare to Compare the Generated Annotation Files to a Reference Annotation

conda create -n ngs-gtf python=3.6 anaconda
source activate ngs-gtf
conda install -c conda-forge pypy3.5

cd ..
wget https://bootstrap.pypa.io/get-pip.py
pypy3 get-pip.py
pypy3 -m pip install gffutils numpy tqdm 'intervaltree<3.0'
mkdir -p gtf-compare/gtfs && cd gtf-compare/gtfs

ln -s /home/maha/Documents/NGS1_Assign./ngs1_project/Raw_Data/bwa_align/ref_sup_*.gtf .
ln -s /home/maha/Documents/NGS1_Assign./ngs1_project/Raw_Data/chr22_with_ERCC92.gtf .mkdir -p /home/maha/Documents/NGS1_Assign./ngs1_project/Raw_Data/gtf-compare/method_one && cd /home/maha/Documents/NGS1_Assign./ngs1_project/Raw_Data/gtf-compare/method_one

~/Documents/NGS1_Assign./ngs1_project/Raw_Data/gtf-compare/method_one$ wget https://raw.githubusercontent.com/abdelrahmanMA/gtf-compare/master/code/comp.py
~/Documents/NGS1_Assign./ngs1_project/Raw_Data/gtf-compare/method_one$ wget https://raw.githubusercontent.com/abdelrahmanMA/gtf-compare/master/code/stat.py

for ((i=1;i<=5;i++)) ;

do
pypy3 comp.py -r ../gtfs/ref_sup_*.gtf ../gtfs/chr22_with_ERCC92.gtf
pypy3 stat.py
done

#8. Differential_expression

cd ../..
mkdir -p assignment/diff_exp && cd assignment/diff_exp
mkdir assignment/ngs1_project/out && cd assignment/ngs1_project/out|mv assignment/ngs1_project/main_reads/out/*.fastq out| mv assignment/ngs1_project/shuff_reads/out/*.fastq out

#wget -c https://0x0.st/zK57.gz -O ref.tar.gz
#tar xvzf ref.tar.gz
#wget -c https://raw.githubusercontent.com/mr-eyes/nu-ngs01/master/Day-6/deseq1.r
#wget -c https://raw.githubusercontent.com/mr-eyes/nu-ngs01/master/Day-6/draw-heatmap.r#1 Setup enviornemnt
#conda activate ngs1

# conda install kallisto
# conda install samtools# Install subread, we will use featureCount : a software program developed for counting reads to genomic features such as genes, exons, promoters and genomic bins.
conda install subread# install r and dependicies

#conda install r
conda install -y bioconductor-deseq r-gplots
RUNLOG=runlog.txt#Step 2 (Quantification)Step
GTF=/home/maha/Documents/NGS1_Assign./ngs1_project/Raw_Data/chr22_with_ERCC92.gtf# Generate the counts.
featureCounts -a $GTF -g gene_name -o counts.txt  /home/maha/Documents/NGS1_Assign./ngs1_project/Raw_Data/bwa_align/main_SRR8797509*.bam  /home/maha/Documents/NGS1_Assign./ngs1_project/Raw_Data/hisat_align/shuff_SRR8797509*.bam

# Simplify the file to keep only the count columns.
cat counts.txt | cut -f 1,7-12 > simple_counts.txt

# Analyze the counts with DESeq1.
cat simple_counts.txt | Rscript deseq1.r 5x5 > results_deseq1.tsv#View only rows with pval < 0.05
cat results_deseq1.tsv | awk ' $8 < 0.05 { print $0 }' > filtered_results_deseq1.tsv
cat filtered_results_deseq1.tsv | Rscript draw-heatmap.r > hisat_output.pdf
