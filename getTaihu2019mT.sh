##########################################################################################
##### FastQC
#!/bin/bash

#SBATCH -J FastqQC # job id
#SBATCH -p kshcnormal # which partition test/cpu
#SBATCH -N 1 # how many nodes
#SBATCH --time=0
#SBATCH -n 32
cd /public/home/hcao/hshcao/GSMnetwork/RNAseq/reads
module load apps/FastQC/0.11.9
for rf in aceA_1 aceB_1 eda_3 frdA_3 ppc_1 pykA_1 pykF_3 sucA_2 WT_2
do
	R1=${rf}.R1.fastq.gz
	R2=${rf}.R2.fastq.gz
	fastqc -t 32 -f fastq $R1 $R2
done
##########################################################################################
#### Trimmomatic
#!/bin/bash

#SBATCH -J trimmomatic
#SBATCH -p kshdnormal
#SBATCH -N 1
#SBATCH --time=0
#SBATCH -n 32

module load apps/Trimmomatic/0.3.8
cd /public/home/hcao/hshcao/GSMnetwork/RNAseq/reads
for rF in *.R1.fastq.gz 
do
	rR=${rF/R1/R2}
	rOut=${rF/R1.fastq.gz/q26tr}
	java -jar  /public/software/apps/Trimmomatic/0.3.8/trimmomatic-0.38.jar PE \
	-threads 32 $rF $rR -baseout $rOut HEADCROP:10 TRAILING:26 SLIDINGWINDOW:4:26 \ 
	MINLEN:115
done
##########################################################################################
############ Gene expression level estimate and differential level calculation ###########
##### RSEM with Bowtie2 #####

#!/bin/bash

#SBATCH -J rsem
#SBATCH -p kshcnormal
#SBATCH -N 1
#SBATCH --time=0
#SBATCH -n 32

conda activate rsem-env # activate rsem installed as an conda environment with bowtie2 2.5
cd /public/home/hcao/hshcao/taihu2019mGmT
### build rsem reference
mkdir rsemref
rsem-prepare-reference --gtf refgenomes/Acir310F/Acir310F_gff.gtf \
                       --bowtie2 refgenomes/Acir310F/Acir310F.fna \
                       rsemref/Acir310Fref
rsem-prepare-reference --gtf refgenomes/doliFACHB1091/doliFACHB1091_gff.gtf \
                       --bowtie2 refgenomes/doliFACHB1091/doliFACHB1091_genome.fna \
                       rsemref/dol1091ref
rsem-prepare-reference --gtf refgenomes/MARtaihu98/MARtaihu98_gff.gtf \
                       --bowtie2 refgenomes/MARtaihu98/MARtaihu98_genome.fna \
                       rsemref/marTaihu98ref                       
##### calcuate gene expression level using RSEM
#!/bin/bash

#SBATCH -J RSEM_acs
#SBATCH -p kshcnormal
#SBATCH -N 1
#SBATCH --time=0
#SBATCH -n 32
cd /public/home/hcao/hshcao/taihu2019mGmT
conda activate rsem-env
for inF in reads/ackA*_q26tr_1P.fq.gz
do
	inR=${inF/_1P/_2P}
	inF_base=${inF%_q26tr_1P.fq.gz}
	inF_out=${inF_base/reads/rsemOUT}
	rsem-calculate-expression -p 32 --paired-end \
					--bowtie2 --append-names \
					$inF $inR \
					refGenome/rsemRef/bw25113ref $inF_out
done
conda deactivate

##########
#!/bin/bash

#SBATCH -J RSEM_G6Acir
#SBATCH -p kshcnormal
#SBATCH -N 1
#SBATCH --time=0
#SBATCH -n 32

conda activate rsem-env
cd /public/home/hcao/hshcao/taihu2019mGmT

for inF in ../reads/mG19_6-15*_1P.fq.gz
do
	inR=${inF/_1P/_2P}
	inF_base1=$(basename $inF)
	inF_base2=Acir.${inF_base1%_q26tr_1P.fq.gz}
	outF=./mGrsemout/Acir310F/$inF_base2
    rsem-calculate-expression -p 32 --paired-end \
                              --bowtie2 --append-names \
                              $inF $inR  rsemref/Acir310Fref $outF
done
conda deactivate
##########################################################################################
########## read count for each clean read file
#!/bin/bash

#SBATCH -J reads_counting
#SBATCH -p kshcnormal
#SBATCH -N 1
#SBATCH --time=0
#SBATCH -n 32

cd /public/home/hcao/hshcao/taihu2019mGmT

for fq in ../reads/mG19_*1P.fq.gz; do echo $fq; zcat $fq|perl -e 'while(<>){chomp; $nn++ if /^\@/;} print $nn, "\n";' ; done
##########################################################################################
########## SNP analysis ########## ########## ########## ########## 
########## step 1: bwa assemble to .sam file
#!/bin/bash

#SBATCH -J bwa6
#SBATCH -p kshcnormal
#SBATCH -N 2
#SBATCH --time=0
#SBATCH -n 64

module load apps/bwa/0.7.17/gcc-7.3.1

cd /public/home/hcao/hshcao/taihu2019mGmT

for inF in ../reads/mG19_6-*_1P.fq.gz
do
	inR=${inF/_1P/_2P}
	inF_base1=$(basename $inF)
	inF_base2=Acir.${inF_base1%_q26tr_1P.fq.gz}.sam
	outF=./mGsnp/Acir310F/$inF_base2
    bwa mem -t 64 -M ./refgenomes/Acir310F/AcirSNPref $inF $inR > $outF
done
########## step 2: PICARD - add read groups, sorting, remove duplicates
#!/bin/bash

#SBATCH -J picard
#SBATCH -p kshcnormal
#SBATCH -N 3
#SBATCH --time=0
#SBATCH -n 96

cd /public/home/hcao/hshcao/taihu2019mGmT/mGsnp/Acir310F
module load apps/picard/2.21.9
module load apps/samtools/1.9/gcc-7.3.1
for inF in Acir.mG19_*[123].sam
do 
	outFile=${inF/\.sam/_sorted.sam}
	java -Xmx2g -jar ~/tools/picard/build/libs/picard.jar SortSam -I $inF \
	 -O $outFile -SORT_ORDER coordinate --TMP_DIR tmp
done
## draw chart
for i in *[123]_sorted.sam; do
	java -Xmx2g -jar ~/tools/picard/build/libs/picard.jar \
	 QualityScoreDistribution REFERENCE_SEQUENCE=/public/home/hcao/hshcao/taihu2019mGmT/refgenomes/Acir310F/AcirSNPref \
	  I=${i} O=${i}_Qdist.txt CHART=${i}_Qdist.pdf
done
## mark duplicates
for ii in *_sorted.sam; do
	filebase=${ii/\.sam/}
	outfile=${filebase}_rmdup.bam
	java -Xmx2g -jar ~/tools/picard/build/libs/picard.jar MarkDuplicates \
	INPUT=$ii OUTPUT=$outfile METRICS_FILE=${filebase}_duplicateMatrix REMOVE_DUPLICATES=true
done
## remove redundancy and index .bam files
for ibm in *_rmdup.bam; do 
	ofile=${ibm/\.bam/}
	java -Xmx2g -jar ~/tools/picard/build/libs/picard.jar \
	AddOrReplaceReadGroups INPUT=$ibm OUTPUT=${ofile}_addgp.bam \
	LB=whatever PL=illumina PU=whatever SM=whatever
	samtools index ${ofile}_addgp.bam
done
########## step 3: Varscan
#!/bin/bash

#SBATCH -J varscan
#SBATCH -p kshcnormal
#SBATCH -N 3
#SBATCH --time=0
#SBATCH -n 96

cd /public/home/hcao/hshcao/taihu2019mGmT/mGsnp/Acir310Fsnp

module load apps/samtools/1.9/gcc-7.3.1
conda activate bioinfo

for inF in *_sorted_rmdup_addgp.bam
do
	outfile=${inF/\.bam/_varscan.vcf}
	samtools mpileup -f ../Acir310FsnpREFgenome/Acir301FsnpRefSamIdx.fasta $inF  \
	 | varscan mpileup2snp --min-coverage 2 --min-var-freq 0.06 --output-vcf 1 --strand-filter 0 --p-value 0.10 > ./VarscanOut/$outfile
done

#####################
## parse VCF files
perl getGtf2pos_geneID.pl Acir310F_gff.gtf > Acir310F_gff_gtf_genPOS-geneID.txt
perl getSort_vcf.pl Acir310F_gff_gtf_genPOS-geneID.txt *_sorted_rmdup_addgp_varscan.vcf
perl getSnp4.pl Acir310F_gff_gtf_genPOS-geneID.txt *_sorted_rmdup_addgp_varscan.sortedVcf
perl getSNP_count_sort.pl Acir310F_gff_gtf_genPOS-geneID.txt *_varscan.snpCount