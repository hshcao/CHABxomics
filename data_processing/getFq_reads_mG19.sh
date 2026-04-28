#!/bin/bash

#SBATCH -J reads_mGcounting
#SBATCH -p kshcnormal
#SBATCH -N 1
#SBATCH --time=0
#SBATCH -n 32

cd /public/home/hcao/hshcao/taihu2019mGmT

for fq in ../reads/mG19_*2P.fq.gz
do
	Total_lines=`zcat $fq|wc -l`
	fq_base=${basename $fq)
	fq_clean=${fq_base%_q26tr_2P.fq.gz}
	outputL=${fq_clean}."ssst".${Total_lines}
	echo $outputL >> mG19_allFiles_readCount.txt
done
