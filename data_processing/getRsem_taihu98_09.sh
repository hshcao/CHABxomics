#!/bin/bash

#SBATCH -J Rs09Mar
#SBATCH -p kshcnormal
#SBATCH -N 1
#SBATCH --time=0
#SBATCH -n 32

conda activate rsem-env
cd /public/home/hcao/hshcao/taihu2019mGmT

for inF in ../reads/mT19_9-10*1P.fq.gz
do
	inR=${inF/_1P/_2P}
	inF_base1=$(basename $inF)
	inF_base2=Mar.${inF_base1%_q26tr_1P.fq.gz}
	outF=./rsemout/Taihu98/$inF_base2
    rsem-calculate-expression -p 32 --paired-end \
                              --bowtie2 --append-names \
                              $inF $inR  rsemref/marTaihu98ref $outF
done
conda deactivate
