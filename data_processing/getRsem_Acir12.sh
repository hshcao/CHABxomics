#!/bin/bash

#SBATCH -J Rs12Acir
#SBATCH -p kshcnormal
#SBATCH -N 2
#SBATCH --time=0
#SBATCH -n 64

conda activate rsem-env
cd /public/home/hcao/hshcao/taihu2019mGmT

for inF in ../reads/mT19_12-*1P.fq.gz
do
	inR=${inF/_1P/_2P}
	inF_base1=$(basename $inF)
	inF_base2=Acir.${inF_base1%_q26tr_1P.fq.gz}
	outF=./rsemout/Acir310F/$inF_base2
    rsem-calculate-expression -p 64 --paired-end \
                              --bowtie2 --append-names \
                              $inF $inR  rsemref/Acir310Fref $outF
done
conda deactivate
