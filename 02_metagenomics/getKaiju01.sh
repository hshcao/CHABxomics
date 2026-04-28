#!/bin/bash
#SBATCH -J seqtk6l          
#SBATCH -p kshcnormal         
#SBATCH -N 1             
#SBATCH --ntasks-per-node=32  

cd /public/home/hcao/hshcao/kaijumG19/OUTequalread

kaiju -z 32 -t /public/home/hcao/toolsNdatabase/kaiju/kaijuDB/nodes.dmp \
      -f /public/home/hcao/toolsNdatabase/kaiju/kaijuDB/kaiju_db_nr_euk.fmi \
      -i /public/home/hcao/hshcao/reads/mG19_6-15-1_qc4m_1P.fq.gz \
      -j /public/home/hcao/hshcao/reads/mG19_6-15-1_qc4m_2P.fq.gz \
      -o mG19_6-15-1_4m_kaiju.out
