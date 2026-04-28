#!/bin/bash
#SBATCH -J kaiju11
#SBATCH -p kshcnormal         
#SBATCH -N 1             
#SBATCH --ntasks-per-node=32  


cd /public/home/hcao/hshcao/kaijumG19/OUTequalread

kaiju2table -t /public/home/hcao/toolsNdatabase/kaiju/kaijuDB/nodes.dmp \
            -n /public/home/hcao/toolsNdatabase/kaiju/kaijuDB/names.dmp \
            -r species -l superkingdom,phylum,class,order,family,genus,species \
            -o kaiju_allSum_species.tsv mG19*_4m_kaiju.out
