# Run with SLURM on IFB core cluster

# Script to run MAFFT which will perform multiple sequence alignment

analysis="mafft_sirtuins_uniprot"
threads=8
memory=64GB
partition="fast"
database="/shared/home/nnalpas/BACTpredict/Sirtuin_conservation/uniref_SIRTs_Sir2_CobB_SirTM_2024_06_25_selected2.fasta"
out="`echo ${database} | sed -e 's/.fasta/_msa/'`"
mafft_options="--amino --auto --treeout --reorder"
#gblocks_options="-t=p -b1=0.3 -b4=5 -s=y -d=y -e=_gbl"
trimal_options="-keepheader -gt 0.8 -cons 1 -w 3 -htmlout ${out}_summary.html"
fasttree_options=""
#iqtree_options="--seqtype AA -m TIM2+I+G -B 1000 --prefix iq_"

basedir=`dirname $database`
mydate=`date +%Y-%m-%d_%H-%M-%S`

echo '#!/bin/bash' > ${basedir}/${analysis}_${mydate}.sh
echo '#' >> ${basedir}/${analysis}_${mydate}.sh
echo "#SBATCH -p ${partition}" >> ${basedir}/${analysis}_${mydate}.sh
echo "#SBATCH --cpus-per-task=${threads}" >> ${basedir}/${analysis}_${mydate}.sh
echo "#SBATCH --mem ${memory}" >> ${basedir}/${analysis}_${mydate}.sh
echo "#SBATCH -o slurm.%N.%j.out" >> ${basedir}/${analysis}_${mydate}.sh
echo "#SBATCH -e slurm.%N.%j.err" >> ${basedir}/${analysis}_${mydate}.sh
echo "" >> ${basedir}/${analysis}_${mydate}.sh
echo "module load mafft/7.515" >> ${basedir}/${analysis}_${mydate}.sh
#echo "module load gblocks/0.91b" >> ${basedir}/${analysis}_${mydate}.sh
echo "module load trimal/1.4.1" >> ${basedir}/${analysis}_${mydate}.sh
echo "module load iqtree/2.3.4" >> ${basedir}/${analysis}_${mydate}.sh
echo "module load fasttree/2.1.10" >> ${basedir}/${analysis}_${mydate}.sh


echo "mafft ${mafft_options} --thread ${threads} ${database} > ${out}.faa" >> ${basedir}/${analysis}_${mydate}.sh
echo "mv ${database}.tree ${out}.tree" >> ${basedir}/${analysis}_${mydate}.sh
#echo "Gblocks ${gblocks_options} ${out}.faa" >> ${basedir}/${analysis}_${mydate}.sh
echo "trimal -in ${out}.faa ${trimal_options} -fasta -out ${out}_trimAl.faa" >> ${basedir}/${analysis}_${mydate}.sh
echo "FastTree ${fasttree_options} < ${out}_trimAl.faa > ${out}_fasttree.tree" >> ${basedir}/${analysis}_${mydate}.sh
#echo "iqtree -s ${out}_trimAl.faa ${iqtree_options} -T ${threads}" >> ${basedir}/${analysis}_${mydate}.sh

chmod 774 ${basedir}/${analysis}_${mydate}.sh

sbatch -A acinetobacterpredict --mail-type ALL --mail-user nicolas.nalpas@univ-rouen.fr -J ${analysis} ${basedir}/${analysis}_${mydate}.sh
#srun -A acinetobacterpredict --mail-type ALL --mail-user nicolas.nalpas@univ-rouen.fr -J ${analysis} ${basedir}/${analysis}_${mydate}.sh


