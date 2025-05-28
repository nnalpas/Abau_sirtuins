# Run with SLURM on IFB core cluster

# Script to run MAFFT which will perform multiple sequence alignment

analysis="mafft_sirtuins_all"
threads=8
memory=64GB
partition="fast"
database="/shared/home/nnalpas/BACTpredict/Sirtuin_conservation/All_consurf_sequences_downloaded_2024-06-12_14.fasta"
mafft_options="--amino --retree 2 --maxiterate 1000 --treeout"

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

out="`echo ${database} | sed -e 's/.fasta/_msa/'`"
echo "mafft ${mafft_options} --thread ${threads} ${database} > ${out}.faa" >> ${basedir}/${analysis}_${mydate}.sh
echo "mv ${database}.tree ${out}.tree" >> ${basedir}/${analysis}_${mydate}.sh

chmod 774 ${basedir}/${analysis}_${mydate}.sh

sbatch -A acinetobacterpredict --mail-type ALL --mail-user nicolas.nalpas@univ-rouen.fr -J ${analysis} ${basedir}/${analysis}_${mydate}.sh
#srun -A acinetobacterpredict --mail-type ALL --mail-user nicolas.nalpas@univ-rouen.fr -J ${analysis} ${basedir}/${analysis}_${mydate}.sh


