# Run with SLURM on IFB core cluster

# Script to run FastTree to combine multiple sequence alignment with NCBI taxonomy tree constraint

analysis="fasttree_ncbi"
threads=7
memory=64GB
partition="fast"
msa="/shared/home/nnalpas/BACTpredict/NCBI_taxonomy/Alignment_concatenation_FIXED.fasta"
tree="/shared/home/nnalpas/BACTpredict/NCBI_taxonomy/NCBI_taxonomy_superkingdom.tree"
options="-slow"

basedir=`dirname $msa`
mydate=`date +%Y-%m-%d_%H-%M-%S`

echo '#!/bin/bash' > ${basedir}/${analysis}_${mydate}.sh
echo '#' >> ${basedir}/${analysis}_${mydate}.sh
echo "#SBATCH -p ${partition}" >> ${basedir}/${analysis}_${mydate}.sh
#echo "#SBATCH --cpus-per-task=${threads}" >> ${basedir}/${analysis}_${mydate}.sh
echo "#SBATCH --ntasks=${threads}" >> ${basedir}/${analysis}_${mydate}.sh
echo "#SBATCH --mem-per-cpu ${memory}" >> ${basedir}/${analysis}_${mydate}.sh
echo "#SBATCH -o slurm.%N.%j.out" >> ${basedir}/${analysis}_${mydate}.sh
echo "#SBATCH -e slurm.%N.%j.err" >> ${basedir}/${analysis}_${mydate}.sh
echo "" >> ${basedir}/${analysis}_${mydate}.sh
echo "module load perl/5.26.2" >> ${basedir}/${analysis}_${mydate}.sh
echo "module load fasttree/2.1.10" >> ${basedir}/${analysis}_${mydate}.sh

constraint="`echo ${tree} | sed -e 's/.tree/_FT.constraint/'`"
out="`echo ${msa} | sed -e 's/.fasta/_FTslow.tree/'`"
echo "perl ${HOME}/Software/TreeToConstraints/TreeToConstraints.pl < ${tree} > ${constraint}" >> ${basedir}/${analysis}_${mydate}.sh
echo "FastTree ${options} -constraints ${constraint} -out ${out} ${msa}" >> ${basedir}/${analysis}_${mydate}.sh

chmod 774 ${basedir}/${analysis}_${mydate}.sh

sbatch -A acinetobacterpredict --mail-type ALL --mail-user nicolas.nalpas@univ-rouen.fr -J ${analysis} ${basedir}/${analysis}_${mydate}.sh
#srun -A acinetobacterpredict --mail-type ALL --mail-user nicolas.nalpas@univ-rouen.fr -J ${analysis} ${basedir}/${analysis}_${mydate}.sh


