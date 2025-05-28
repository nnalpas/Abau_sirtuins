# Run with SLURM on IFB core cluster

# Script to run RAxML to combine multiple sequence alignment with NCBI taxonomy tree profile

analysis="raxml_ncbi"
threads=28
memory=128GB
partition="long"
msa="/shared/home/nnalpas/BACTpredict/NCBI_taxonomy/Alignment_concatenation_FIXED.phy"
tree="/shared/home/nnalpas/BACTpredict/NCBI_taxonomy/NCBI_taxonomy_superkingdom.tree"
options="--msa-format PHYLIP --data-type AA --model LG+G8+F --tree pars{4} --bs-trees 4 -all"

basedir=`dirname $msa`
mydate=`date +%Y-%m-%d_%H-%M-%S`

echo '#!/bin/bash' > ${basedir}/${analysis}_${mydate}.sh
echo '#' >> ${basedir}/${analysis}_${mydate}.sh
echo "#SBATCH -p ${partition}" >> ${basedir}/${analysis}_${mydate}.sh
echo "#SBATCH --cpus-per-task=${threads}" >> ${basedir}/${analysis}_${mydate}.sh
echo "#SBATCH --mem ${memory}" >> ${basedir}/${analysis}_${mydate}.sh
echo "#SBATCH -o slurm.%N.%j.out" >> ${basedir}/${analysis}_${mydate}.sh
echo "#SBATCH -e slurm.%N.%j.err" >> ${basedir}/${analysis}_${mydate}.sh
echo "" >> ${basedir}/${analysis}_${mydate}.sh
echo "module load raxml-ng/1.1.0" >> ${basedir}/${analysis}_${mydate}.sh

echo "raxml-ng --msa ${msa} --tree-constraint ${tree} ${options} --threads ${threads}" >> ${basedir}/${analysis}_${mydate}.sh

chmod 774 ${basedir}/${analysis}_${mydate}.sh

sbatch -A acinetobacterpredict --mail-type ALL --mail-user nicolas.nalpas@univ-rouen.fr -J ${analysis} ${basedir}/${analysis}_${mydate}.sh
#srun -A acinetobacterpredict --mail-type ALL --mail-user nicolas.nalpas@univ-rouen.fr -J ${analysis} ${basedir}/${analysis}_${mydate}.sh


