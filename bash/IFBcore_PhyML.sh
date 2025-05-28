# Run with SLURM on IFB core cluster

# Script to run phyML to combine multiple sequence alignment with NCBI taxonomy tree profile

analysis="phyml_test"
threads=4
memory=32GB
partition="fast"
msa="/shared/home/nnalpas/BACTpredict/NCBI_taxonomy/concatenation.phy"
tree="/shared/home/nnalpas/BACTpredict/NCBI_taxonomy/phyliptree_unranked.nwk"
options=" -d aa -b 4 -m LG -c 4 -a e -o lr"

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
echo "module load phyml/3.3.20190909" >> ${basedir}/${analysis}_${mydate}.sh

echo "mpirun phyml-mpi -i ${msa} -u ${tree} ${options}" >> ${basedir}/${analysis}_${mydate}.sh

chmod 774 ${basedir}/${analysis}_${mydate}.sh

sbatch -A acinetobacterpredict --mail-type ALL --mail-user nicolas.nalpas@univ-rouen.fr -J ${analysis} ${basedir}/${analysis}_${mydate}.sh
#srun -A acinetobacterpredict --mail-type ALL --mail-user nicolas.nalpas@univ-rouen.fr -J ${analysis} ${basedir}/${analysis}_${mydate}.sh


