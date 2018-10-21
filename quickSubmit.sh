for i in {0..29}
do
  cd /gscratch/afibProject/slurmFiles/homogeneous/set${i}
  sbatch afib_hom${i}.pbs
  cd /gscratch/afibProject/slurmFiles/heterogeneous/set${i}
  sbatch afib_het${i}.pbs
done
