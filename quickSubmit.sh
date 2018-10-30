for i in {0..29}
do
  cd /gscratch/gdicker1/CMAES_runner/slurmFiles/homogeneous/set${i}
  sbatch afib_hom${i}.pbs
  cd /gscratch/gdicker1/CMAES_runner/slurmFiles/heterogeneous/set${i}
  sbatch afib_het${i}.pbs
done
