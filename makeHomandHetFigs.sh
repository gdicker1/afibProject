touch /extHDD/slurmFiles/homFigMakeLog.txt
touch /extHDD/slurmFiles/hetFigMakeLog.txt
sleep 2
python multiFigMake.py /extHDD/slurmFiles/homogeneous/ /extHDD/slurmFiles/analysis /gscratch/gdicker1/CMAES_runner/Tissues/ /gscratch/gdicker1/CMAES_runner/Batchtool/VepCore/ 0 29 --gif > /extHDD/slurmFiles/homFigMakeLog.txt
python multiFigMake.py /extHDD/slurmFiles/heterogeneous/ /extHDD/slurmFiles/analysis /gscratch/gdicker1/CMAES_runner/Tissues/ /gscratch/gdicker1/CMAES_runner/Batchtool/VepCore/ 0 29 --het --gif > /extHDD/slurmFiles/hetFigMakeLog.txt

