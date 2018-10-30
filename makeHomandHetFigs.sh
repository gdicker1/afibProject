touch $PWD/homFigMakeLog.txt
touch $PWD/hetFigMakeLog.txt
sleep 1
python multiFigMake.py $PWD/slurmFiles/homogeneous/ $PWD/analysis $PWD/Tissues/ $PWD/Batchtool/VepCore/ 0 29 --gif --bndTest > $PWD/homFigMakeLog.txt
python multiFigMake.py $PWD/slurmFiles/heterogeneous/ $PWD/analysis $PWD/Tissues/ $PWD/Batchtool/VepCore/ 0 29 --het --gif --bndTest > $PWD/hetFigMakeLog.txt

