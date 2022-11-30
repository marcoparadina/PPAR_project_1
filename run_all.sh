!#/bin/bash

for host in $(cat "$OAR_NODEFILE" | uniq);
do
    oarsh "$host" sudo-g5k apt install libscalapack-openmpi-dev -y &
done;

make

mpirun -n 961 -machinefile "$OAR_NODEFILE" ./pmodel3_clean --data /srv/storage/ppar@storage1.nancy.grid5000.fr/ETOPO1_ultra.csv --model out.csv --npoint 120000 --lmax 124