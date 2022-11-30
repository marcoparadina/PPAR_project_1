!#/bin/bash

for host in $(cat "$OAR_NODEFILE" | uniq);
do
    oarsh "$host" sudo-g5k apt install libscalapack-openmpi-dev -y &
done;

