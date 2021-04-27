#
#        FILE: reconstruct.sh
# DESCRIPTION: To enable parallel running, each simulation writes its own .traj and .log files, with
#              the simulation number appended, i.e. foo.log.1, foo.log.2, etc. This script
#              concatenates the files together to generate an equivalent file to if it had been run
#              sequentially.
#       USAGE: /path/to/reconstruct.sh <file root> <number of simulations>
#              For example, if you had generated a set of 4 experiments called foo.log.x and
#              foo.traj.x (1 <= x <= 4) then you would call
#                /path/to/reconstruct.sh foo 4
#              to reconstruct the file.
#

# Get command line options
f=${1}
n=${2}

# Perform reconstruction
for i in $(seq 1 ${n})
do
    cat ${f}.log.${i} >> ${f}.log.par
    cat ${f}.traj.${i} >> ${f}.traj.par
done
