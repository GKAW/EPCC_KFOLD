#!/bin/bash

getNumBasePairs () {

NAME=$1

if [ "$NAME" = "20mer" ]; then
        echo "20"
elif [ "$NAME" = "spliced" ]; then
        echo "56"
elif [ "$NAME" = "100test" ]; then
        echo "100"
elif [ "$NAME" = "stmv" ]; then
        echo "1058"
fi

}

module load valgrind

#set up array of times to run 
TIMES=( 1000 2500 3000 4500 7000 8500 10000 )

#set up array of base-pairs to run 
EXAMPLES=( "20mer" "spliced" "100test" "stmv" )

FILE="tmp"
NRUNS=5
EXAMPLE_DIR="../../examples/"

for SAMPLE in "${EXAMPLES[@]}"; do
	for TIME in "${TIMES[@]}"; do

                BASEPAIRS=$( getNumBasePairs $SAMPLE )

                # collect memory bandwidth 
                MEMORY=$( valgrind ../src/Kfold.x -i $EXAMPLE_DIR$SAMPLE".start" -o $SAMPLE".traj" -t $TIME -s 81231 2>&1 | grep "total heap usage" | awk '{ print $(NF-2) }' | sed s/,//g ) 

                # iterate over NRUNS times
                for ((i=1;i<=NRUNS;i++)); do
                        # collect time value
                        RUNTIME=$( { time ../src/Kfold.x -i $EXAMPLE_DIR$SAMPLE".start" -o $SAMPLE".traj" -t $TIME -s 81231; } 2>&1 | grep user | awk '{ print $NF}' | awk -F'm|s' '{ print 60*$1 "+" $2 }' | bc )

                        # append to file
                        echo $RUNTIME >> $FILE

                done

                MEDIAN=$( sort -n $FILE | awk -f median.awk )
                rm $FILE

                echo $SAMPLE","$BASEPAIRS","$TIME","$MEMORY","$MEDIAN

        done
done

