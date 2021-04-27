#!/bin/bash

find ../src_omp/ -name "*f90" | xargs sed -i '/DO.*=.*,.*/i!\$acc kernels loop'
DELLINE=$(grep -n "DO i=1,nn" ../src_omp/kfold.f90 | awk -F: '{ print $1 - 1}' | bc) 
sed -i -e "$DELLINE"'d' ../src_omp/kfold.f90
