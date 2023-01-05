#!/bin/bash
DIR_IN=$1
DIR_OUT=$2

# echo arguments
echo DIR_IN, DIR_OUT = $DIR_IN, $DIR_OUT

# 
for f in $(ls $DIR_IN); do
    # awk -F "\t" -v OFS="\t" '$5 > 10 && $9 > 0.1 && $14 == 0.00 {print}' ${DIR_IN}/${f} | head
    # awk -F "\t" -v OFS="\t" '$5 > 10 && $9 > 0.1 && $14 == 0.00 {print}' ${DIR_IN}/${f} > ${DIR_OUT}/${f} &
    awk \
        -v FS='\t' \
        -v OFS='\t' \
        '$5 > 10 && $14 == 0.00 {print}' ${DIR_IN}/${f} \
        > ${DIR_OUT}/${f} &
    echo PID:$!, $f
done
wait

echo Done