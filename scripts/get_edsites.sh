#!/bin/bash
echo $0
echo start
PIN=pyout/archive-20230104/
POUT=shout/edsites/

for f in $(find $PIN -name "*.txt"); do
    awk \
        -v FS='\t' \
        -v OFS='\t' \
        '$5 > 10 && $9 > 0.1 && $14 == 0.00 {print $1,$2}' \
        $f \
        > $POUT/$(basename $f).tmp &
    echo PID: $! $f
done
wait


cat $POUT/ptr-*.tmp \
    | sort -k1,1 -k2n,2 \
    | uniq \
    > $POUT/ptr.txt &

cat $POUT/egr-*.tmp \
    | sort -k1,1 -k2n,2 \
    | uniq \
    > $POUT/egr.txt &

wait
echo finished