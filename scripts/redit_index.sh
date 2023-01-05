#!/bin/bash
PIN=pyout/archive-20230104/
echo $0
echo start
for f in $(find $PIN -name "*.txt"); do
    bgzip $f \
        && tabix -s 1 -b 2 -e 2 $f.gz &
    echo PID: $! $f

done
wait
echo done