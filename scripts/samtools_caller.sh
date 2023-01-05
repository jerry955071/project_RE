DIR_IN=processed/sorted_bam/
DIR_OUT=shout/samtools_out/
MAX_CONTAINER=16

# for f in $(find $DIR_IN -name "*r.sorted.bam")
# do
#     # calculate idxstats
#     # docker run \
#     #     --rm \
#     #     --cpus .5 \
#     #     -v $(pwd):/local:rw \
#     #     -w /local \
#     #     --name $(basename $f)_idxstats \
#     #     ccc/reditools2:latest \
#     #         samtools idxstats $f > $DIR_OUT/idxstats/$(basename $f).idxstats &
#     # echo $f: $!

#     # calcualte coverage
#     # docker run \
#     #     --rm \
#     #     --cpus .5 \
#     #     -v $(pwd):/local:rw \
#     #     -w /local \
#     #     --name $(basename $f)_coverage \
#     #     ccc/reditools2:latest \
#     #         samtools coverage $f > $DIR_OUT/coverage/$(basename $f).coverage &
#     # echo $f: $!

#     # calculate depth
#     echo $f: samtools
#     docker run \
#         --rm \
#         --cpus 1 \
#         -v $(pwd):/local:rw \
#         -w /local \
#         --name $(basename $f)_depth \
#         ccc/reditools2:latest \
#             samtools depth $f > $DIR_OUT/depth/$(basename $f).tmp &
#     sleep 5
#     while [[ $(docker ps | grep samtools | wc -l) -ge $MAX_CONTAINER ]]; do
#         sleep 300
#     done

# done
# wait

echo $0 $!
for f in $(find $DIR_IN -name "*r.sorted.bam")
do
    echo $f
    samtools depth $f | \
        awk -v OFS="\t" '$3 > 0 {print}' > \
        $DIR_OUT/depth/$(basename $f).depth &
done
wait