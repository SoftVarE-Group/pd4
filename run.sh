cd  build
ninja
./d4 -m proj-counting   --partitioning-heuristic decomposition-static-proj-dual\
    --partitioning-heuristic-partitioner kahypar   \
    --partitioning-heuristic-simplification-equivalence true \
    --partitioning-heuristic-partitioner-np-cost 100 \
    --partitioning-heuristic-max-cut-ratio 0.5 \
    -i ~/Documents/MC2022_track3-pmc_private/mc2022_track3_117/mc2022_track3_117.dimacs\
    --crs none \
    -p proj\
    --sm vsads2\
    --cache-alloc std\
    --cache-method lru-prob\
    --cache-fixed-size 8 \
    --preproc-equiv true \
    --preproc-ve-check false\
    --preproc-ve-only-simpical true\
    --preproc-ve-prefer-simpical false\
    --preproc-ve-limit 4



 



