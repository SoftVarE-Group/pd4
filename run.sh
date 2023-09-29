cd  build
ninja
./d4 -m proj-ddnnf-compiler   --partitioning-heuristic decomposition-static-proj-dual\
    --partitioning-heuristic-partitioner kahypar   \
    --partitioning-heuristic-simplification-equivalence true \
    --partitioning-heuristic-partitioner-np-cost 100 \
    -i ~/projects/projected-ddnnf-compilation-eval/instances/fm-hard/automotive01.automotive01/17.dimacs\
    -p proj \
    --proj-backup none \
    --crs none \
    --cache-alloc std\
    --cache-method lru-prob\
    --preproc-equiv true \
    --preproc-ve-check true\
    --preproc-ve-only-simpical false\
    --preproc-ve-prefer-simpical true\
    --preproc-ve-limit 4



 
    #-i ~/projects/projected-ddnnf-compilation-eval/instances/fm-gen/KConfig.linux-2.6.33.3/16.dimacs\



