cd  build
ninja
./d4 -m proj-ddnnf-compiler   --partitioning-heuristic decomposition-static-proj-dual\
    --partitioning-heuristic-partitioner kahypar   \
    --partitioning-heuristic-simplification-equivalence true \
    --partitioning-heuristic-partitioner-np-cost 100 \
    -i ~/projects/projected-ddnnf-compilation-eval/instances/fm-hard/automotive02.automotive2_4/28.dimacs\
    -p gpmc \
    --crs none \
    --cache-method lru\
    --proj-backup none\
    --preproc-equiv true\
    --preproc-ve-check true\



 
    #-i ~/projects/projected-ddnnf-compilation-eval/instances/fm-gen/KConfig.linux-2.6.33.3/16.dimacs\



