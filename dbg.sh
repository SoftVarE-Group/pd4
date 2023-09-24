cd build 
cmake .. -DCMAKE_BUILD_TYPE=Debug -G Ninja
ninja
gdb --args ./d4 -m proj-ddnnf-compiler   --partitioning-heuristic decomposition-static-dual\
    --partitioning-heuristic-partitioner kahypar   \
    --partitioning-heuristic-simplification-equivalence true \
    --partitioning-heuristic-partitioner-np-cost 100 \
    -i ~/projects/projected-ddnnf-compilation-eval/instances/fm-hard/automotive01.automotive01/95.dimacs\
    -p gpmc \
    --crs none \
    --cache-method lru\
    --proj-backup none\
    --preproc-equiv true\
    --preproc-ve-check true\
