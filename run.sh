cd  build
ninja
./d4 -m proj-counting  --partitioning-heuristic decomposition-static-proj-dual\
    --partitioning-heuristic-partitioner kahypar   \
    --partitioning-heuristic-simplification-equivalence true \
    --partitioning-heuristic-partitioner-np-cost 100 \
    --partitioning-heuristic-max-cut-ratio 0.6 \
    -i ~/projects/projected-ddnnf-compilation-eval/instances/fm-real/cars/0.dimacs \
    --crs none \
    -p proj\
    --sm vsads2\
    --cache-alloc std\
    --cache-method lru-prob\
    --cache-fixed-size 8 \
    --preproc-equiv true \
    --preproc-ve-check true\
    --preproc-ve-only-simpical false\
    --preproc-ve-prefer-simpical true\
    --preproc-ve-limit 4 \
    --projddnnf-pure-lit-elim true \
    --scoring-method-decay-freq 128




 



