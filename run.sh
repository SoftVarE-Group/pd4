cd  build
ninja
./d4 -m proj-ddnnf-compiler   --partitioning-heuristic decomposition-static-proj-dual \
    --partitioning-heuristic-partitioner mtkahypar  --partitioning-heuristic-simplification-equivalence true \
    --partitioning-heuristic-partitioner-np-cost 100  \
    -i  ~/projects/projected-ddnnf-compilation-eval/instances/bb/Smarch.embtoolkit/0.dimacs\
    -p proj \
    --cache-impl conf \
    --proj-backup none\
    --proj-backup-min-cover 0.4 \
    --proj-backup-scoring-method 0 \
    --occurrence-manager dynamic \
    --scoring-method vsads \
    --scoring-method-x 1 \
    --scoring-method-y 128 \
    --scoring-method-z 128 \

