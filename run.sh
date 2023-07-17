sh build.sh
cd build
./d4 -i ../0.dimacs  -m proj-ddnnf-compiler   --partitioning-heuristic bipartition-dual-proj --partitioning-heuristic-partitioner kahypar  --partitioning-heuristic-simplification-equivalence false --partitioning-heuristic-partitioner-np-cost 101 --partitioning-heuristic-bipartite-phase-dynamic 0.1 --partitioning-heuristic-bipartite-phase proj
 

./d4 -i ../0.dimacs  -m proj-ddnnf-compiler   --partitioning-heuristic decomposition-static-proj-dual --partitioning-heuristic-partitioner kahypar  --partitioning-heuristic-simplification-equivalence false --partitioning-heuristic-partitioner-np-cost 101 
