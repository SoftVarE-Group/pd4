#include "PartitionerPatoh.hpp"
namespace d4 {
class PartitionerPatohProj : public PartitionerPatoh {
private:
  int *m_cost;

public:
  PartitionerPatohProj(unsigned maxNodes, unsigned maxEdges,
                       unsigned maxSumEdgeSize, std::ostream &out);
  ~PartitionerPatohProj();
  void computePartition(HyperGraph &hypergraph, Level level,
                        std::vector<int> &partition);
};
} // namespace d4
