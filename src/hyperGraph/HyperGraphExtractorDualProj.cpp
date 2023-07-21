
#include "HyperGraphExtractorDualProj.hpp"
namespace d4 {
HyperGraphExtractorDualProj::HyperGraphExtractorDualProj(unsigned nbVar,
                                                         unsigned nbClause,
                                                         int nProjCost)
    : HyperGraphExtractorDual(nbVar, nbClause), m_nProjCost(nProjCost) {}
void HyperGraphExtractorDualProj::constructHyperGraph(
    SpecManagerCnf &om, std::vector<Var> &component,
    std::vector<Var> &equivClass, std::vector<std::vector<Var>> &equivVar,
    bool reduceFormula, std::vector<Var> &considered, HyperGraph &hypergraph) {
  unsigned pos = 0;
  m_idxClauses.resize(0);
  hypergraph.setSize(0);
  // std::cout << "Cost " << m_nProjCost << std::endl;

  // first considere the equivalence.

  for (auto &vec : equivVar) {
    unsigned &size = hypergraph[pos++];
    size = 0;
    bool has_proj = false;

    for (auto &v : vec) {
      if (om.varIsAssigned(v))
        continue;

      has_proj |= om.isSelected(v);
      assert(!m_markedVar[v]);
      for (auto l : {Lit::makeLitFalse(v), Lit::makeLitTrue(v)}) {
        IteratorIdxClause listIdx = om.getVecIdxClauseNotBin(l);
        for (int *ptr = listIdx.start; ptr != listIdx.end; ptr++) {
          int idx = *ptr;
          if (!m_markedClauses[idx]) {
            m_markedClauses[idx] = true;
            m_unmarkSet.push_back(idx);
            hypergraph[pos++] = idx;
            size++;
          }
        }

        listIdx = om.getVecIdxClauseBin(l);
        for (int *ptr = listIdx.start; ptr != listIdx.end; ptr++) {
          int idx = *ptr;
          if (!m_markedClauses[idx]) {
            m_markedClauses[idx] = true;
            m_unmarkSet.push_back(idx);
            hypergraph[pos++] = idx;
            size++;
          }
        }
      }

      m_markedVar[v] = true;
    }

    for (auto &idx : m_unmarkSet) {
      if (!m_keepClause[idx]) {
        m_keepClause[idx] = true;
        m_idxClauses.push_back(idx);
      }
      m_markedClauses[idx] = false;
    }
    m_unmarkSet.resize(0);
    assert(equivClass[vec.back()] == vec.back());

    if (!size)
      pos--;
    else {

      hypergraph.getCost()[hypergraph.getSize()] = has_proj ? 1 : m_nProjCost;
      hypergraph.incSize();
      considered.push_back(vec.back());
    }
  }

  // next consider the remaining variables (unmarked).
  for (auto &v : component) {
    if (m_markedVar[v] || om.varIsAssigned(v))
      continue;
    m_markedVar[v] = true;

    unsigned &size = hypergraph[pos++];
    size = 0;

    for (auto l : {Lit::makeLitFalse(v), Lit::makeLitTrue(v)}) {
      IteratorIdxClause listIdx = om.getVecIdxClauseNotBin(l);
      for (int *ptr = listIdx.start; ptr != listIdx.end; ptr++) {
        int idx = *ptr;
        if (!m_keepClause[idx]) {
          m_keepClause[idx] = true;
          m_idxClauses.push_back(idx);
        }
        hypergraph[pos++] = idx;
        size++;
      }

      listIdx = om.getVecIdxClauseBin(l);
      for (int *ptr = listIdx.start; ptr != listIdx.end; ptr++) {
        int idx = *ptr;
        if (!m_keepClause[idx]) {
          m_keepClause[idx] = true;
          m_idxClauses.push_back(idx);
        }
        hypergraph[pos++] = idx;
        size++;
      }
    }

    if (!size)
      pos--;
    else {
      bool sel = om.isSelected(v);
      hypergraph.getCost()[hypergraph.getSize()] = sel ? 1 : m_nProjCost;
      hypergraph.incSize();
      considered.push_back(v);
    }
  }

  // unmark.
  for (auto &v : component)
    m_markedVar[v] = false;

  // remove useless edges.
  if (reduceFormula)
    reduceHyperGraph(om, hypergraph, considered, m_idxClauses, equivClass);

  // unmark.
  for (auto &idx : m_idxClauses)
    m_keepClause[idx] = false;
}
} // namespace d4
