#include "ProjBackupHeuristic.hpp"
#include "src/hyperGraph/HyperGraphExtractorDualProj.hpp"

namespace d4 {
ProjBackupHeuristic::ProjBackupHeuristic(po::variables_map &vm, SpecManager &om,
                                         WrapperSolver &s, std::ostream &out)
    : m_om(dynamic_cast<SpecManagerCnf &>(om)), m_s(s) {
  m_nbVar = m_om.getNbVariable();
  m_nbClause = m_om.getNbClause();

  m_em.initEquivExtractor(m_nbVar + 1);
  unsigned sumSize = m_om.getSumSizeClauses();

  m_partition.resize(m_nbClause + 1, 0);

  m_pm = PartitionerManager::makePartitioner(vm, m_nbClause, m_nbVar, sumSize,
                                             out);

  m_hypergraph.init(m_nbVar + m_nbClause + sumSize + 1, m_nbClause + 1);

  // initialize the vectors.
  m_markedVar.resize(m_nbVar + 1, false);
  m_equivClass.resize(m_nbVar + 1, 0);

  // get the options.
  m_reduceFormula =
      vm["partitioning-heuristic-simplification-hyperedge"].as<bool>();
  m_equivSimp = vm["proj-use-equiv"].as<bool>();
  m_hypergraphExtractor = new HyperGraphExtractorDualProj(
      m_nbVar, m_nbClause,
      vm["partitioning-heuristic-partitioner-np-cost"].as<int>());

} // constructor
ProjBackupHeuristic::~ProjBackupHeuristic() {
  if (m_hypergraphExtractor) {
    delete m_hypergraphExtractor;
  }
  if(m_pm){
      delete m_pm;
  }
}

void ProjBackupHeuristic::computeEquivClass(
    std::vector<Var> &component, std::vector<Lit> &unitEquiv,
    std::vector<Var> &equivClass, std::vector<std::vector<Var>> &equivVar) {
  if (m_equivSimp) {
    for (auto &v : component) {
      assert(equivClass.size() >= (unsigned)v);
      equivClass[v] = v;
    }

    m_em.searchEquiv(m_s, component, equivVar);
    m_s.whichAreUnits(component, unitEquiv);

    for (auto &c : equivVar) {
      Var vi = c.back();
      for (auto &v : c)
        equivClass[v] = vi;
    }

  } else {
    for (auto &v : component)
      equivClass[v] = v;
  }
} // computeEquivclass

bool ProjBackupHeuristic::computeCutSetDyn(ProjVars &component,
                                           std::vector<Var> &cutSet) {
  m_counter++;

  // search for equiv class if requiered.
  std::vector<Lit> unitEquiv;
  std::vector<std::vector<Var>> equivVar;

  computeEquivClass(component.vars, unitEquiv, m_equivClass, equivVar);

  // synchronize the SAT solver and the spec manager.
  m_om.preUpdate(unitEquiv);

  // construct the hypergraph
  std::vector<Var> considered;
  m_hypergraphExtractor->constructHyperGraph(m_om, component.vars, m_equivClass,
                                             equivVar, m_reduceFormula,
                                             considered, m_hypergraph);

  if (m_hypergraph.getSize() < 5)
    cutSet = component.vars;
  else {
    // set the level.
    PartitionerManager::Level level = PartitionerManager::Level::NORMAL;
    if (m_hypergraph.getSize() >= 200)
      level = PartitionerManager::Level::QUALITY;

    m_pm->computePartition(m_hypergraph, level, m_partition);
    m_hypergraphExtractor->extractCutFromHyperGraph(m_hypergraph, considered,
                                                    m_partition, cutSet);

    // extend with equivalence literals.
    for (auto &v : cutSet)
      m_markedVar[v] = true;
    for (auto &v : component.vars) {
      if (m_markedVar[v])
        continue;
      if (m_markedVar[m_equivClass[v]])
        cutSet.push_back(v);
    }
    for (auto &v : cutSet)
      m_markedVar[v] = false;
    if (!cutSet.size())
      for (auto l : unitEquiv)
        cutSet.push_back(l.var());
  }

  m_om.postUpdate(unitEquiv);
  return true;
}
} // namespace d4
