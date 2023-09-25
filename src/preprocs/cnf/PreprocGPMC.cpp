
#include "PreprocGPMC.hpp"
#include "gpmc.hpp"
#include "src/heuristics/cnf/PartitioningHeuristicStaticSingleProjDual.hpp"
#include "src/hyperGraph/HyperGraphExtractorDualProj.hpp"
#include "src/partitioner/PartitionerManager.hpp"
#include "src/partitioner/PartitionerPatoh.hpp"
#include "src/utils/EquivExtractor.hpp"

namespace d4 {

PreprocGPMC::PreprocGPMC(po::variables_map &vm, std::ostream &out)
    : m_vm(vm), backbone(vm, out) {
  ws = WrapperSolver::makeWrapperSolverPreproc(vm, out);
  equiv = vm["preproc-equiv"].as<bool>();
  ve_check = vm["preproc-ve-check"].as<bool>();
}

PreprocGPMC::~PreprocGPMC() { delete ws; }
ProblemManager *PreprocGPMC::run(ProblemManager *pin,
                                 LastBreathPreproc &lastBreath) {

  // For testing stuff..
  if (pin->getSelectedVar().empty()) {
    for (int i = 1; i <= pin->getNbVar(); i++) {
      pin->getSelectedVar().push_back(i);
    }
  }

  ProblemManagerCnf *cnf = (ProblemManagerCnf *)pin;

  int vars = cnf->getNbVar();
  int pvars = cnf->getSelectedVar().size();
  int freevars;
  std::vector<Lit> asignes;
  std::vector<Lit> gmap;
  std::vector<std::vector<Lit>> clauses = cnf->getClauses();
  auto make_problem = [&]() {
    std::vector<double> weight(vars + freevars + 2, 2.0);
    std::vector<double> weightLit((vars + freevars + 2) << 2, 1.0);
    for (unsigned i = 0; i <= vars; i++)
      weight[i] = weightLit[i << 1] + weightLit[(i << 1) + 1];
    std::vector<Var> selected;
    for (Var i = 1; i <= pvars; i++) {
      selected.push_back(i);
    }
    ProblemManagerCnf *out =
        new ProblemManagerCnf(vars, weightLit, weight, selected, freevars);
    for(auto cl:clauses){
        for(auto v:cl){

            assert(v.var()>0);
            
        }

    }
    out->getClauses() = std::move(clauses);
    return out;
  };
  if (!GPMC::simplify(clauses, pin->getSelectedVar(), lastBreath.learnt, vars,
                      pvars, freevars, asignes, gmap, equiv,ve_check)) {
    return pin->getUnsatProblem();
  }

  ProblemManagerCnf *out = make_problem();
  /*std::vector<int> bad_vars = find_bad_proj_vars(*out);
  for (int i = 0; i < bad_vars.size(); i++) {
    if (bad_vars[i] >= 0) {
      std::cout << "BadVar: " << i + 1 << " at " << bad_vars[i] << std::endl;
    }
  }*/
  // delete out;

  ws->initSolver(*out);
  lastBreath.panic = 0;
  lastBreath.countConflict.resize(out->getNbVar() + 1, 0);

  if (!ws->solve())
    return out->getUnsatProblem();
  lastBreath.panic = ws->getNbConflict() > 100000;

  // get the activity given by the solver.
  for (unsigned i = 1; i <= out->getNbVar(); i++)
    lastBreath.countConflict[i] = ws->getCountConflict(i);

  std::vector<Lit> units;
  ws->getUnits(units);
  auto final = out->getConditionedFormula(units);

  delete out;
  return final;
}

} // namespace d4
