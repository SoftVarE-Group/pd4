#pragma once
#include <vector>

#include "src/preprocs/PreprocManager.hpp"
#include "3rdParty/glucose-3.0/core/Solver.h"
#include "src/preprocs/cnf/PreprocBackboneCnf.hpp"
#include "src/problem/cnf/ProblemManagerCnf.hpp"
#include "src/heuristics/cnf/PartitioningHeuristicStaticSingleProjDual.hpp"



namespace d4{

class PreprocGPMC:public PreprocManager{

  WrapperSolver *ws;
  po::variables_map& m_vm;
  PreprocBackboneCnf backbone;
  bool equiv = true;
  bool ve_check = true;


public: 
  PreprocGPMC(po::variables_map &vm, std::ostream &out);
  virtual ~PreprocGPMC();
  
  virtual ProblemManager *run(ProblemManager *pin,
                              LastBreathPreproc &lastBreath) override;

};
}




namespace gpmc {





}
