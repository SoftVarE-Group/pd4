
/*
 * d4
 * Copyright (C) 2020  Univ. Artois & CNRS
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
#pragma once

#include "src/heuristics/cnf/PartitioningHeuristicBipartiteDual.hpp"
#include <algorithm>
#include <boost/program_options.hpp>
#include <ctime>
#include <iomanip>
#include <iostream>
#include <vector>

#include "Counter.hpp"
#include "DataBranch.hpp"
#include "MethodManager.hpp"
#include "src/caching/CacheManager.hpp"
#include "src/caching/CachedBucket.hpp"
#include "src/caching/TmpEntry.hpp"
#include "src/heuristics/PartitioningHeuristic.hpp"
#include "src/heuristics/PhaseHeuristic.hpp"
#include "src/heuristics/ProjBackupHeuristic.hpp"
#include "src/heuristics/ScoringMethod.hpp"
#include "src/preprocs/PreprocManager.hpp"
#include "src/problem/ProblemManager.hpp"
#include "src/problem/ProblemTypes.hpp"
#include "src/solvers/WrapperSolver.hpp"
#include "src/specs/SpecManager.hpp"
#include "src/utils/MemoryStat.hpp"
#include "src/utils/Proj.hpp"

#define NB_SEP_MC 104
#define MASK_SHOWRUN_MC ((2 << 13) - 1)
#define WIDTH_PRINT_COLUMN_MC 12
#define MASK_HEADER 1048575

#include "CountingOperation.hpp"
#include "DecisionDNNFOperation.hpp"
#include "OperationManager.hpp"
#include <chrono>
template <class result_t = std::chrono::milliseconds,
          class clock_t = std::chrono::steady_clock,
          class duration_t = std::chrono::milliseconds>
auto since(std::chrono::time_point<clock_t, duration_t> const &start) {
  return std::chrono::duration_cast<result_t>(clock_t::now() - start);
}

namespace d4 {
namespace po = boost::program_options;
template <class T> class Counter;

template <class T, class U>
class ProjDpllStyleMethod : public MethodManager, public Counter<T> {
private:
  bool optDomConst;
  bool optReversePolarity;
  unsigned m_failed_cuts = 0;
  unsigned m_nbCallCall;
  unsigned m_nbSplit;
  unsigned m_callPartitioner;
  unsigned m_nbDecisionNode;
  unsigned m_optCached;
  unsigned m_stampIdx;
  bool m_isProjectedMode;

  std::vector<unsigned> m_stampVar;
  std::vector<std::vector<Lit>> m_clauses;

  std::vector<unsigned long> nbTestCacheVarSize;
  std::vector<unsigned long> nbPosHitCacheVarSize;

  std::vector<bool> m_currentPrioritySet;

  ProblemManager *m_problem;
  WrapperSolver *m_solver;
  SpecManager *m_specs;
  ScoringMethod *m_hVar;
  PhaseHeuristic *m_hPhase;
  PartitioningHeuristic *m_hCutSet;
  std::unique_ptr<ProjBackupHeuristic> m_hBackUp;

  TmpEntry<U> NULL_CACHE_ENTRY;
  CacheManager<U> *m_cache;
  std::ostream m_out;
  bool m_panicMode;

  Operation<T, U> *m_operation;
  std::vector<int> m_var_sign;

public:
  /**
     Constructor.

     @param[in] vm, the list of options.
   */
  ProjDpllStyleMethod(po::variables_map &vm, std::string &meth, bool isFloat,
                      ProblemManager *initProblem, std::ostream &out,
                      LastBreathPreproc &lastBreath)
      : m_problem(initProblem), m_out(nullptr) {
    // init the output stream
    m_out.copyfmt(out);
    m_out.clear(out.rdstate());
    m_out.basic_ios<char>::rdbuf(out.rdbuf());

    // we create the SAT solver.
    m_solver = WrapperSolver::makeWrapperSolver(vm, m_out);
    assert(m_solver);
    m_panicMode = lastBreath.panic;

    m_solver->initSolver(*m_problem, lastBreath.learnt);
    m_solver->setCountConflict(lastBreath.countConflict, 1,
                               m_problem->getNbVar());
    m_solver->setNeedModel(true);

    // we initialize the object that will give info about the problem.
    m_specs = SpecManager::makeSpecManager(vm, *m_problem, m_out);
    assert(m_specs);

    // we initialize the object used to compute score and partition.
    m_hVar = ScoringMethod::makeScoringMethod(vm, *m_specs, *m_solver, m_out);
    m_hPhase =
        PhaseHeuristic::makePhaseHeuristic(vm, *m_specs, *m_solver, m_out);

    m_isProjectedMode = true;

    // select the partitioner regarding if it projected model counting or not.
    /*
    if ((m_isProjectedMode = m_problem->getNbSelectedVar())) {
      m_out << "c [MODE] projected\n";
      m_hCutSet = PartitioningHeuristic::makePartitioningHeuristicNone(m_out);
    } else {
      m_out << "c [MODE] classic\n";
      m_hCutSet = PartitioningHeuristic::makePartitioningHeuristic(
          vm, *m_specs, *m_solver, m_out);
    }
    */

    m_currentPrioritySet.resize(m_problem->getNbVar() + 1, false);
    if ((double)m_specs->nbSelected() / (double)m_specs->getNbVariable() <
            0.10 &&
        (m_specs->getNbVariable() > 1000)) {
      m_hCutSet = PartitioningHeuristic::makePartitioningHeuristicNone(m_out);
    } else {
      m_hCutSet = PartitioningHeuristic::makePartitioningHeuristic(
          vm, *m_specs, *m_solver, m_out);
    }

    assert(m_hVar && m_hPhase && m_hCutSet);
    m_cache = CacheManager<U>::makeCacheManager(vm, m_problem->getNbVar(),
                                                m_specs, m_out);

    // init the clock time.
    initTimer();

    m_optCached = vm["cache-activated"].as<bool>();
    m_callPartitioner = 0;
    m_nbDecisionNode = m_nbSplit = m_nbCallCall = 0;
    m_stampIdx = 0;
    m_stampVar.resize(m_specs->getNbVariable() + 1, 0);
    nbTestCacheVarSize.resize(m_specs->getNbVariable() + 1, 0);
    nbPosHitCacheVarSize.resize(m_specs->getNbVariable() + 1, 0);

    void *op = Operation<T, U>::makeOperationManager(meth, isFloat, m_problem,
                                                     m_specs, m_solver, m_out);
    m_operation = static_cast<Operation<T, U> *>(op);
    m_out << "c\n";
    m_var_sign.resize(m_specs->getNbVariable() + 1, 0);
    m_hBackUp = ProjBackupHeuristic::make(vm, *m_specs, *m_solver, m_out);

  } // constructor

  /**
     Destructor.
   */
  ~ProjDpllStyleMethod() {
    delete m_operation;
    delete m_problem;
    delete m_solver;
    delete m_specs;
    delete m_hVar;
    delete m_hPhase;
    delete m_hCutSet;
    delete m_cache;
  } // destructor
  int m_frequency = 0;

private:
  int reduceNProjVars(ProjVars &vars) {
    if (m_frequency % 100000 != 0 )
      return 0;
    m_frequency++;
    int assume_npoj = 0;

    SpecManagerCnf *specs = ((SpecManagerCnf *)m_specs);
    // specs->computeCleaness();
    for (auto v : vars.iter_nproj()) {
      int p = specs->getNbOccurrence(Lit::makeLit(v, true));
      if (p == 0) {
        m_solver->pushAssumption(Lit::makeLit(v, false));
        assume_npoj++;
      } else {
        int n = specs->getNbOccurrence(Lit::makeLit(v, false));
        if (n == 0) {
          m_solver->pushAssumption(Lit::makeLit(v, true));
          assume_npoj++;
        }
      }
    }
    return assume_npoj;
  }
  void reduceNProjVarsPost(int assume_nproj) {
    if (assume_nproj > 0) {
      std::cout << "Nuked " << assume_nproj << std::endl;
    }
    m_solver->popAssumption(assume_nproj);
  }
  void completeCutset(std::vector<Var> &vars) {}
  /**
     Expel from a set of variables the ones they are marked as being decidable.

     @param[out] vars, the set of variables we search to filter.

     @param[in] isDecisionvariable, a boolean vector that marks as true decision
     variables.
   */
  void expelNoDecisionVar(std::vector<Var> &vars) {

    unsigned j = 0;
    for (unsigned i = 0; i < vars.size(); i++)
      if (m_specs->isSelected(vars[i]))
        vars[j++] = vars[i];
    vars.resize(j);
  } // expelNoDecisionVar

  /**
     Expel from a set of variables the ones they are marked as being decidable.

     @param[out] lits, the set of literals we search to filter.

     @param[in] isDecisionvariable, a boolean vector that marks as true decision
     variables.
   */
  void expelNoDecisionLit(std::vector<Lit> &lits) {
    if (!m_isProjectedMode)
      return;

    unsigned j = 0;
    for (unsigned i = 0; i < lits.size(); i++)
      if (m_specs->isSelected(lits[i].var()))
        lits[j++] = lits[i];
    lits.resize(j);
  } // expelNoDecisionLit

  /**
     Compute the current priority set.

     @param[in] connected, the current component.
     @param[in] priorityVar, the current priority variables.
     @param[out] currPriority, the intersection of the two previous sets.
  */
  inline void computePrioritySubSet(std::vector<Var> &connected,
                                    std::vector<Var> &priorityVar,
                                    std::vector<Var> &currPriority) {
    currPriority.resize(0);
    m_stampIdx++;
    for (auto &v : connected)
      m_stampVar[v] = m_stampIdx;
    for (auto &v : priorityVar)
      if (m_stampVar[v] == m_stampIdx && !m_specs->varIsAssigned(v))
        currPriority.push_back(v);
  } // computePrioritySet

  /**
     Print out information about the solving process.

     @param[in] out, the stream we use to print out information.
  */
  inline void showInter(std::ostream &out) {
    out << "c " << std::fixed << std::setprecision(2) << "|"
        << std::setw(WIDTH_PRINT_COLUMN_MC) << getTimer() << "|"
        << std::setw(WIDTH_PRINT_COLUMN_MC) << m_cache->getNbPositiveHit()
        << "|" << std::setw(WIDTH_PRINT_COLUMN_MC)
        << m_cache->getNbNegativeHit() << "|"
        << std::setw(WIDTH_PRINT_COLUMN_MC) << m_cache->usedMemory() << "|"
        << std::setw(WIDTH_PRINT_COLUMN_MC) << m_nbSplit << "|"
        << std::setw(WIDTH_PRINT_COLUMN_MC) << MemoryStat::memUsedPeak() << "|"
        << std::setw(WIDTH_PRINT_COLUMN_MC) << m_nbDecisionNode << "|"
        << std::setw(WIDTH_PRINT_COLUMN_MC) << m_callPartitioner
        << std::setw(WIDTH_PRINT_COLUMN_MC) << m_failed_cuts << "|\n";
  } // showInter

  /**
     Print out a line of dashes.

     @param[in] out, the stream we use to print out information.
   */
  inline void separator(std::ostream &out) {
    out << "c ";
    for (int i = 0; i < NB_SEP_MC; i++)
      out << "-";
    out << "\n";
  } // separator

  /**
     Print out the header information.

     @param[in] out, the stream we use to print out information.
  */
  inline void showHeader(std::ostream &out) {
    separator(out);
    out << "c "
        << "|" << std::setw(WIDTH_PRINT_COLUMN_MC) << "time"
        << "|" << std::setw(WIDTH_PRINT_COLUMN_MC) << "#posHit"
        << "|" << std::setw(WIDTH_PRINT_COLUMN_MC) << "#negHit"
        << "|" << std::setw(WIDTH_PRINT_COLUMN_MC) << "memory"
        << "|" << std::setw(WIDTH_PRINT_COLUMN_MC) << "#split"
        << "|" << std::setw(WIDTH_PRINT_COLUMN_MC) << "mem(MB)"
        << "|" << std::setw(WIDTH_PRINT_COLUMN_MC) << "#dec. Node"
        << "|" << std::setw(WIDTH_PRINT_COLUMN_MC) << "#cutter"
        << "|" << std::setw(WIDTH_PRINT_COLUMN_MC) << "#bad cuts"
        << "|\n";
    separator(out);
  } // showHeader

  /**
     Print out information when it is requiered.

     @param[in] out, the stream we use to print out information.
   */
  inline void showRun(std::ostream &out) {
    if (!(m_nbCallCall & (MASK_HEADER)))
      showHeader(out);
    if (m_nbCallCall && !(m_nbCallCall & MASK_SHOWRUN_MC))
      showInter(out);
  } // showRun

  /**
     Print out the final stat.

     @param[in] out, the stream we use to print out information.
   */
  inline void printFinalStats(std::ostream &out) {
    separator(out);
    out << "c\n";
    out << "c \033[1m\033[31mStatistics \033[0m\n";
    out << "c \033[33mCompilation Information\033[0m\n";
    out << "c Number of recursive call: " << m_nbCallCall << "\n";
    out << "c Number of split formula: " << m_nbSplit << "\n";
    out << "c Number of decision: " << m_nbDecisionNode << "\n";
    out << "c Number of paritioner calls: " << m_callPartitioner << "\n";
    out << "c Number of bad cuts " << m_failed_cuts << "\n";

    out << "c\n";
    m_cache->printCacheInformation(out);
    if (m_hCutSet) {
      out << "c\n";
      m_hCutSet->displayStat(out);
    }
    out << "c Final time: " << getTimer() << "\n";
    out << "c\n";
  } // printFinalStat

  /**
     Initialize the assumption in order to compute compiled formula under this
     one.

     @param[in] assums, the assumption
  */
  inline void initAssumption(std::vector<Lit> &assums) {
    m_solver->restart();
    m_solver->popAssumption(m_solver->getAssumption().size());
    m_solver->setAssumption(assums);
  } // initAssumption

  /**
     Decide if the cache is realized or not.
   */
  bool cacheIsActivated(std::vector<Var> &connected) {
    if (!m_optCached)
      return false;
    return m_cache->isActivated(connected.size());
  } // cacheIsActivated

  /**
     Call the CNF formula into a FBDD.

     @param[in] setOfVar, the current set of considered variables
     @param[in] unitsLit, the set of unit literal detected at this level
     @param[in] freeVariable, the variables which become free
     @param[in] out, the stream we use to print out information.

     \return an element of type U that sums up the given CNF sub-formula using a
     DPLL style algorithm with an operation manager.
  */
  bool onlyExtra(const std::vector<Var> &connected,
                 const std::vector<bool> &isDecisionVariable) {
    bool ok = true;
    for (int i = 0; i < connected.size(); i++) {
      ok &= !isDecisionVariable[connected[i]];
    }

    return ok;
  }
  U compute_(ProjVars &setOfVar, bool emergencyMode, std::vector<Lit> &unitsLit,
             std::vector<Var> &freeVariable, std::ostream &out) {
    showRun(out);
    m_nbCallCall++;
    int assume_nproj = reduceNProjVars(setOfVar);

    if (!m_solver->solve(setOfVar.vars)) {
      reduceNProjVarsPost(assume_nproj);
      return m_operation->manageBottom();
    }

    m_solver->whichAreUnits(setOfVar.vars, unitsLit); // collect unit literals
    m_specs->preUpdate(unitsLit);

    // compute the connected composant
    std::vector<ProjVars> varConnected;
    int nbComponent = m_specs->computeConnectedComponent(
        varConnected, setOfVar.vars, freeVariable);
    expelNoDecisionVar(freeVariable);
    varConnected.erase(
        std::partition(varConnected.begin(), varConnected.end(),
                       [&](const ProjVars &comp) { return comp.nbProj > 0; }),
        varConnected.end());

    std::sort(varConnected.begin(), varConnected.end(),
              [](ProjVars &a, ProjVars &b) {
                bool clean_a = a.vars.size() == a.nbProj;
                bool clean_b = b.vars.size() == b.nbProj;
                return clean_a > clean_b;
              });

    nbComponent = varConnected.size();

    // consider each connected component.
    if (nbComponent) {
      U tab[nbComponent];
      m_nbSplit += (nbComponent > 1) ? nbComponent : 0;
      for (int cp = 0; cp < nbComponent; cp++) {
        ProjVars &connected = varConnected[cp];
        bool cacheActivated = cacheIsActivated(connected.vars);

        TmpEntry<U> cb = cacheActivated ? m_cache->searchInCache(connected.vars)
                                        : NULL_CACHE_ENTRY;

        if (cacheActivated)
          nbTestCacheVarSize[connected.vars.size()]++;
        if (cacheActivated && cb.defined) {
          nbPosHitCacheVarSize[connected.vars.size()]++;
          tab[cp] = cb.getValue();
        } else {
          // recursive call
          tab[cp] = computeDecisionNode(connected, emergencyMode, out);

          if (cacheActivated)
            m_cache->addInCache(cb, tab[cp]);
        }
      }

      m_specs->postUpdate(unitsLit);
      expelNoDecisionLit(unitsLit);
      reduceNProjVarsPost(assume_nproj);

      return m_operation->manageDecomposableAnd(tab, nbComponent);
    } // else we have a tautology

    m_specs->postUpdate(unitsLit);
    expelNoDecisionLit(unitsLit);
    reduceNProjVarsPost(assume_nproj);

    return m_operation->createTop();
  } // compute_

  /**
   * @brief Set the Current Priority.
   *
   * @param cutSet is the set of variables that become decision variables.
   */
  inline void setCurrentPriority(std::vector<Var> &cutSet) {
    for (auto &v : cutSet)
      if (m_specs->isSelected(v))
        m_currentPrioritySet[v] = true;
  } // setCurrentPriority

  /**
   * @brief Unset the Current Priority.

   * @param cutSet is the set of variables that become decision variables.
   */
  inline void unsetCurrentPriority(std::vector<Var> &cutSet) {
    for (auto &v : cutSet)
      if (m_specs->isSelected(v))
        m_currentPrioritySet[v] = false;
  } // setCurrentPriority

  /**
     This function select a variable and compile a decision node.

     @param[in] connected, the set of variable present in the current problem.
     @param[in] out, the stream whare are printed out the logs.

     \return the compiled formula.
  */
  int test = true;
  U computeDecisionNode(ProjVars &connected, bool emergencyMode,
                        std::ostream &out) {
    std::vector<Var> cutSet;
    bool hasPriority = false, hasVariable = false;
    for (auto v : connected.vars) {
      if (m_specs->varIsAssigned(v) || !m_specs->isSelected(v))
        continue;
      hasVariable = true;
      if ((hasPriority = m_currentPrioritySet[v]))
        break;
    }
    bool cut_valid = false;

    if (hasVariable && !hasPriority && m_hCutSet->isReady(connected.vars)) {
      m_hCutSet->computeCutSet(connected.vars, cutSet);
      if(test){
          std::cout<<"Cut:"<<std::endl;
          for(auto i:cutSet){
              std::cout<<i<<" "<<std::endl;
          }
          std::cout<<std::endl;
          test = false;

      }

      for (int i = 0; i < cutSet.size(); i++) {
        cut_valid |= m_specs->isSelected(cutSet[i]) &&
                     !m_specs->varIsAssigned(cutSet[i]);
      }

      // cut_valid = cutset_ok(connected, cutSet);
      if (cut_valid) {
        setCurrentPriority(cutSet);

      } else if (!emergencyMode) {
        cutSet.clear();
        if (m_hBackUp->computeCutSetDyn(connected, cutSet)) {
          setCurrentPriority(cutSet);
          emergencyMode = true;
        }
      }
      m_callPartitioner++;
    }

    // search the next variable to branch on
    Var v;
    if (cut_valid || hasPriority) {
      v = m_hVar->selectVariable(connected.vars, *m_specs,
                                 m_currentPrioritySet);
    } else {
      m_failed_cuts += 1;
      v = m_hVar->selectVariable(connected.vars, *m_specs);
    }

    if (v == var_Undef) {
      unsetCurrentPriority(cutSet);
      return m_operation->manageTop(connected.vars);
    }

    assert(m_specs->isSelected(v));

    Lit l = Lit::makeLit(v, m_hPhase->selectPhase(v));
    m_nbDecisionNode++;

    // compile the formula where l is assigned to true
    DataBranch<U> b[2];

    assert(!m_solver->isInAssumption(l.var()));
    m_solver->pushAssumption(l);
    b[0].d =
        compute_(connected, emergencyMode, b[0].unitLits, b[0].freeVars, out);
    m_solver->popAssumption();

    if (m_solver->isInAssumption(l))
      b[1].d = m_operation->manageBottom();
    else if (m_solver->isInAssumption(~l))
      b[1].d =
          compute_(connected, emergencyMode, b[1].unitLits, b[1].freeVars, out);
    else {
      m_solver->pushAssumption(~l);
      b[1].d =
          compute_(connected, emergencyMode, b[1].unitLits, b[1].freeVars, out);
      m_solver->popAssumption();
    }

    unsetCurrentPriority(cutSet);
    return m_operation->manageDeterministOr(b, 2);
  } // computeDecisionNode

  /**
     Compute U using the trace of a SAT solver.

     @param[in] setOfVar, the set of variables of the considered problem.
     @param[in] out, the stream are is print out the logs.
     @param[in] warmStart, to activate/deactivate the warm start strategy.
     /!\ When the warm start is activated the assumptions are reset.

     \return an element of type U that sums up the given CNF formula using a
     DPLL style algorithm with an operation manager.
  */
  U compute(ProjVars &setOfVar, std::ostream &out, bool warmStart = true) {
    if (m_problem->isUnsat() ||
        (warmStart && !m_panicMode &&
         !m_solver->warmStart(29, 11, setOfVar.vars, m_out)))
      return m_operation->manageBottom();

    DataBranch<U> b;
    auto start = std::chrono::steady_clock::now();
    b.d = compute_(setOfVar, false, b.unitLits, b.freeVars, out);
    std::cout << "Elapsed(ms)=" << since(start).count() << std::endl;
    return m_operation->manageBranch(b);
  } // compute

public:
  /**
     Given an assumption, we compute the number of models.  That is different
     from the query strategy, where we first compute and then condition the
     computed structure.

     @param[in] setOfVar, the set of variables of the considered problem.
     @param[in] assumption, the set of literals we want to assign.
     @param[in] out, the stream where are print out the log.

     \return the number of models when the formula is simplified by the given
     assumption.
   */
  T count(std::vector<Var> &setOfVar, std::vector<Lit> &assumption,
          std::ostream &out) {
    initAssumption(assumption);

    // get the unit not in setOfVar.
    std::vector<Lit> shadowUnits;
    m_stampIdx++;
    for (auto &v : setOfVar)
      m_stampVar[v] = m_stampIdx;
    for (auto &l : assumption)
      if (m_stampVar[l.var()] != m_stampIdx)
        shadowUnits.push_back(l);

    m_specs->preUpdate(shadowUnits);

    ProjVars setOfVar_{setOfVar, m_specs->nbSelected()};
    U result = compute(setOfVar_, out, false);
    m_specs->postUpdate(shadowUnits);

    return m_operation->count(result);
  } // count

  /**
     Run the DPLL style algorithm with the operation manager.

     @param[in] vm, the set of options.
   */
  void run(po::variables_map &vm) {
    ProjVars setOfVar;
    setOfVar.nbProj = m_specs->nbSelected();
    for (int i = 1; i <= m_specs->getNbVariable(); i++)
      setOfVar.vars.push_back(i);

    U result = compute(setOfVar, m_out);
    printFinalStats(m_out);
    m_operation->manageResult(result, vm, m_out);
  } // run
};
} // namespace d4
