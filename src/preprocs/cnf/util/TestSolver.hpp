#pragma once
#include <vector>
#include "3rdParty/glucose-3.0/core/Solver.h"
//This code was stolen straight from gpmc
namespace d4 {
class Identifier {
public:
	Identifier(int vars) { cidx.resize(2*vars, -1); num_elem = 0; }

	void identify(Glucose::Lit l1, Glucose::Lit l2);

	bool hasEquiv() { return num_elem > 0; }
	std::vector<std::vector<Glucose::Lit>>& getEquivClasses() { return eqc; }
	Glucose::Lit delegateLit(Glucose::Lit l) { return (cidx[toInt(l)] == -1) ? l : eqc[cidx[toInt(l)]][0]; }
	int getIndex(Glucose::Lit l) { return cidx[toInt(l)]; }
	void removeEquivClass(Glucose::Lit l);

private:
	void MergeEquivClasses(int c1, int c2);

	std::vector<std::vector<Glucose::Lit>> eqc;
	std::vector<int> cidx;
	int num_elem;
};
class TestSolver : public Glucose::Solver {
public:
  TestSolver(int nvars) { newVars(nvars); }
  TestSolver(int nvars, std::vector<std::vector<Glucose::Lit>> clauses,
             std::vector<std::vector<Glucose::Lit>> learnts,
             std::vector<Glucose::Lit> assignedLits);

  void addClauseWith(const std::vector<Glucose::Lit> &ps, bool learnt = false);
  void resetClauses(std::vector<std::vector<Glucose::Lit>> &clauses);

  Glucose::lbool Solve() { return solve_(); }
  Glucose::lbool Solve(const Glucose::vec<Glucose::Lit> &assumptions) {
    budgetOff();
    setConfBudget(clauses.size() * 10);
    return solveLimited(assumptions);
  }

  bool falsifiedBy(Glucose::Lit l);
  bool falsifiedBy(Glucose::Lit l1, Glucose::Lit l2);
  bool falsifiedBy(Glucose::vec<Glucose::Lit> &assump);
  bool FailedLiterals();

  void exportLearnts(std::vector<std::vector<Glucose::Lit>> &learnts);

  void assign(Glucose::Lit l) { enqueue(l); }
  bool bcp() { return propagate() == Glucose::CRef_Undef; }

  Glucose::vec<Glucose::lbool> &getAssigns() { return assigns; }
  Glucose::vec<Glucose::Lit> &getTrail() { return trail; }

private:
  void newVars(int nvars);
  void backTo(int pos);
};
} // namespace d4
