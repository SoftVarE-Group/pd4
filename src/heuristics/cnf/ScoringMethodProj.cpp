
#include "ScoringMethodProj.hpp"

namespace d4 {

/**
   Constructor.

   @param[in] o, the specification of a CNF problem.
   @param[in] a, an activity manager linked to a solver.
 */
ScoringMethodProj::ScoringMethodProj(SpecManagerCnf &o, ActivityManager &a)
    : om(o), activity(a) {}  // constructor

/**
   The well-known VSADS heuristic.

   Sang, T.; Beame, P.; and Kautz, H. Heuristics for Fast
   Exact Model Counting. In Proceedings of the 8th International
   Conference on Theory and Applications of Satisfiability Testing.

   @param[in] v, the variable we want the score.
 */
double ScoringMethodProj::computeScore(Var v) {
  return activity.getCountConflict(v) + (double)(om.getNbClause(v) << 7);
}

}  // namespace d4
