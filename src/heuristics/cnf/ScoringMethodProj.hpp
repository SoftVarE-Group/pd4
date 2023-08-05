#include <src/solvers/ActivityManager.hpp>
#include <src/specs/cnf/SpecManagerCnf.hpp>

#include "../ScoringMethod.hpp"

namespace d4 {
class ScoringMethodProj : public ScoringMethod {
 private:
  SpecManagerCnf &om;
  ActivityManager &activity;

 public:
  ScoringMethodProj(SpecManagerCnf &o, ActivityManager &a);
  double computeScore(Var v);
};
}  // namespace d4
