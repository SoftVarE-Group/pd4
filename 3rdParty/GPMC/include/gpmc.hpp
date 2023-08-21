#include <span>
#include <vector>
#include "src/problem/ProblemTypes.hpp"
namespace GPMC{
    bool simplify(std::vector<std::vector<d4::Lit>>& clauses,std::vector<d4::Var>& sel,std::vector<std::vector<d4::Lit>>& learnts,int& vars,int& pvars,
            std::vector<d4::Lit>& assigns,std::vector<d4::Lit>&gmap,std::vector<double>& activity);


}
