#include <algorithm>

#pragma once

#include <boost/program_options.hpp>
#include <vector>

#include "../PreprocManager.hpp"
#include "util/PreprocGPMC.hpp"
#include "src/solvers/WrapperSolver.hpp"
//Preprocessor based on gpmc and sharpsat-td

namespace d4 {
namespace po = boost::program_options;


class PreprocProj : public PreprocManager {
private:
    WrapperSolver *ws;
    PRE::ConfigPreprocessor config;

public:
  PreprocProj(po::variables_map &vm, std::ostream &out);
  ~PreprocProj();
  virtual ProblemManager *run(ProblemManager *pin,
                              LastBreathPreproc &lastBreath) override;

};
} // namespace d4
