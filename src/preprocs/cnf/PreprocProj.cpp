
#include "PreprocProj.hpp"

#include "3rdParty/glucose-3.0/utils/System.h"
#include "src/problem/cnf/ProblemManagerCnf.hpp"

const Glucose::Var var_Determined = {-2};
namespace d4 {

void PreprocProj::merge_equiv() {
  int befor_cl = clauses.size();
  int befor_var = vars;
  Identifier Id(vars);
  {
    std::vector<sspp::Bitset> adjmat;
    adjmat.resize(vars);
    for (int i = 0; i < vars; i++)
      adjmat[i] = sspp::Bitset(vars);

    for (const auto &cls : {clauses, learnts}) {
      for (const auto &clause : cls)
        for (int i = 0; i < clause.size(); i++)
          for (int j = i + 1; j < clause.size(); j++) {
            Var v1 = var(clause[i]);
            Var v2 = var(clause[j]);
            if (!adjmat[v1].Get(v2)) {
              adjmat[v1].SetTrue(v2);
              adjmat[v2].SetTrue(v1);
            }
          }
    }

    //
    // assert: any projected var idx < any non-projceted var idx
    //

    TestSolver S(vars, clauses, learnts, assignedLits);
    Glucose::vec<Glucose::Lit> &trail = S.getTrail();

    for (Var v1 = 0; v1 < vars; v1++) {
      for (Var v2 = v1 + 1; v2 < vars; v2++) {
        if (!adjmat[v1].Get(v2))
          continue;

        Glucose::Lit l1 = Glucose::mkLit(v1);
        Glucose::Lit l2 = Glucose::mkLit(v2);
        int idx1, idx2;

        if (S.value(l1) != Glucose::l_Undef)
          break;
        if (S.value(l2) != Glucose::l_Undef)
          continue;

        int pos = trail.size();

        if (S.falsifiedBy(l1, l2) && S.falsifiedBy(~l1, ~l2)) {
          Id.identify(l1, ~l2);
        } else if (S.falsifiedBy(l1, ~l2) && S.falsifiedBy(~l1, l2)) {
          Id.identify(l1, l2);
        }

        if (trail.size() - pos > 0) {
          int ppos = pos;
          int cpos = trail.size();

          while (ppos < cpos) {
            for (int i = ppos; i < cpos; i++) {
              int idx = Id.getIndex(trail[i]);
              if (idx == -1) {
                S.assign(trail[i]);
              } else {
                for (Glucose::Lit l : Id.getEquivClasses()[idx])
                  S.assign(l);
              }
              S.bcp();

              Id.removeEquivClass(trail[i]);
            }
            ppos = cpos;
            cpos = trail.size();
          }
          for (int i = pos; i < trail.size(); i++)
            assignedLits.push_back(trail[i]);
        }
      }
    }
  }

  bool subsump = false;
  if (Id.hasEquiv()) {
    std::vector<Glucose::Lit> map(2 * vars);
    int new_id = 0;
    for (int i = 0; i < vars; i++) {
      Glucose::Lit dl = Id.delegateLit(Glucose::mkLit(i));
      Glucose::Lit l = Glucose::mkLit(i);

      if (Glucose::mkLit(i) == dl) {
        Glucose::Lit newlit = Glucose::mkLit(new_id);
        map[toInt(l)] = newlit;
        map[toInt(~l)] = ~newlit;
        new_id++;
      } else {
        Glucose::Lit newlit = map[toInt(dl)];
        map[toInt(l)] = newlit;
        map[toInt(~l)] = ~newlit;
      }
      if (i == pvars - 1)
        pvars = new_id;
    }
    vars = new_id;
    ispvars.clear();
    ispvars.resize(pvars, true);
    ispvars.resize(vars, false);

    if (assignedLits.size() > 0) {
      std::vector<Glucose::Lit> assignedLits_new;
      for (Glucose::Lit l : assignedLits)
        assignedLits_new.push_back(map[toInt(l)]);
      assignedLits = assignedLits_new;
    }

    rewrite_claues(clauses, map);
    rewrite_claues(learnts, map);

    // update the variable map
    for (int i = 0; i < gmap.size(); i++) {
      Glucose::Lit l = gmap[i];
      if (l == Glucose::lit_Undef)
        continue;
      gmap[i] = map[toInt(l)];
    }

    subsump = true;
  }

  if (assignedLits.size() > 0) {
    sspp::SortAndDedup(assignedLits);

    TestSolver S(vars, clauses, learnts, assignedLits);
    compact(S.getAssigns());
    sspp::SortAndDedup(clauses);
    sspp::SortAndDedup(learnts);

    subsump = true;
  }

  if (subsump)
    subsume();
  stats.elim_adj_var += befor_var - vars;
  stats.elim_adj += befor_cl - clauses.size();
}

/**
   The constructor.

   @param[in] vm, the options used (solver).
 */
PreprocProj::PreprocProj(po::variables_map &vm, std::ostream &out) {

  ws = WrapperSolver::makeWrapperSolverPreproc(vm, out);
} // constructor

/**
   Destructor.
 */
PreprocProj::~PreprocProj() { delete ws; } // destructor

/**
 * @brief The preprocessing itself.
 * @param[out] p, the problem we want to preprocess.
 * @param[out] lastBreath gives information about the way the    preproc sees
 * the problem.
 */

bool PreprocProj::sat_fle() {
  TestSolver S(vars, clauses, learnts, assignedLits);
  if (!S.okay())
    return false;

  if (!S.FailedLiterals() || S.Solve() == Glucose::l_False) {
    unsat = true;
    return false;
  }

  S.exportLearnts(learnts);
  compact(S.getAssigns());
  sspp::SortAndDedup(clauses);
  sspp::SortAndDedup(learnts);
  subsume();
  return true;
}
void PreprocProj::strenthen() {
  TestSolver S(vars, clauses, learnts, assignedLits);
  int befor = clauses.size();
  bool removecl = false;
  bool removelit = false;
  for (int i = clauses.size() - 1; i >= 0; i--) {
    Glucose::vec<Glucose::Lit> assump;

    for (int j = clauses[i].size() - 1; j >= 0; j--) {
      assump.clear();
      removelit = false;

      for (int k = 0; k < clauses[i].size(); k++) {
        Glucose::Lit l = clauses[i][k];
        if (S.value(l) == Glucose::l_True) {
          removecl = true;
          break;
        }
        if (k == j) {
          if (S.value(l) == Glucose::l_False) {
            removelit = true;
            break;
          }
        } else if (S.value(l) == Glucose::l_Undef) {
          assump.push(~l);
        }
      }
      if (removecl)
        break;

      if (removelit || S.falsifiedBy(assump)) {
        for (int k = j + 1; k < clauses[i].size(); k++)
          clauses[i][k - 1] = clauses[i][k];
        clauses[i].pop_back();

        if (clauses[i].size() == 1) {
          S.assign(clauses[i][0]);
          S.bcp(); // no conflict because SAT test was already passed.
          removecl = true;
          break;
        }
      }
    }
    if (removecl) {
      sspp::SwapDel(clauses, i);
      removecl = false;
    }
  }

  S.exportLearnts(learnts);
  compact(S.getAssigns());
  sspp::SortAndDedup(clauses);
  sspp::SortAndDedup(learnts);
  subsume();
  stats.elim_stren += befor - clauses.size();
}
bool PreprocProj::elim_mono() {
  int befor = clauses.size();
  std::vector<int> monolit(vars);
  for (auto &cl : clauses) {
    for (auto l : cl) {
      int ll = sign(l) + 1;
      if (monolit[var(l)] == 0) {
        monolit[var(l)] = ll;
      } else if (monolit[var(l)] != ll) {
        monolit[var(l)] = -1;
      }
    }
  }
  std::vector<bool> is_mono_var(vars, false);
  int cnt = 0;
  for (int i = pvars; i < vars; i++) {
    if (monolit[i] == 1 || monolit[i] == 2) {
      is_mono_var[i] = true;
      cnt++;
    }
  }
  if (cnt == 0) {
    return false;
  }
  for (int i = clauses.size() - 1; i >= 0; i--) {
    bool removecl = false;
    for (auto l : clauses[i]) {
      if (is_mono_var[var(l)]) {
        removecl = true;
        break;
      }
    }
    if (removecl) {
      sspp::SwapDel(clauses, i);
    }
  }
  stats.elim_mono += befor - clauses.size();
  return true;
}

static inline Glucose::lbool val(const Glucose::vec<Glucose::lbool> &assigns,
                                 Glucose::Lit p) {
  return assigns[Glucose::var(p)] ^ Glucose::sign(p);
}
void PreprocProj::compact_clauses(const Glucose::vec<Glucose::lbool> &assigns,
                                  std::vector<std::vector<Glucose::Lit>> &cls,
                                  std::vector<bool> &occurred, int &varnum) {
  int i1, i2;
  int j1, j2;

  for (i1 = 0, i2 = 0; i1 < cls.size(); i1++) {
    std::vector<Glucose::Lit> &c = cls[i1];
    for (j1 = 0, j2 = 0; j1 < c.size(); j1++) {
      if (val(assigns, c[j1]) == Glucose::l_Undef)
        c[j2++] = c[j1];
      else if (val(assigns, c[j1]) == Glucose::l_True) {
        goto NEXTC;
      }
    }
    c.resize(j2);
    assert(c.size() == j2);
    assert(c.size() > 1);

    for (auto l : c) {
      Var v = var(l);
      if (!occurred[v]) {
        occurred[v] = true;
        varnum++;
      }
    }

    cls[i2++] = cls[i1];
  NEXTC:;
  }
  cls.resize(i2);
}
void PreprocProj::compact(const Glucose::vec<Glucose::lbool> &assigns,
                          const std::vector<Glucose::Var> &elimvars) {
  int varnum = 0;
  std::vector<bool> occurred;
  occurred.resize(vars, false);

  // Compact Clauses
  compact_clauses(assigns, clauses, occurred, varnum);
  compact_clauses(assigns, learnts, occurred, varnum);

  // Compact Variables
  int new_idx = 0;
  std::vector<Var> map;
  std::vector<Var> nonpvars;

  map.clear();
  map.resize(vars, var_Undef);
  for (auto v : elimvars)
    if (v < pvars)
      map[v] = var_Determined;

  nonpvars.reserve(vars - pvars);

  for (Var v = 0; v < vars; v++) {
    if (occurred[v]) {
      if (ispvars[v]) {
        map[v] = new_idx;
        new_idx++;
      } else
        nonpvars.push_back(v);
    } else {
      if (ispvars[v]) {
        if (assigns[v] == Glucose::l_Undef) {
          freevars++;
        } else {
        }
      }
    }
  }
  for (int i = 0; i < pvars; i++)
    gmap.push_back(Glucose::mkLit(i));
  pvars = new_idx;
  for (int i = 0; i < nonpvars.size(); i++) {
    map[nonpvars[i]] = new_idx;
    new_idx++;
  }
  vars = new_idx;

  ispvars.clear();
  ispvars.resize(pvars, true);
  ispvars.resize(vars, false);

  // Replace literals according to map
  rewrite_claues(clauses, map);
  rewrite_claues(learnts, map);

  // update the variable map
  std::unordered_map<Glucose::Var, int> table;

  for (int i = 0; i < gmap.size(); i++) {
    Glucose::Lit l = gmap[i];
    if (l == Glucose::lit_Undef)
      continue;

    Var x = var(gmap[i]);
    if (assigns[x] == Glucose::l_Undef) {
      Var newv = map[Glucose::var(l)];
      if (newv >= 0) {
        gmap[i] = Glucose::mkLit(newv, sign(l));
      } else {
        assert(newv == var_Undef);
        auto itr = table.find(x);
        if (itr == table.end()) {
          table[x] = freeLitClasses.size();
          freeLitClasses.push_back({Glucose::mkLit(i, sign(l))});
        } else {
          freeLitClasses[itr->second].push_back(Glucose::mkLit(i, sign(l)));
        }
        gmap[i] = Glucose::lit_Undef;
      }
    } else if (map[var(l)] == var_Determined) {
      definedVars.push_back(i);
      gmap[i] = Glucose::lit_Undef;
    } else {
      bool lsign = ((assigns[x] == Glucose::l_False) != sign(l));
      fixedLits.push_back(Glucose::mkLit(i, lsign));
      gmap[i] = Glucose::lit_Undef;
    }
  }
  assignedLits.clear();
}

void PreprocProj::subsume() {
  {
    sspp::Subsumer subsumer;
    clauses = subsumer.Subsume(clauses);
    sspp::SortAndDedup(clauses);
  }

  if (learnts.empty())
    return;
  for (const auto &clause : clauses) {
    learnts.push_back(clause);
  }

  {
    sspp::Subsumer subsumer_lrnt;
    learnts = subsumer_lrnt.Subsume(learnts);
    sspp::SortAndDedup(learnts);
  }

  for (int i = 0; i < learnts.size(); i++) {
    if (std::binary_search(clauses.begin(), clauses.end(), learnts[i])) {
      sspp::SwapDel(learnts, i);
      i--;
    }
  }
}
void PreprocProj::rewrite_claues(std::vector<std::vector<Glucose::Lit>> &cls,
                                 const std::vector<Glucose::Lit> &map) {
  // map may not be injective, i.e., this may strengthen clauses.
  for (int i = 0; i < cls.size(); i++) {
    std::vector<Glucose::Lit> &c = cls[i];
    for (int j = 0; j < c.size(); j++)
      c[j] = map[toInt(c[j])];

    sspp::SortAndDedup(c);

    bool unit = false;
    bool taut = false;
    if (c.size() == 1) {
      assigns[var(c[0])] =
          Glucose::sign(c[0]) ? Glucose::l_False : Glucose::l_True;
      assignedLits.push_back(c[0]);
      unit = true;
    } else {
      for (int j = 1; j < c.size(); j++)
        if (var(c[j]) == var(c[j - 1])) {
          taut = true;
          break;
        }
    }
    if (unit || taut) {
      sspp::SwapDel(cls, i);
      i--;
    }
  }
}
void PreprocProj::rewrite_claues(std::vector<std::vector<Glucose::Lit>> &cls,
                                 const std::vector<Glucose::Var> &map) {

  for (int i = 0; i < cls.size(); i++) {
    std::vector<Glucose::Lit> &c = cls[i];
    for (int j = 0; j < c.size(); j++)
      c[j] = Glucose::mkLit(map[var(c[j])], sign(c[j]));
    sort(c.begin(), c.end());
  }
}
std::vector<std::vector<Lit>>
conv_cl(std::vector<std::vector<Glucose::Lit>> &cl) {
  std::vector<std::vector<Lit>> out(cl.size());
  for (int i = 0; i < out.size(); i++) {
    out[i].resize(cl[i].size());
    for (int k = 0; k < cl[i].size(); k++) {
      out[i][k] =
          Lit::makeLit(Glucose::var(cl[i][k]) + 1, Glucose::sign(cl[i][k]));
    }
  }
  return out;
}
ProblemManager *PreprocProj::run(ProblemManager *pin,
                                 LastBreathPreproc &lastBreath) {
  ProblemManagerCnf *in = (ProblemManagerCnf *)pin;
  freevars = 0;
  unsat = in->isUnsat();
  vars = in->getNbVar();
  pvars = in->getNbSelectedVar();
  ispvars.resize(vars, false);
  assigns.resize(vars, Glucose::l_Undef);
  for (auto v : in->getSelectedVar()) {
    ispvars[v - 1] = true;
  }
  for (auto &cl : in->getClauses()) {
    std::vector<Glucose::Lit> lits(cl.size());
    for (int i = 0; i < cl.size(); i++) {
      lits[i] = Glucose::mkLit(cl[i].var() - 1, cl[i].sign());
    }
    add_clause(lits);
  }
  simplify();
  std::vector<Lit> units;
  std::vector<double> weight(vars + 1, 1.0);
  std::vector<double> weightLit((vars + 1) << 2, 1.0);
  for (unsigned i = 0; i <= vars; i++)
    weight[i] = weightLit[i << 1] + weightLit[(i << 1) + 1];
  std::vector<Var> selected;
  for (Var i = 1; i <= pvars; i++) {
    selected.push_back(i);
  }
  ProblemManagerCnf *out =
      new ProblemManagerCnf(vars, weightLit, weight, selected);
  out->getClauses() = conv_cl(clauses);
  lastBreath.learnt = conv_cl(learnts);

  ws->initSolver(*out);
  lastBreath.panic = 0;
  lastBreath.countConflict.resize(out->getNbVar() + 1, 0);
  if (!ws->solve())
    return out->getUnsatProblem();
  lastBreath.panic = ws->getNbConflict() > 100000;

  // get the activity given by the solver.
  for (unsigned i = 1; i <= out->getNbVar(); i++)
    lastBreath.countConflict[i] = ws->getCountConflict(i);
  return out;
} // run
Glucose::lbool PreprocProj::value(Glucose::Lit p) {
  return assigns[var(p)] ^ sign(p);
}
bool PreprocProj::add_clause(std::vector<Glucose::Lit> &ps) {
  sort(ps.begin(), ps.end());
  Glucose::vec<Glucose::Lit> oc;
  oc.clear();
  Glucose::Lit p;
  int i, j = 0;
  for (i = j = 0, p = Glucose::lit_Undef; i < ps.size(); i++)
    if (value(ps[i]) == Glucose::l_True || ps[i] == ~p)
      return true;
    else if (value(ps[i]) != Glucose::l_False && ps[i] != p)
      ps[j++] = p = ps[i];
  ps.resize(j);

  if (ps.size() == 0) {
    unsat = true;
    return false;
  } else if (ps.size() == 1) {
    assigns[var(ps[0])] = Glucose::lbool(!sign(p));
    assignedLits.push_back(ps[0]);
    // No check by UP here. Preprocessor will do later.
  } else {
    clauses.push_back({});
    copy(ps.begin(), ps.end(), back_inserter(clauses.back()));
  }
  return true;
}
void PreprocProj::simplify() {
  double start = Glucose::cpuTime();
  if (unsat || !sat_fle()) {
    return;
  }
  elim_mono();
  strenthen();
  for (int i = 0; i < config.reps; i++) {
    if (config.max_time > (start - Glucose::cpuTime())) {
      break;
    }
    int last_cls = clauses.size();
    int last_vars = vars;
    elim_mono();
    strenthen();
    //merge_equiv();
    std::cout << "Elim " << last_cls - clauses.size() << std::endl;
    if ((last_cls == clauses.size()) && (last_vars == vars)) {
      break;
    }
  }

  int clearnts = learnts.size();
  for (int i = learnts.size() - 1; i >= 0; i--) {
    std::vector<Glucose::Lit> &clause = learnts[i];
    if (clause.size() > 3)
      sspp::SwapDel(learnts, i);
  }

  if (config.max_time < (start - Glucose::cpuTime())) {
    sat_fle();
  }
  if (clearnts * 3 / 2 < learnts.size()) {
    if (vars > 0) {
      int m = clauses.size() / vars;
      int len = std::max(m + 1, 3);
      for (int i = learnts.size() - 1; i >= 0; i--) {
        std::vector<Glucose::Lit> &clause = learnts[i];
        if (clause.size() > len)
          sspp::SwapDel(learnts, i);
      }
    }
  }
  std::cout << "Elim Mono: " << stats.elim_mono << std::endl;
  std::cout << "Elim Strenthen: " << stats.elim_stren << std::endl;
  std::cout << "Elim Adj: " << stats.elim_adj << std::endl;
  std::cout << "Elim Adj Vars: " << stats.elim_adj_var << std::endl;
}
} // namespace d4
