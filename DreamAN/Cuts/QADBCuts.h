#ifndef QADBCUTS_H_
#define QADBCUTS_H_

#include <mutex>
#include <set>
#include <string>
#include <vector>

// For Long64_t
#include "RtypesCore.h"

// Needed for ROOT::VecOps::RVec (common with hipoDS + RDataFrame)
#include "ROOT/RVec.hxx"

#include "QADB.h"  // or "qadb/QADB.h"

// QADB-based event filter with configurable defects and run veto.
// Compatible with:
//   - hipo2root output (branches are std::vector<T>)
//   - hipoDS direct RDF (often ROOT::VecOps::RVec<T>, sometimes scalars)
//
// Important: ROOT's RDataFrame CallableTraits cannot deduce return types from
// a purely-templated operator(). Therefore we provide concrete overloads for
// the common column types.
class QADBCuts {
 public:
  QADBCuts();
  QADBCuts(const QADBCuts& other);
  ~QADBCuts();

  // --- Functor overloads for RDataFrame ---
  bool operator()(int run, int ev) const;

  bool operator()(const std::vector<int>& run_vec,
                  const std::vector<int>& ev_vec) const;
  bool operator()(const ROOT::VecOps::RVec<int>& run_vec,
                  const ROOT::VecOps::RVec<int>& ev_vec) const;

  bool operator()(const std::vector<Long64_t>& run_vec,
                  const std::vector<Long64_t>& ev_vec) const;
  bool operator()(const ROOT::VecOps::RVec<Long64_t>& run_vec,
                  const ROOT::VecOps::RVec<Long64_t>& ev_vec) const;

  // --- Configuration ---
  void SetDefects(const std::vector<std::string>& defects);

  // Convenience overload for const char* array:
  template <std::size_t N>
  void SetDefects(const char* const (&defs)[N]) {
    std::vector<std::string> v;
    v.reserve(N);
    for (std::size_t i = 0; i < N; ++i) v.emplace_back(defs[i]);
    SetDefects(v);
  }

  void AddDefect(const std::string& name);
  void AddExcludedRun(int run);
  void SetExcludedRuns(const std::set<int>& runs);
  void ClearExcludedRuns();

  // Allow misc-bit for runs
  static void SetAllowedMiscRuns(const std::vector<int>& runs);

  void SetAccumulateCharge(bool on);
  bool GetAccumulateCharge() const { return fAccumulateCharge; }

  static double GetAccumulatedCharge();
  static void ResetAccumulatedCharge();

  static bool Pass(int run, int ev);               // no charge
  static bool PassAndAccumulate(int run, int ev);  // with charge
  static void AccumulateCharge();

 private:
  bool fAccumulateCharge{true};

  static QA::QADB& GetQADB();
  static std::set<int>& GetExcludedRuns();
  static std::set<std::string>& GetDefectSet();
  static std::mutex& GetMutex();

  // Helpers to interpret "scalar-like" vectors (size 0 or 1 typical).
  // Empty -> -1, then Pass() remains permissive (as your code intended).
  static int FirstOrMinus1(const std::vector<int>& v);
  static int FirstOrMinus1(const ROOT::VecOps::RVec<int>& v);
  static int FirstOrMinus1(const std::vector<Long64_t>& v);
  static int FirstOrMinus1(const ROOT::VecOps::RVec<Long64_t>& v);
};

// Optional helper, similar to the old QAOk(run,ev) style.
bool QAOk(int run, int ev);

#endif  // QADBCUTS_H_
