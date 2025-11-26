#ifndef QADBCUTS_H_
#define QADBCUTS_H_

#include <mutex>
#include <set>
#include <string>
#include <vector>

#include "QADB.h"  // or "qadb/QADB.h"

// QADB-based event filter with configurable defects and run veto.
// Internally uses a single shared QA::QADB instance (static).
class QADBCuts {
 public:
  QADBCuts();
  QADBCuts(const QADBCuts& other);
  ~QADBCuts();

  // Functor used by RDataFrame: returns true if event passes QA.
  bool operator()(int run, int ev) const;

  // ---- configuration (call from RunDVCSAnalysis.C, etc.) ----

  // Replace the entire set of defects to check (intended to be called once).
  void SetDefects(const std::vector<std::string>& defects);

  // Convenience overload for const char* array:
  template <std::size_t N>
  void SetDefects(const char* const (&defs)[N]) {
    std::vector<std::string> v;
    v.reserve(N);
    for (std::size_t i = 0; i < N; ++i) {
      v.emplace_back(defs[i]);
    }
    SetDefects(v);
  }

  // Add one more defect to check (on top of those already configured).
  void AddDefect(const std::string& name);

  // Exclude a specific run entirely.
  void AddExcludedRun(int run);

  // Replace the full set of excluded runs.
  void SetExcludedRuns(const std::set<int>& runs);

  // Clear all excluded runs.
  void ClearExcludedRuns();

  // Set runs for which the 'Misc' bit defect should be explicitly allowed to pass.
  // This is equivalent to calling QA::QADB::allowMiscBit() for each run.
  static void SetAllowedMiscRuns(const std::vector<int>& runs); // sometimes misc are allowed

  // Control whether successful events also accumulate charge.
  void SetAccumulateCharge(bool on);
  bool GetAccumulateCharge() const { return fAccumulateCharge; }

  // ---- global charge helpers for the shared QADB instance ----

  static double GetAccumulatedCharge();
  static void ResetAccumulatedCharge();

  // Static QA helpers (if you prefer to call them directly).
  static bool Pass(int run, int ev);               // no charge
  static bool PassAndAccumulate(int run, int ev);  // with charge
  static void AccumulateCharge();

 private:
  bool fAccumulateCharge;  // if true, operator() also accumulates charge

  // Internal helpers to access shared state.
  static QA::QADB& GetQADB();
  static std::set<int>& GetExcludedRuns();
  static std::set<std::string>& GetDefectSet();
  static std::mutex& GetMutex();
};

// Optional helper, similar to the old QAOk(run,ev) style.
bool QAOk(int run, int ev);

#endif  // QADBCUTS_H_
