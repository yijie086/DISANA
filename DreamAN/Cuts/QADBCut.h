#ifndef QADBCUT_H_
#define QADBCUT_H_

#include <mutex>
#include "QADB.h"  // or "qadb/QADB.h"

namespace qadbcfg {

// Thread-safe controller around one shared QADB instance.
// We lock Pass+Accumulate in ONE critical section to preserve bin context.
class QADBctl {
private:
  static QA::QADB& qa() {
    static QA::QADB inst("latest");   // single shared instance
    static bool inited = false;
    if (!inited) {
      inst.CheckForDefect("TotalOutlier");
      inst.CheckForDefect("TerminalOutlier");
      inst.CheckForDefect("MarginalOutlier");
      inst.CheckForDefect("SectorLoss");
      inst.CheckForDefect("LowLiveTime");
      inst.CheckForDefect("Misc");
      inst.CheckForDefect("TotalOutlierFT");
      inst.CheckForDefect("TerminalOutlierFT");
      inst.CheckForDefect("MarginalOutlierFT");
      inst.CheckForDefect("LossFT");
      inst.CheckForDefect("BSAWrong");
      inst.CheckForDefect("BSAUnknown");
      inst.CheckForDefect("TSAWrong");
      inst.CheckForDefect("TSAUnknown");
      inst.CheckForDefect("DSAWrong");
      inst.CheckForDefect("DSAUnknown");
      inst.CheckForDefect("ChargeHigh");
      inst.CheckForDefect("ChargeNegative");
      inst.CheckForDefect("ChargeUnknown");
      inst.CheckForDefect("PossiblyNoBeam");
      inited = true;
    }
    return inst;
  }
  static std::mutex& mtx() {
    static std::mutex m;
    return m;
  }

public:
  // Atomic version: Pass AND (if pass) AccumulateCharge in one lock.
  static inline bool PassAndAccumulate(int run, int ev) {
    if (run == 5740 || run == 3262) return false;           // your veto
    if (run <= 0 || ev <= 0) return true;    // permissive on missing meta
    std::lock_guard<std::mutex> lock(mtx());
    if (qa().Pass(run, ev)) {
      qa().AccumulateCharge();               // same bin context as Pass()
      return true;
    }
    return false;
  }

  // (Optional) separate ops — do NOT use them split under MT.
  static inline bool Pass(int run, int ev) {
    if (run == 5740 || run == 3262) return false;
    if (run <= 0 || ev <= 0) return true;
    std::lock_guard<std::mutex> lock(mtx());
    return qa().Pass(run, ev);
  }
  static inline void AccumulateCharge() {
    std::lock_guard<std::mutex> lock(mtx());
    qa().AccumulateCharge();
  }

  static inline double GetAccumulatedCharge() {
    std::lock_guard<std::mutex> lock(mtx());
    return qa().GetAccumulatedCharge();      // typically nC
  }
  static inline void ResetAccumulatedCharge() {
    std::lock_guard<std::mutex> lock(mtx());
    qa().ResetAccumulatedCharge();
  }
};

// Backward-compatible wrapper if you still call QAOk(run,ev)
inline bool QAOk(int run, int ev) { return QADBctl::Pass(run, ev); }

} // namespace qadbcfg

// EventCut-style functor (unchanged)
class QADBCut {
public:
  inline bool operator()(int run, int ev) const { return qadbcfg::QAOk(run, ev); }
};

#endif  // QADBCUT_H_
