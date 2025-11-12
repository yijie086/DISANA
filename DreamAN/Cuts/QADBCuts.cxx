#include <utility>

#include "QADBCuts.h"

// ---------- static helpers for shared state ----------

QA::QADB& QADBCuts::GetQADB() {
  static QA::QADB instance("latest");  // single shared instance
  return instance;
}

std::set<int>& QADBCuts::GetExcludedRuns() {
  static std::set<int> runs;
  return runs;
}

std::set<std::string>& QADBCuts::GetDefectSet() {
  static std::set<std::string> defects;
  return defects;
}

std::mutex& QADBCuts::GetMutex() {
  static std::mutex m;
  return m;
}

// ---------- constructors / destructor ----------

QADBCuts::QADBCuts() : fAccumulateCharge(true) {}

QADBCuts::QADBCuts(const QADBCuts& other) : fAccumulateCharge(other.fAccumulateCharge) {}

QADBCuts::~QADBCuts() = default;

// ---------- configuration methods ----------

void QADBCuts::SetDefects(const std::vector<std::string>& defects) {
  std::lock_guard<std::mutex> lock(GetMutex());
  std::set<std::string>& defset = GetDefectSet();
  QA::QADB& qa = GetQADB();

  // Add new defects; we do NOT try to remove old ones (QADB has no "unset").
  // Intended: call once at startup with your full list.
  for (const auto& d : defects) {
    if (defset.insert(d).second) {
      qa.CheckForDefect(d.c_str());
    }
  }
}

void QADBCuts::AddDefect(const std::string& name) {
  std::lock_guard<std::mutex> lock(GetMutex());
  std::set<std::string>& defset = GetDefectSet();
  if (defset.insert(name).second) {
    GetQADB().CheckForDefect(name.c_str());
  }
}

void QADBCuts::AddExcludedRun(int run) {
  std::lock_guard<std::mutex> lock(GetMutex());
  GetExcludedRuns().insert(run);
}

void QADBCuts::SetExcludedRuns(const std::set<int>& runs) {
  std::lock_guard<std::mutex> lock(GetMutex());
  GetExcludedRuns() = runs;
}

void QADBCuts::ClearExcludedRuns() {
  std::lock_guard<std::mutex> lock(GetMutex());
  GetExcludedRuns().clear();
}

void QADBCuts::SetAccumulateCharge(bool on) { fAccumulateCharge = on; }

// ---------- static QA / charge helpers ----------

bool QADBCuts::Pass(int run, int ev) {
  // no charge accumulation
  if (GetExcludedRuns().count(run)) return false;
  if (run <= 0 || ev <= 0) return true;  // permissive on missing meta
  std::lock_guard<std::mutex> lock(GetMutex());
  return GetQADB().Pass(run, ev);
}

bool QADBCuts::PassAndAccumulate(int run, int ev) {
  if (GetExcludedRuns().count(run)) return false;
  if (run <= 0 || ev <= 0) return true;
  std::lock_guard<std::mutex> lock(GetMutex());
  if (GetQADB().Pass(run, ev)) {
    GetQADB().AccumulateCharge();
    return true;
  }
  return false;
}

void QADBCuts::AccumulateCharge() {
  std::lock_guard<std::mutex> lock(GetMutex());
  GetQADB().AccumulateCharge();
}

double QADBCuts::GetAccumulatedCharge() {
  std::lock_guard<std::mutex> lock(GetMutex());
  return GetQADB().GetAccumulatedCharge();
}

void QADBCuts::ResetAccumulatedCharge() {
  std::lock_guard<std::mutex> lock(GetMutex());
  GetQADB().ResetAccumulatedCharge();
}

// ---------- functor implementation ----------

bool QADBCuts::operator()(int run, int ev) const { return fAccumulateCharge ? PassAndAccumulate(run, ev) : Pass(run, ev); }

// ---------- free helper ----------

bool QAOk(int run, int ev) {
  // By default, don't accumulate charge in this helper.
  return QADBCuts::Pass(run, ev);
}
