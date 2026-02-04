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
  if (GetExcludedRuns().count(run)) return false;
  if (run <= 0 || ev <= 0) return true;  // permissive on missing meta
  std::lock_guard<std::mutex> lock(GetMutex());
  return GetQADB().Pass(run, ev);
}

void QADBCuts::SetAllowedMiscRuns(const std::vector<int>& runs) {
  std::lock_guard<std::mutex> lock(GetMutex());
  QA::QADB& qa = GetQADB();
  for (int run : runs) qa.AllowMiscBit(run);
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

// ---------- helpers for "scalar-like" vectors ----------

int QADBCuts::FirstOrMinus1(const std::vector<int>& v) {
  return v.empty() ? -1 : v.front();
}

int QADBCuts::FirstOrMinus1(const ROOT::VecOps::RVec<int>& v) {
  return v.empty() ? -1 : v[0];
}

int QADBCuts::FirstOrMinus1(const std::vector<Long64_t>& v) {
  return v.empty() ? -1 : static_cast<int>(v.front());
}

int QADBCuts::FirstOrMinus1(const ROOT::VecOps::RVec<Long64_t>& v) {
  return v.empty() ? -1 : static_cast<int>(v[0]);
}

// ---------- functor overloads (RDataFrame compatible) ----------

bool QADBCuts::operator()(int run, int ev) const {
  return fAccumulateCharge ? PassAndAccumulate(run, ev) : Pass(run, ev);
}

bool QADBCuts::operator()(const std::vector<int>& run_vec,
                          const std::vector<int>& ev_vec) const {
  const int run = FirstOrMinus1(run_vec);
  const int ev  = FirstOrMinus1(ev_vec);
  return fAccumulateCharge ? PassAndAccumulate(run, ev) : Pass(run, ev);
}

bool QADBCuts::operator()(const ROOT::VecOps::RVec<int>& run_vec,
                          const ROOT::VecOps::RVec<int>& ev_vec) const {
  const int run = FirstOrMinus1(run_vec);
  const int ev  = FirstOrMinus1(ev_vec);
  return fAccumulateCharge ? PassAndAccumulate(run, ev) : Pass(run, ev);
}

bool QADBCuts::operator()(const std::vector<Long64_t>& run_vec,
                          const std::vector<Long64_t>& ev_vec) const {
  const int run = FirstOrMinus1(run_vec);
  const int ev  = FirstOrMinus1(ev_vec);
  return fAccumulateCharge ? PassAndAccumulate(run, ev) : Pass(run, ev);
}

bool QADBCuts::operator()(const ROOT::VecOps::RVec<Long64_t>& run_vec,
                          const ROOT::VecOps::RVec<Long64_t>& ev_vec) const {
  const int run = FirstOrMinus1(run_vec);
  const int ev  = FirstOrMinus1(ev_vec);
  return fAccumulateCharge ? PassAndAccumulate(run, ev) : Pass(run, ev);
}

// ---------- free helper ----------

bool QAOk(int run, int ev) {
  return QADBCuts::Pass(run, ev);
}
