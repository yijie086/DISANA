#ifndef HE3DVCSANALYSIS_H
#define HE3DVCSANALYSIS_H

#include "DVCSAnalysis.h"

// Dedicated He3 DVCS task type. It intentionally reuses only the common
// event-selection pipeline; pDVCS momentum/W corrections are not configured.
class He3DVCSAnalysis : public DVCSAnalysis {
 public:
  He3DVCSAnalysis(bool isMC = false, bool isReproc = false, bool isMinBook = false);
  ~He3DVCSAnalysis() override = default;
};

#endif
