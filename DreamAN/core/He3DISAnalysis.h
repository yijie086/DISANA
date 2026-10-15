#ifndef HE3DISANALYSIS_H
#define HE3DISANALYSIS_H

#include "DVCSAnalysis.h"

// Dedicated He3 DIS task. It reuses the established RDF, fiducial-cut and
// snapshot pipeline without enabling proton-DVCS momentum corrections.
class He3DISAnalysis : public DVCSAnalysis {
 public:
  He3DISAnalysis(bool isMC = false, bool isReproc = false,
                 bool isMinBook = false);
  ~He3DISAnalysis() override = default;
};

#endif
