#include "He3DVCSAnalysis.h"

He3DVCSAnalysis::He3DVCSAnalysis(bool isMC, bool isReproc, bool isMinBook)
    : DVCSAnalysis(isMC, isReproc, isMinBook) {
  // He3 analyses do not currently apply a momentum correction, so retain the
  // fiducial-selected RDF even in minimal-booking mode.
  SetWriteAfterFidWhenNoCorrection(true);
}
