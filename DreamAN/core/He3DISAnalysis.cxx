#include "He3DISAnalysis.h"

He3DISAnalysis::He3DISAnalysis(bool isMC, bool isReproc, bool isMinBook)
    : DVCSAnalysis(isMC, isReproc, isMinBook) {
  // He3 DIS has no momentum-correction stage at present, but its
  // fiducial-selected RDF is a required minimal-booking output.
  SetWriteAfterFidWhenNoCorrection(true);
}
