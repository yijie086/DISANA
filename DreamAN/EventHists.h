#define EVENTHIST_H_
#include "Rtypes.h"
#include "TList.h"
#include "TH1F.h"
#include "TH2F.h"
class EventHist {
 public:
  EventHist();
  EventHist(bool centVsMultPlot);
  EventHist(const EventHist& hists);
  EventHist& operator=(const EventHist& hists);
  virtual ~EventHist();
  void FillEvtCounter(int iBin) {
    fEvtCounter->Fill(iBin);
  }
  ;
  void FillEvtVtxX(int i, float val) {
    fEvtVtxX[i]->Fill(val);
  }
  ;
  void FillEvtVtxY(int i, float val) {
    fEvtVtxY[i]->Fill(val);
  }
  ;
  void FillEvtVtxZ(int i, float val) {
    fEvtVtxZ[i]->Fill(val);
  }
;

  TList *GetHistList() {
    return fEventCutList;
  }
  ;
 private:
  TH1F *fEvtCounter;        //!
  TH1F *fEvtVtxX[2];        //!
  TH1F *fEvtVtxY[2];        //!
  TH1F *fEvtVtxZ[2];        //!
};

#endif /* EVENTHIST_H_ */