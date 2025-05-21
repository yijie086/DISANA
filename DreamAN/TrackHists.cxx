#include "TrackHists.h"
#include "TMath.h"

//ClassImp(TrackHists)

TrackHists::TrackHists()
    : fMinimalBooking(false),
      fpmin(0),
      fpmax(0),
      fHistList(0),
      fConfig(0) {
  for (int i = 0; i < 2; ++i) {
    fTrackCutQA[i] = nullptr;
    fpDist[i] = nullptr;
    fphiDist[i] = nullptr;
    fTrackChi2[i] = nullptr;
    fVz[i] = nullptr;
  }
}

TrackHists::TrackHists(bool VzDist, float pmin , float pmax )
    : fMinimalBooking(false),
      fpmin(pmin),
      fpmax(pmax) {

  TString sName[2] = { "before", "after" };
  int pBins = 120;
  int twoDBins = 400;
  fHistList = new TList();
  fHistList->SetName("TrackCuts");
  fHistList->SetOwner();

  fConfig = new TProfile("TrackCutConfig", "Track Cut Config", 30, 0, 30);
  fConfig->SetStats(0);
  fConfig->GetXaxis()->SetBinLabel(1, "p_{min}");
  fConfig->GetXaxis()->SetBinLabel(2, "p_{max}");
  fConfig->GetXaxis()->SetBinLabel(3, "Charge");
  fConfig->GetXaxis()->SetBinLabel(4, "DCA_{xy}");
  fConfig->GetXaxis()->SetBinLabel(5, "DCA_{z}");
  fConfig->GetXaxis()->SetBinLabel(6, "#chi^{2} min");
  fConfig->GetXaxis()->SetBinLabel(7, "#chi^{2} max");
  fHistList->Add(fConfig);

  for (int i = 0; i < 2; ++i) {
    fTrackCutQA[i] = new TList();
    fTrackCutQA[i]->SetName(sName[i].Data());
    fTrackCutQA[i]->SetOwner();
    fHistList->Add(fTrackCutQA[i]);

    TString pITSName = Form("pDist_%s", sName[i].Data());
    fpDist[i] = new TH1F(pITSName.Data(), pITSName.Data(), 1000 * pmax, 0,
                            pmax);
    fpDist[i]->GetXaxis()->SetTitle("#it{p}_{} (GeV/#it{c})");
    fpDist[i]->GetYaxis()->SetTitle("Entries");
    fTrackCutQA[i]->Add(fpDist[i]);


    TString phiName = Form("phiDist_%s", sName[i].Data());
    fphiDist[i] = new TH1F(phiName.Data(), phiName.Data(), 200, 0.,
                           2 * TMath::Pi());
    fphiDist[i]->GetXaxis()->SetTitle("#phi");
    fphiDist[i]->GetYaxis()->SetTitle("Entries");
    fTrackCutQA[i]->Add(fphiDist[i]);

    TString VzName = Form("Vz_%s", sName[i].Data());
    fVz[i] = new TH2F(VzName.Data(), VzName.Data(), pBins, pmin, pmax,
                        2.5 * twoDBins, -5., 5.);
    fVz[i]->GetXaxis()->SetTitle("#it{p} (GeV/#it{c})");
    fVz[i]->GetYaxis()->SetTitle("V_{z}");
    fTrackCutQA[i]->Add(fVz[i]);

    TString ChiSquareName = Form("TrackChi2_%s", sName[i].Data());
    fTrackChi2[i] = new TH2F(ChiSquareName.Data(), ChiSquareName.Data(), pBins,
                             pmin, pmax, 100, 0, 20);
    fTrackChi2[i]->GetXaxis()->SetTitle("#it{p} (GeV/#it{c})");
    fTrackChi2[i]->GetYaxis()->SetTitle("#chi^{2}");
    fTrackCutQA[i]->Add(fTrackChi2[i]);
  }
}

// for memory mngmnt its necessary to delete all the new members created
TrackHists::TrackHists(TString MinimalBooking)
    : fMinimalBooking(true),
      fConfig(0)
      {
  for (int i = 0; i < 2; ++i) {
    fTrackCutQA[i] = nullptr;
    fpDist[i] = nullptr;
    fphiDist[i] = nullptr;
    fVz[i] = nullptr;
  }
  fHistList = new TList();
  fHistList->SetOwner();
  fHistList->SetName(MinimalBooking.Data());

  fpDist[0] = 0;
  TString ptName = Form("pDist_%s", "after");
  fpDist[1] = new TH1F(ptName.Data(), ptName.Data(), 200, 0, 5);
  fpDist[1]->GetXaxis()->SetTitle("#it{p} (GeV/#it{c})");
  fpDist[1]->GetYaxis()->SetTitle("Entries");
  fHistList->Add(fpDist[1]);
}

TrackHists::~TrackHists() {
  delete fHistList;
}

// we can include more histograms here
