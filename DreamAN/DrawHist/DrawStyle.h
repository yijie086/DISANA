#ifndef DRAWSTYLE_H
#define DRAWSTYLE_H

#include <TPad.h>
#include <TH1.h>
#include <TH2.h>
#include <TStyle.h>
#include <TGaxis.h>

class DrawStyle {
 public:
  DrawStyle(double titleSize = 0.05, double labelSize = 0.04,
            double xTitleOffset = 1.1, double yTitleOffset = 1.6,
            int font = 42, int maxDigits = 5, int nDivisions = 510,
            double leftMargin = 0.16, double rightMargin = 0.07,
            double bottomMargin = 0.13, double topMargin = 0.06)
      : titleSize_(titleSize),
        labelSize_(labelSize),
        xTitleOffset_(xTitleOffset),
        yTitleOffset_(yTitleOffset),
        font_(font),
        maxDigits_(maxDigits),
        nDivisions_(nDivisions),
        leftMargin_(leftMargin),
        rightMargin_(rightMargin),
        bottomMargin_(bottomMargin),
        topMargin_(topMargin) {}

  void ApplyGlobalStyle() const {
    gStyle->SetTitleSize(titleSize_, "XYZ");
    gStyle->SetLabelSize(labelSize_, "XYZ");
    gStyle->SetTitleFont(font_, "XYZ");
    gStyle->SetLabelFont(font_, "XYZ");

    gStyle->SetTitleOffset(xTitleOffset_, "X");
    gStyle->SetTitleOffset(yTitleOffset_, "Y");

    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
    gStyle->SetOptStat(0);
    gStyle->SetCanvasColor(0);
    gStyle->SetFrameBorderMode(0);
    gStyle->SetLegendBorderSize(0);
    gStyle->SetPalette(kRainBow);

    TGaxis::SetMaxDigits(maxDigits_);
  }

  void NormalizeHistogram(TH1* hist) const {
  if (!hist) return;
  double integral = hist->Integral();
  if (integral > 0)
    hist->Scale(1.0 / integral);
  }

  void StylePad(TPad* pad) const {
    if (!pad) return;
    pad->SetLeftMargin(leftMargin_);
    pad->SetRightMargin(rightMargin_);
    pad->SetBottomMargin(bottomMargin_);
    pad->SetTopMargin(topMargin_);
    pad->SetTicks(1, 1);
  }

  void StyleTH1(TH1* hist, const char* title = "") const {
    if (!hist) return;

    hist->SetStats(0);
    hist->SetLineWidth(2);
    hist->SetTitle(title);
    hist->SetTitleFont(font_, "");
    hist->SetTitleSize(titleSize_, "");

    hist->GetXaxis()->SetTitleSize(titleSize_);
    hist->GetXaxis()->SetLabelSize(labelSize_);
    hist->GetXaxis()->SetTitleOffset(xTitleOffset_);
    hist->GetXaxis()->SetTitleFont(font_);
    hist->GetXaxis()->SetLabelFont(font_);
    hist->GetXaxis()->SetNdivisions(nDivisions_);

    hist->GetYaxis()->SetTitleSize(titleSize_);
    hist->GetYaxis()->SetLabelSize(labelSize_);
    hist->GetYaxis()->SetTitleOffset(yTitleOffset_);
    hist->GetYaxis()->SetTitleFont(font_);
    hist->GetYaxis()->SetLabelFont(font_);
    hist->GetYaxis()->SetNdivisions(nDivisions_);
  }

  void StyleTH2(TH2* hist, const char* title = "", double zLabelSize = -1.0) const {
    if (!hist) return;

    hist->SetStats(0);
    hist->SetTitle(title);
    hist->SetTitleFont(font_, "");
    hist->SetTitleSize(titleSize_, "");

    hist->GetXaxis()->SetTitleSize(titleSize_);
    hist->GetXaxis()->SetLabelSize(labelSize_);
    hist->GetXaxis()->SetTitleOffset(xTitleOffset_);
    hist->GetXaxis()->SetTitleFont(font_);
    hist->GetXaxis()->SetLabelFont(font_);
    hist->GetXaxis()->SetNdivisions(nDivisions_);

    hist->GetYaxis()->SetTitleSize(titleSize_);
    hist->GetYaxis()->SetLabelSize(labelSize_);
    hist->GetYaxis()->SetTitleOffset(yTitleOffset_);
    hist->GetYaxis()->SetTitleFont(font_);
    hist->GetYaxis()->SetLabelFont(font_);
    hist->GetYaxis()->SetNdivisions(nDivisions_);

    hist->GetZaxis()->SetLabelSize(zLabelSize > 0 ? zLabelSize : labelSize_);
    hist->GetZaxis()->SetTitleFont(font_);
    hist->GetZaxis()->SetNdivisions(nDivisions_);
  }

 private:
  double titleSize_;
  double labelSize_;
  double xTitleOffset_;
  double yTitleOffset_;
  int font_;
  int maxDigits_;
  int nDivisions_;
  double leftMargin_, rightMargin_, bottomMargin_, topMargin_;
};

#endif  // DRAWSTYLE_H
