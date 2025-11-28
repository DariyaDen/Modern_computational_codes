// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
// This code processes ROOT histograms for ttbar analysis,
// (produced by ttbarMakeHist.cxx), and makes final plots and numbers
// (more precisely, control plots to be compared to TOP-11-013 Fig. 4, 
// normalised cross sections to be compared to TOP-11-013 Fig. 10 
// and the total cross section to be compared to TOP-13-004).
// Run: ./ttbarMakePlots
// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

// additional files from this analysis 
#include "plots.h"
#include "settings.h"
// C++ library or ROOT header files
#include <TStyle.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TGaxis.h>
#include <TFile.h>
#include <TLegend.h>
#include <TGraphAsymmErrors.h>

// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
// >>>>>>>>>>>>>>>>>>>> Prepare plot style >>>>>>>>>>>>>>>>>>>>>>>>>>>>>
// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
//
// Modify as you want, if needed consult 
// https://root.cern.ch/doc/master/classTStyle.html 
//
void Style()
{
	gStyle->SetOptStat(000000000);
	gStyle->SetTitle(0);
	gStyle->SetFrameFillColor(0);
	gStyle->SetPadColor(0);
	gStyle->SetCanvasColor(0);
	gStyle->SetStatColor(0);
	gStyle->SetCanvasBorderMode(0);
	gStyle->SetCanvasBorderSize(0);
	gStyle->SetFrameBorderMode(0);
	gStyle->SetFrameBorderSize(0);
	gStyle->SetPadTickX(1);
	gStyle->SetPadTickY(1);
	//gStyle->SetPadGridX(1);
	//gStyle->SetPadGridY(1);
	gStyle->SetLegendBorderSize(0);
	gStyle->SetEndErrorSize(5);
    TGaxis::SetMaxDigits(3);
    gStyle->SetErrorX(0.0);
    gStyle->SetPadLeftMargin(0.18);
    gStyle->SetPadBottomMargin(0.12);
    gStyle->SetPadTopMargin(0.06);
    gStyle->SetPadRightMargin(0.08);
    //gStyle->SetNdivisions(206, "xyz");
}

// Additional tuning for plotted histograms (font sizes etc.)
void SetCPHRange(TH2* h)
{
  h->GetXaxis()->SetTitleSize(0.045);
  h->GetXaxis()->SetLabelSize(0.045);
  h->GetYaxis()->SetTitleSize(0.045);
  h->GetYaxis()->SetLabelSize(0.045);
  h->GetXaxis()->SetTitleOffset(1.20);
  h->GetYaxis()->SetTitleOffset(1.70);
}
// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
// >>>>>>>>>>>>>>>>>>>>> Main function >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
int main(int argc, char** argv)
{
  // set user style
  Style();

  // directory with input histograms
  TString baseDir = gHistDir;
  // directory for output plots (must exist)
  TString plotDir = gPlotsDir;
  
  // name patterns for decay channels
  TString suf[3] = {"ee", "mumu", "emu"};
  
  // MC samples for control plots
  std::vector<std::vector<TString> > vecMCName;
  std::vector<int> vecMCColor;
  std::vector<TString> vecMCtitle;

  // Drell-Yan
  std::vector<TString> DYNames;
  DYNames.push_back("DYlm");
  DYNames.push_back("DYhm");
  vecMCName.push_back(DYNames);
  vecMCColor.push_back(kBlue);
  vecMCtitle.push_back("Z / #gamma*");
  // W+jets
  vecMCName.push_back(std::vector<TString>(1, "Wjets"));
  vecMCColor.push_back(kGreen - 2);
  vecMCtitle.push_back("W+Jets");
  // single top
  vecMCName.push_back(std::vector<TString>(1, "SingleTop"));
  vecMCColor.push_back(kMagenta);
  vecMCtitle.push_back("Single Top");
  // ttbar other
  vecMCName.push_back(std::vector<TString>(1, "SigOther"));
  vecMCColor.push_back(kRed - 7);
  vecMCtitle.push_back("t#bar{t} Other");
  // ttbar signal
  vecMCName.push_back(std::vector<TString>(1, "Sig"));
  vecMCColor.push_back(kRed);
  vecMCtitle.push_back("t#bar{t} Signal");

  // *** make control plots ***
  std::vector<TH2F*> cpHR;
  std::vector<TString> cpVar;

  // pT(top)
  TH2F* hr_cp_ptt = new TH2F("hr_cp_ptt", "", 1, 0, 400, 1, 0, 1000.);
  hr_cp_ptt->GetXaxis()->SetTitle("p_{T}(t) [GeV]");
  hr_cp_ptt->GetYaxis()->SetTitle("Top quarks / 20 GeV");
  SetCPHRange(hr_cp_ptt);
  cpHR.push_back(hr_cp_ptt);
  cpVar.push_back("ptt");

  // pT(ttbar)
  TH2F* hr_cp_pttt = new TH2F("hr_cp_pttt", "", 1, 0, 300, 1, 0, 1000.);
  hr_cp_pttt->GetXaxis()->SetTitle("p_{T}(t#bar{t}) [GeV]");
  hr_cp_pttt->GetYaxis()->SetTitle("Top quarks / 20 GeV");
  SetCPHRange(hr_cp_pttt);
  cpHR.push_back(hr_cp_pttt);
  cpVar.push_back("pttt");

  // rapidity(top)
  TH2F* hr_cp_yt = new TH2F("hr_cp_yt", "", 1, -2.6, 2.6, 1, 0, 800.);
  hr_cp_yt->GetXaxis()->SetTitle("y(t)");
  hr_cp_yt->GetYaxis()->SetTitle("Top quarks / 0.2");
  SetCPHRange(hr_cp_yt);
  cpHR.push_back(hr_cp_yt);
  cpVar.push_back("yt");

  // rapidity(ttbar)
  TH2F* hr_cp_ytt = new TH2F("hr_cp_ytt", "", 1, -2.6, 2.6, 1, 0, 1000.);
  hr_cp_ytt->GetXaxis()->SetTitle("y(t#bar{t})");
  hr_cp_ytt->GetYaxis()->SetTitle("Top quarks / 0.2");
  SetCPHRange(hr_cp_ytt);
  cpHR.push_back(hr_cp_ytt);
  cpVar.push_back("ytt");

  // ** NEW: phi(ttbar) **
  TH2F* hr_cp_phi_tt = new TH2F("hr_cp_phi_tt", "", 1, -3.2, 3.2, 1, 0, 1000.);
  hr_cp_phi_tt->GetXaxis()->SetTitle("#phi(t#bar{t})");
  hr_cp_phi_tt->GetYaxis()->SetTitle("Top quarks / 0.2");
  SetCPHRange(hr_cp_phi_tt);
  cpHR.push_back(hr_cp_phi_tt);
  cpVar.push_back("phi_tt");

  // ... the rest of main() remains unchanged, just include "phi_tt" in loops where cpVar is used

  // For cross sections, also add phi_tt
  ZPlotCSInput csIn;
  csIn.Norm = true;
  csIn.Paper = true;
  csIn.baseDir = baseDir;
  csIn.plotDir = plotDir;

  // channels
  csIn.VecColor.push_back(1);
  csIn.VecStyle.push_back(20);
  csIn.VecTitle.push_back("Dilepton");
  csIn.VecColor.push_back(kBlue);
  csIn.VecStyle.push_back(26);
  csIn.VecTitle.push_back("ee");
  csIn.VecColor.push_back(kGreen + 2);
  csIn.VecStyle.push_back(32);
  csIn.VecTitle.push_back("#mu#mu");
  csIn.VecColor.push_back(kRed);
  csIn.VecStyle.push_back(24);
  csIn.VecTitle.push_back("e#mu");

  // MC background
  csIn.VecMCBackgr.push_back("SigOther");
  csIn.VecMCBackgr.push_back("SingleTop");
  csIn.VecMCBackgr.push_back("DYlm");
  csIn.VecMCBackgr.push_back("DYhm");
  csIn.VecMCBackgr.push_back("Wjets");

  // Variables for cross sections
  TH2F* hr_cs_phi_tt = new TH2F("hr_cs_phi_tt", "", 1, -3.2, 3.2, 1, 0, 0.8);
  hr_cs_phi_tt->GetXaxis()->SetTitle("#phi(t#bar{t})");
  hr_cs_phi_tt->GetYaxis()->SetTitle("#frac{1}{#sigma} #frac{d#sigma}{d#phi(t#bar{t})}");
  SetCPHRange(hr_cs_phi_tt);
  csIn.VecHR.push_back(hr_cs_phi_tt);
  csIn.VecVar.push_back("phi_tt");

  // ... rest of main() remains unchanged: PlotCS(csIn), return 0, etc.
  
  return 0;
}
