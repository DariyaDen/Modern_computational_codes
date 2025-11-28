// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
// >>>>>>>>>>>>>>>>>>>> Helper for plotter >>>>>>>>>>>>>>>>>>>>>>>>>>>>>
// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

// additional files from this analysis 
#include "settings.h"
// C++ library or ROOT header files
#include <TCanvas.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TGaxis.h>
#include <TFile.h>
#include <TLegend.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
#include <TMath.h>
#include <TVector2.h>

// produce graph from histogram
TGraphAsymmErrors* HtoGraph(const TH1* h, double xpos = 0.5)
{
    TGraphAsymmErrors* g = new TGraphAsymmErrors;
    for(int b = 0; b < h->GetNbinsX(); b++)
    {
        double xmin = h->GetXaxis()->GetBinLowEdge(b + 1);
        double xmax = h->GetXaxis()->GetBinUpEdge(b + 1);
        double bw = xmax - xmin;
        double x = xmin + xpos * bw;
        double y = h->GetBinContent(b + 1);
        double yer = h->GetBinError(b + 1);
        int point = g->GetN();
        g->SetPoint(point, x, y);
        g->SetPointError(point, 0.0, 0.0, yer, yer);
    }
    return g;
}

// Helper class for cross-section plotting
class ZPlotCSInput
{
public:
    TString baseDir;
    TString plotDir;
    std::vector<TString> VecMCBackgr;
    std::vector<TString> VecVar;
    std::vector<TH2F*> VecHR;
    std::vector<TString> VecTitle;
    std::vector<int> VecColor;
    std::vector<int> VecStyle;
    bool Norm;
    bool Paper;
};

// === NEW FUNCTION: Fill histogram for ΔR between leptons ===
void FillDeltaR(TH1D* hDeltaR, const TLorentzVector& lep1, const TLorentzVector& lep2)
{
    double deltaEta = lep1.Eta() - lep2.Eta();
    double deltaPhi = TVector2::Phi_mpi_pi(lep1.Phi() - lep2.Phi());
    double deltaR = sqrt(deltaEta*deltaEta + deltaPhi*deltaPhi);
    hDeltaR->Fill(deltaR);
}

// routine to calculate and plot x-sections
void PlotCS(const ZPlotCSInput& in)
{
    TCanvas* c_cs = new TCanvas("ccs", "", 1200, 800);
    c_cs->Divide(3, 2, 0.0001);

    // loop over variables
    for(int v = 0; v < in.VecVar.size(); v++)
    {
        c_cs->cd(v + 1);
        if(v == 4) gPad->SetLogy();
        in.VecHR[v]->Draw();

        TString var = in.VecVar[v];
        TH1D* hcombsig = nullptr;
        TH1D* hcombreco = nullptr;
        TH1D* hcombgen = nullptr;
        TLegend* leg = new TLegend(0.46, 0.68, 0.92, 0.90);
        leg->SetTextSize(0.045);
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);

        for(int ch = 1; ch < 4; ch++)
        {
            // MC background
            TH1D* hbackgr = nullptr;
            for(int mc = 0; mc < in.VecMCBackgr.size(); mc++)
            {
                TFile* f = TFile::Open(TString::Format("%s/mc%sReco-c%d.root", in.baseDir.Data(), in.VecMCBackgr[mc].Data(), ch));
                TH1D* h = (TH1D*)f->Get(TString::Format("h_%s_cs", var.Data()));
                if(mc == 0)
                    hbackgr = new TH1D(*h);
                else
                    hbackgr->Add(h);
            }

            // data
            TFile* fdata = TFile::Open(TString::Format("%s/data-c%d.root", in.baseDir.Data(), ch));
            TH1D* hsig = new TH1D(*(TH1D*)fdata->Get(TString::Format("h_%s_cs", var.Data())));
            hsig->Add(hbackgr, -1.0);

            if(ch == 1) hcombsig = new TH1D(*hsig);
            else hcombsig->Add(hsig);

            // MC reconstruction level
            fdata = TFile::Open(TString::Format("%s/mcSigReco-c%d.root", in.baseDir.Data(), ch));
            TH1D* hacc = new TH1D(*(TH1D*)fdata->Get(TString::Format("h_%s_cs", var.Data())));
            if(ch == 1) hcombreco = new TH1D(*hacc);
            else hcombreco->Add(hacc);

            // MC generator level
            fdata = TFile::Open(TString::Format("%s/mcSigGen-c%d.root", in.baseDir.Data(), ch));
            TH1D* hgen = (TH1D*)fdata->Get(TString::Format("h_%s_cs", var.Data()));
            if(ch == 1) hcombgen = new TH1D(*hgen);
            else hcombgen->Add(hgen);

            // ACCEPTANCE
            hacc->Divide(hgen);
            hsig->Divide(hacc);

            if(in.Norm) hsig->Scale(1.0 / hsig->Integral(), "width");

            // convert histogram to TGraph
            TGraphAsymmErrors* gcs = HtoGraph(hsig, 0.10 + 0.1*ch);
            gcs->SetLineColor(in.VecColor[ch]);
            gcs->SetMarkerStyle(in.VecStyle[ch]);
            gcs->SetMarkerColor(in.VecColor[ch]);
            gcs->SetMarkerSize(0.9);
            leg->AddEntry(gcs, in.VecTitle[ch], "p");
            gcs->Draw("pz0");

            // === Example: Fill ΔR histogram for control plots ===
            // TH1D* hDeltaR = new TH1D("hDeltaR", "DeltaR between leptons;#DeltaR;Events", 50, 0, 5);
            // for each event: FillDeltaR(hDeltaR, lep1, lep2);
        }

        leg->Draw();
    }

    TString name = TString::Format("%s/cs%s", in.plotDir.Data(), (in.Norm) ? "_norm" : "");
    c_cs->SaveAs(name + ".eps");
    c_cs->SaveAs(name + ".pdf");
}
