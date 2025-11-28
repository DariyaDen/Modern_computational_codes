#ifndef TTBAR_EVENTRECO_H
#define TTBAR_EVENTRECO_H

#include "tree.h"
#include "kinReco.h"
#include "selection.h"
#include "settings.h"

#include <map>
#include <vector>
#include <TChain.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TH1D.h>
#include <TLorentzVector.h>
#include <TMath.h>

class ZVarHisto
{
private:
    TH1D* zHisto;
    TString zVar;

public:
    ZVarHisto(const TString& str, TH1D* h)
    {
        zHisto = h;
        zVar = str;
    }

    ZVarHisto(const ZVarHisto& old)
    {
        zHisto = new TH1D(*(old.zHisto));
        zVar = old.zVar;
    }

    TH1* H() { return zHisto; }
    TString V() { return zVar; }
};

// ----------------- FillHistos with new variable Δφ -----------------
void FillHistos(std::vector<ZVarHisto>& VecVarHisto, double w,
                TLorentzVector* t, TLorentzVector* tbar,
                TLorentzVector* vecLepM = NULL, TLorentzVector* vecLepP = NULL)
{
    TLorentzVector ttbar = *t + *tbar;

    for (int h = 0; h < VecVarHisto.size(); h++)
    {
        TString var = VecVarHisto[h].V();
        TH1* histo = VecVarHisto[h].H();

        if (var == "ptt") histo->Fill(t->Pt(), w);
        else if (var == "ptat") histo->Fill(tbar->Pt(), w);
        else if (var == "pttat") { histo->Fill(t->Pt(), w); histo->Fill(tbar->Pt(), w); }
        else if (var == "pttt") histo->Fill(ttbar.Pt(), w);
        else if (var == "yt") histo->Fill(t->Rapidity(), w);
        else if (var == "yat") histo->Fill(tbar->Rapidity(), w);
        else if (var == "ytat") { histo->Fill(t->Rapidity(), w); histo->Fill(tbar->Rapidity(), w); }
        else if (var == "ytt") histo->Fill(ttbar.Rapidity(), w);
        else if (var == "mtt") histo->Fill(ttbar.M(), w);
        else if (var == "ptl") {
            histo->Fill(vecLepM->Pt(), w);
            histo->Fill(vecLepP->Pt(), w);
        }
        // ---------------- NEW VARIABLE: Δφ between top and antitop ----------------
        else if (var == "dphi_tt") {
            double dphi = t->DeltaPhi(*tbar); // Δφ(t, tbar)
            histo->Fill(dphi, w);
        }
        else continue;
    }
}

// ----------------- StoreHistos -----------------
void StoreHistos(std::vector<ZVarHisto>& VecVarHisto)
{
    for (int h = 0; h < VecVarHisto.size(); h++)
    {
        TH1* histo = VecVarHisto[h].H();
        TString name = histo->GetName();
        TString title = histo->GetTitle();
        histo->SetNameTitle(name, title);
        histo->Write();
    }
}

// ----------------- Input parameters -----------------
class ZEventRecoInput
{
public:
    TString Name;
    std::vector<ZVarHisto> VecVarHisto;
    int Channel;
    int Type;
    bool Gen;
    std::vector<TString> VecInFile;
    double Weight;
    long MaxNEvents;

    ZEventRecoInput() {
        Weight = 1.0;
        MaxNEvents = 100e10;
        Gen = false;
    }

    void AddToChain(const TString& str) { VecInFile.push_back(str); }
    void ClearChain() { VecInFile.clear(); }
};

// ----------------- Main reconstruction -----------------
void eventreco(ZEventRecoInput in)
{
    printf("****** EVENTRECO ******\n");
    printf("input sample: %s\n", in.Name.Data());
    printf("type: %d   channel: %d\n", in.Type, in.Channel);

    const double bTagDiscrL = 0.244;
    TString outDir = gHistDir;
    bool flagMC = (in.Type == 2 || in.Type == 3);
    
    TFile* fout = TFile::Open(TString::Format("%s/%s-c%d.root", outDir.Data(), in.Name.Data(), in.Channel), "recreate");
    TChain* chain = new TChain("tree");
    for (int f = 0; f < in.VecInFile.size(); f++) chain->Add(in.VecInFile[f]);
    
    ZTree* preselTree = new ZTree(flagMC);
    preselTree->Init(chain);

    if (in.Gen)
    {
        chain->SetBranchStatus("*", 0);
        chain->SetBranchStatus("mcEventType", 1);
        chain->SetBranchStatus("mcT", 1);
        chain->SetBranchStatus("mcTbar", 1);
    }

    long nSel = 0, nReco = 0;
    long nEvents = chain->GetEntries();
    if (nEvents > in.MaxNEvents) nEvents = in.MaxNEvents;
    printf("nEvents: %ld\n", nEvents);

    TH1D* hInacc = new TH1D("hInacc", "KinReco inaccuracy", 1000, 0.0, 100.0);
    TH1D* hAmbig = new TH1D("hAmbig", "KinReco ambiguity", 100, 0.0, 100.0);

    for (int e = 0; e < nEvents; e++)
    {
        chain->GetEntry(e);

        if (flagMC)
        {
            if (in.Type == 2 && preselTree->mcEventType != in.Channel) continue;
            if (in.Type == 3 && preselTree->mcEventType == in.Channel) continue;
        }

        if (in.Gen)
        {
            TLorentzVector t, tbar;
            t.SetXYZM(preselTree->mcT[0], preselTree->mcT[1], preselTree->mcT[2], preselTree->mcT[3]);
            tbar.SetXYZM(preselTree->mcTbar[0], preselTree->mcTbar[1], preselTree->mcTbar[2], preselTree->mcTbar[3]);
            FillHistos(in.VecVarHisto, in.Weight, &t, &tbar);
            continue;
        }

        if (preselTree->Npv < 1 || preselTree->pvNDOF < 4 || preselTree->pvRho > 2.0 || TMath::Abs(preselTree->pvZ) > 24.0)
            continue;

        TLorentzVector vecLepM, vecLepP;
        double maxPtDiLep = -1.0;
        bool trig = false;

        if (in.Channel == 3)
        {
            for (int bit = 12; bit < 17; bit++)
                if ((preselTree->Triggers >> bit) & 1) { trig = true; break; }
            if (trig) SelectDilepEMu(preselTree, vecLepM, vecLepP, maxPtDiLep);
        }
        if (in.Channel == 1)
        {
            for (int bit = 6; bit < 11; bit++)
                if ((preselTree->Triggers >> bit) & 1) { trig = true; break; }
            double met = TMath::Sqrt(preselTree->metPx*preselTree->metPx + preselTree->metPy*preselTree->metPy);
            if (trig && met > 30.0) SelectDilepEE(preselTree, vecLepM, vecLepP, maxPtDiLep);
        }
        if (in.Channel == 2)
        {
            for (int bit = 0; bit < 5; bit++)
                if ((preselTree->Triggers >> bit) & 1) { trig = true; break; }
            double met = TMath::Sqrt(preselTree->metPx*preselTree->metPx + preselTree->metPy*preselTree->metPy);
            if (trig && met > 30.0) SelectDilepMuMu(preselTree, vecLepM, vecLepP, maxPtDiLep);
        }

        if (maxPtDiLep < 0.0) continue;

        std::vector<TLorentzVector> vecJets;
        bool oneBTagJet = false;
        for (int j = 0; j < preselTree->Njet; j++)
        {
            if (TMath::Abs(preselTree->jetEta[j]) > 2.4) continue;
            TLorentzVector vecJet;
            vecJet.SetPtEtaPhiM(preselTree->jetPt[j], preselTree->jetEta[j], preselTree->jetPhi[j], preselTree->jetMass[j]);
            double corrE = vecJet.E() - preselTree->jetMuEn[j] - preselTree->jetElEn[j];
            double corrPt = preselTree->jetPt[j] * corrE / vecJet.E();
            if (corrPt < 30.0) continue;
            TLorentzVector corrVec;
            corrVec.SetPtEtaPhiE(corrPt, preselTree->jetEta[j], preselTree->jetPhi[j], corrE);
            if (preselTree->jetBTagDiscr[j] > bTagDiscrL)
            {
                corrVec.SetPtEtaPhiM(corrVec.Pt(), corrVec.Eta(), corrVec.Phi(), -1*corrVec.M());
                oneBTagJet = true;
            }
            vecJets.push_back(corrVec);
        }

        if (vecJets.size() < 2 || !oneBTagJet) continue;
        nSel++;

        TLorentzVector t, tbar;
        int status = KinRecoDilepton(vecLepM, vecLepP, vecJets, preselTree->metPx, preselTree->metPy, t, tbar, hInacc, hAmbig);

        if (status > 0)
        {
            nReco++;
            FillHistos(in.VecVarHisto, in.Weight, &t, &tbar, &vecLepM, &vecLepP);
        }
    }

    printf("nSel  : %ld\n", nSel);
    printf("nReco : %ld\n", nReco);

    fout->cd();
    StoreHistos(in.VecVarHisto);
    fout->Close();
}

#endif
