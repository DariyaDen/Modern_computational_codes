// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
// >>>>>>>>>>>> Helper for ttbar event selection >>>>>>>>>>>>>>>>>>>>>>>
// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
// Consult analysis documentation (papers, description-ttbar.pdf) for 
// better description of applied cuts etc.

// additional files from this analysis 
#include "tree.h"
// C++ library or ROOT header files
#include <TLorentzVector.h>
#include <TMath.h>

// constants: electron and muon masses
const double massEl = 0.000511;
const double massMu = 0.105658;

// Routine for electron selection
bool SelectEl(const ZTree* preselTree, const int el)
{
    if(TMath::Abs(preselTree->elPt[el]) < 20.0) return false;
    if(TMath::Abs(preselTree->elEta[el]) > 2.4) return false;
    if(preselTree->elIso03[el] > 0.17) return false;
    if(preselTree->elMissHits[el] > 0) return false;
    return true;
}

// Routine for muon selection
bool SelectMu(const ZTree* preselTree, const int mu)
{
    if(TMath::Abs(preselTree->muPt[mu]) < 20.0) return false;
    if(TMath::Abs(preselTree->muEta[mu]) > 2.4) return false;
    if(preselTree->muIso03[mu] > 0.20) return false;
    if(preselTree->muHitsValid[mu] < 12 || preselTree->muHitsPixel[mu] < 2) return false;
    if(preselTree->muDistPV0[mu] > 0.02 || preselTree->muDistPVz[mu] > 0.5 || preselTree->muTrackChi2NDOF[mu] > 10) return false;
    return true;
}

// Generic routine for dilepton selection (e-mu, ee, or mu-mu)
void SelectDilepPair(const ZTree* preselTree,
                     TLorentzVector& vecLepM,
                     TLorentzVector& vecLepP,
                     double& maxPtDiLep,
                     int lepType = 0) // 0 = e-mu, 1 = ee, 2 = mu-mu
{
    for(int i = 0; i < ((lepType==2)?preselTree->Nmu:preselTree->Nel); i++)
    {
        bool pass1 = (lepType==2) ? SelectMu(preselTree,i) : SelectEl(preselTree,i);
        if(!pass1) continue;

        TLorentzVector lep1;
        if(lepType==2)
            lep1.SetPtEtaPhiM(TMath::Abs(preselTree->muPt[i]), preselTree->muEta[i], preselTree->muPhi[i], massMu);
        else
            lep1.SetPtEtaPhiM(TMath::Abs(preselTree->elPt[i]), preselTree->elEta[i], preselTree->elPhi[i], massEl);

        for(int j = (lepType==1 || lepType==2)?i+1:0; j < ((lepType==1)?preselTree->Nel:preselTree->Nmu); j++)
        {
            bool pass2 = (lepType==2) ? SelectMu(preselTree,j) : SelectEl(preselTree,j);
            if(!pass2) continue;

            TLorentzVector lep2;
            if(lepType==2)
                lep2.SetPtEtaPhiM(TMath::Abs(preselTree->muPt[j]), preselTree->muEta[j], preselTree->muPhi[j], massMu);
            else if(lepType==1)
                lep2.SetPtEtaPhiM(TMath::Abs(preselTree->elPt[j]), preselTree->elEta[j], preselTree->elPhi[j], massEl);
            else
                lep2.SetPtEtaPhiM(TMath::Abs(preselTree->muPt[j]), preselTree->muEta[j], preselTree->muPhi[j], massMu);

            // require opposite signs
            if((lepType==0 && preselTree->elPt[i]*preselTree->muPt[j]>0) ||
               (lepType==1 && preselTree->elPt[i]*preselTree->elPt[j]>0) ||
               (lepType==2 && preselTree->muPt[i]*preselTree->muPt[j]>0)) continue;

            TLorentzVector vecDiLep = lep1 + lep2;
            if(vecDiLep.M() < 12.0) continue;
            if(lepType==1 || lepType==2)
            {
                if(vecDiLep.M() > 76.0 && vecDiLep.M() < 106.0) continue; // DY veto
            }

            double sumPt = lep1.Pt() + lep2.Pt();
            if(sumPt < maxPtDiLep) continue;
            maxPtDiLep = sumPt;
            vecLepM = (lep1.Pt()<0)?lep1:lep2;
            vecLepP = (lep1.Pt()<0)?lep2:lep1;

            // === NEW KINEMATIC VARIABLE ===
            // Δφ between two leptons (radians)
            double deltaPhiLL = TVector2::Phi_mpi_pi(vecLepM.Phi() - vecLepP.Phi());
            // save or fill histogram for Δφ_ll here
            // e.g. hDeltaPhiLL->Fill(deltaPhiLL);
        }
    }
}
