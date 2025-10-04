#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TRandom3.h>
#include <Math/Vector4D.h>
#include <Math/Boost.h>
#include <vector>
#include <cmath>
#include <iostream>

using namespace ROOT::Math;
using namespace std;

// Smear particle momentum
PxPyPzMVector smear_momentum(const PxPyPzMVector &particle, TRandom3 &randGen) {
    double fluct = randGen.Gaus(1., 0.05);
    return PxPyPzMVector(particle.Px()*fluct, particle.Py()*fluct, particle.Pz()*fluct, particle.M());
}

// Generate two-body decay in parent rest frame
pair<PxPyPzMVector, PxPyPzMVector> two_body_decay(const PxPyPzMVector &parent, double m1, double m2, TRandom3 &randGen) {
    double m0 = parent.M();
    double p12 = sqrt((m0*m0 - pow(m1+m2,2))*(m0*m0 - pow(m1-m2,2)))/(2.*m0);
    double E1 = (m0*m0 + m1*m1 - m2*m2)/(2.*m0);
    double p = sqrt(E1*E1 - m1*m1);

    double theta = randGen.Uniform(0, M_PI);
    double phi = randGen.Uniform(0, 2*M_PI);

    double p1x = p*sin(theta)*cos(phi);
    double p1y = p*sin(theta)*sin(phi);
    double p1z = p*cos(theta);

    PxPyPzMVector daughter1(p1x, p1y, p1z, m1);
    PxPyPzMVector daughter2(-p1x, -p1y, -p1z, m2);

    Boost boost(parent.BoostToCM());
    daughter1 = boost(daughter1);
    daughter2 = boost(daughter2);

    return make_pair(daughter1, daughter2);
}

// Generate parent particle and decay into two daughters
vector<PxPyPzMVector> generate_and_decay(double parent_mass, double m1, double m2, TRandom3 &randGen) {
    double momentum = randGen.Exp(1.0);
    double theta = randGen.Uniform(0, M_PI);
    double phi = randGen.Uniform(0, 2*M_PI);

    double px = momentum*sin(theta)*cos(phi);
    double py = momentum*sin(theta)*sin(phi);
    double pz = momentum*cos(theta);

    PxPyPzMVector parent(px, py, pz, parent_mass);

    auto daughters = two_body_decay(parent, m1, m2, randGen);

    vector<PxPyPzMVector> tracks;
    tracks.push_back(smear_momentum(daughters.first, randGen));
    tracks.push_back(smear_momentum(daughters.second, randGen));

    return tracks;
}

int invmass() {
    gROOT->SetBatch(kTRUE);  // No GUI windows
    TRandom3 randGen(42);

    int nParticles = 1000;

    double mass_pi_ch = 0.13957;
    double mass_k_zero = 0.497611;
    double mass_d_zero = 1.86484;

    TH1F *hInvMass = new TH1F("hInvMass", "Invariant Mass", 300, 0, 3);

    TFile *fileout = new TFile("tracks.root", "RECREATE");
    TTree *tree = new TTree("tree", "Tree with tracks");
    vector<PxPyPzMVector> tracks_vec;
    tree->Branch("tracks", &tracks_vec);

    for(int i=0; i<nParticles; i++) {
        tracks_vec.clear();
        auto tracks1 = generate_and_decay(mass_k_zero, mass_pi_ch, mass_pi_ch, randGen);
        auto tracks2 = generate_and_decay(mass_d_zero, mass_pi_ch, mass_pi_ch, randGen);

        tracks_vec.insert(tracks_vec.end(), tracks1.begin(), tracks1.end());
        tracks_vec.insert(tracks_vec.end(), tracks2.begin(), tracks2.end());

        // Fill histogram
        for(size_t j=0; j<tracks_vec.size(); j++) {
            for(size_t k=j+1; k<tracks_vec.size(); k++) {
                hInvMass->Fill((tracks_vec[j]+tracks_vec[k]).M());
            }
        }

        tree->Fill();
    }

    fileout->cd();
    tree->Write();
    fileout->Close();

    // Fit function: two Gaussians + 5th-degree polynomial
    TF1 *fitFunc = new TF1("fitFunc",
                           "[0]/(sqrt(2*TMath::Pi())*[2])*TMath::Exp(-0.5*((x-[1])/[2])^2)"
                           "+ [3]/(sqrt(2*TMath::Pi())*[5])*TMath::Exp(-0.5*((x-[4])/[5])^2)"
                           "+ [6] + [7]*x + [8]*x^2 + [9]*x^3 + [10]*x^4 + [11]*x^5",
                           -10, 10);

    Double_t params[12] = {1, 0.5, 0.01, 1, 1.85, 0.05, 0, 0, 0, 0, 0, 0};
    fitFunc->SetParameters(params);

    fitFunc->SetParLimits(1, 0.48, 0.52); // K0 mass
    fitFunc->SetParLimits(2, 0.0, 0.1);   // K0 width
    fitFunc->SetParLimits(4, 1.8, 1.9);   // D0 mass
    fitFunc->SetParLimits(5, 0.0, 0.2);   // D0 width

    hInvMass->Fit(fitFunc);

    Double_t bin_width = hInvMass->GetBinLowEdge(2) - hInvMass->GetBinLowEdge(1);
    cout << "Number of signal events #1 = " << fitFunc->GetParameter(0)/bin_width
         << " +- " << fitFunc->GetParError(0)/bin_width << endl;
    cout << "Number of signal events #2 = " << fitFunc->GetParameter(3)/bin_width
         << " +- " << fitFunc->GetParError(3)/bin_width << endl;

    TCanvas *canvas = new TCanvas("canvas", "Invariant Mass", 600, 600);
    hInvMass->GetXaxis()->SetTitle("M(#pi^{+}#pi^{-}) [GeV]");
    hInvMass->GetYaxis()->SetTitle("Events");
    hInvMass->Draw();

    canvas->SaveAs("invmass.pdf");
    canvas->SaveAs("invmass.png");

    return 0;
}

