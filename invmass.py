import ROOT
import math
import random
import tqdm
import sys
from array import array

def generate_and_decay(parent_mass, daughter_mass_1, daughter_mass_2):
    # Generate random values for momentum and direction
    momentum = random.expovariate(1.0)
    theta = random.uniform(0, math.pi)
    phi = random.uniform(0, 2 * math.pi)

    px = momentum * math.sin(theta) * math.cos(phi)
    py = momentum * math.sin(theta) * math.sin(phi)
    pz = momentum * math.cos(theta)

    parent = ROOT.Math.PxPyPzMVector(px, py, pz, parent_mass)
    daughter1, daughter2 = two_body_decay(parent, daughter_mass_1, daughter_mass_2)

    # Smear momentum of daughter tracks
    daughter1 = smear_momentum(daughter1)
    daughter2 = smear_momentum(daughter2)

    return [daughter1, daughter2]

def smear_momentum(particle):
    fluct = random.normalvariate(1., 0.05)
    return ROOT.Math.PxPyPzMVector(particle.X() * fluct,
                                   particle.Y() * fluct,
                                   particle.Z() * fluct,
                                   particle.M())

def two_body_decay(parent, m1, m2):
    m0 = parent.M()
    p12 = math.sqrt(max(0.0, (m0**2 - (m1 + m2)**2) * (m0**2 - (m1 - m2)**2))) / (2 * m0)
    E1 = (m0**2 + m1**2 - m2**2) / (2 * m0)
    p = math.sqrt(max(0.0, E1**2 - m1**2))

    theta = random.uniform(0, math.pi)
    phi = random.uniform(0, 2 * math.pi)

    p1x = p * math.sin(theta) * math.cos(phi)
    p1y = p * math.sin(theta) * math.sin(phi)
    p1z = p * math.cos(theta)

    daughter1 = ROOT.Math.PxPyPzMVector(p1x, p1y, p1z, m1)
    daughter2 = ROOT.Math.PxPyPzMVector(-p1x, -p1y, -p1z, m2)

    boost_vector = ROOT.Math.Boost(parent.BoostToCM())
    daughter1 = ROOT.Math.PxPyPzMVector(boost_vector(daughter1))
    daughter2 = ROOT.Math.PxPyPzMVector(boost_vector(daughter2))

    return daughter1, daughter2

if __name__ == '__main__':
    nParticles = 1000
    mass_pi_ch = 0.13957
    mass_k_zero = 0.497611
    mass_d_zero = 1.86484
    mass_b_zero = 5.27958

    random.seed(42)
    hInvMass = ROOT.TH1F("hInvMass", "Invariant Mass", 300, 0, 6)

    fileout = ROOT.TFile("tracks.root", "recreate")
    tree = ROOT.TTree('tree', 'Tree with tracks')
    tracks_vec = ROOT.std.vector(ROOT.Math.PxPyPzMVector)()
    branch = tree.Branch('tracks', tracks_vec)

    for i in tqdm.tqdm(range(nParticles)):
        tracks = []
        tracks += generate_and_decay(mass_k_zero, mass_pi_ch, mass_pi_ch)
        tracks += generate_and_decay(mass_d_zero, mass_pi_ch, mass_pi_ch)
        tracks += generate_and_decay(mass_b_zero, mass_pi_ch, mass_pi_ch)
        assert len(tracks) == 6

        for itr1 in range(len(tracks)):
            for itr2 in range(itr1 + 1, len(tracks)):
                hInvMass.Fill((tracks[itr1] + tracks[itr2]).M())

        tracks_vec.clear()
        for tr in tracks:
            tracks_vec.push_back(tr)
        tree.Fill()

    fileout.cd()
    tree.Write()
    fileout.Close()

 
    fitFunc = ROOT.TF1("fitFunc",
                       "[0]/(sqrt(2*TMath::Pi())*[2])*TMath::Exp(-0.5*((x-[1])/[2])^2)"   # K0
                       " + [3]/(sqrt(2*TMath::Pi())*[5])*TMath::Exp(-0.5*((x-[4])/[5])^2)" # D0
                       " + [6]/(sqrt(2*TMath::Pi())*[8])*TMath::Exp(-0.5*((x-[7])/[8])^2)" # B0
                       " + [9] + [10]*x + [11]*x^2 + [12]*x^3 + [13]*x^4 + [14]*x^5",     # Polynomial
                       0, 6)

    params = array('d', [
        200., 0.4976, 0.02,   # K0
        500., 1.865, 0.03,    # D0
        300., 5.2796, 0.05,   # B0
        0, 0, 0, 0, 0, 0      # Polynomial
    ])
    fitFunc.SetParameters(params)

    fitFunc.SetParLimits(1, 0.48, 0.52)   # K0 mean
    fitFunc.SetParLimits(2, 0., 0.1)      # K0 sigma
    fitFunc.SetParLimits(4, 1.80, 1.90)   # D0 mean
    fitFunc.SetParLimits(5, 0., 0.2)      # D0 sigma
    fitFunc.SetParLimits(7, 5.20, 5.35)   # B0 mean
    fitFunc.SetParLimits(8, 0., 0.5)      # B0 sigma

    hInvMass.Fit(fitFunc)

    bin_width = hInvMass.GetBinWidth(1)
    print(f'Number of signal events K0 = {fitFunc.GetParameter(0)/bin_width:.1f} ± {fitFunc.GetParError(0)/bin_width:.1f}')
    print(f'Number of signal events D0 = {fitFunc.GetParameter(3)/bin_width:.1f} ± {fitFunc.GetParError(3)/bin_width:.1f}')
    print(f'Number of signal events B0 = {fitFunc.GetParameter(6)/bin_width:.1f} ± {fitFunc.GetParError(6)/bin_width:.1f}')

    canvas = ROOT.TCanvas("canvas", "Invariant Mass", 600, 600)
    hInvMass.GetXaxis().SetTitle('M(π⁺π⁻) [GeV]')
    hInvMass.GetYaxis().SetTitle('Events')
    hInvMass.Draw()
    fitFunc.Draw("same")

    canvas.SaveAs("invmass_1.pdf")
    canvas.SaveAs("invmass_1.png")
