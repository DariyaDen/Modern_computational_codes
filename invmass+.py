import ROOT
import math
import random
import tqdm
from array import array

# -------------------------------
# Функції генерації і розпаду
# -------------------------------
def generate_and_decay(parent_mass, daughter_mass_1, daughter_mass_2):
    momentum = random.expovariate(1.0)
    theta = random.uniform(0, math.pi)
    phi = random.uniform(0, 2 * math.pi)

    px = momentum * math.sin(theta) * math.cos(phi)
    py = momentum * math.sin(theta) * math.sin(phi)
    pz = momentum * math.cos(theta)

    parent = ROOT.Math.PxPyPzMVector(px, py, pz, parent_mass)
    daughter1, daughter2 = two_body_decay(parent, daughter_mass_1, daughter_mass_2)

    daughter1 = smear_momentum(daughter1)
    daughter2 = smear_momentum(daughter2)
    return [daughter1, daughter2]

def smear_momentum(p):
    fluct = random.normalvariate(1., 0.05)
    return ROOT.Math.PxPyPzMVector(p.X()*fluct, p.Y()*fluct, p.Z()*fluct, p.M())

def two_body_decay(parent, m1, m2):
    m0 = parent.M()
    E1 = (m0**2 + m1**2 - m2**2)/(2*m0)
    p = math.sqrt(max(0.0, E1**2 - m1**2))

    theta = random.uniform(0, math.pi)
    phi = random.uniform(0, 2*math.pi)

    p1x = p*math.sin(theta)*math.cos(phi)
    p1y = p*math.sin(theta)*math.sin(phi)
    p1z = p*math.cos(theta)

    d1 = ROOT.Math.PxPyPzMVector(p1x,p1y,p1z,m1)
    d2 = ROOT.Math.PxPyPzMVector(-p1x,-p1y,-p1z,m2)

    boost_vector = ROOT.Math.Boost(parent.BoostToCM())
    d1 = ROOT.Math.PxPyPzMVector(boost_vector(d1))
    d2 = ROOT.Math.PxPyPzMVector(boost_vector(d2))
    return d1, d2

# -------------------------------
# Основна програма
# -------------------------------
if __name__ == '__main__':
    random.seed(42)
    nParticles = 1000

    mass_pi_ch = 0.13957
    mass_k_zero = 0.497611
    mass_d_zero = 1.86484
    mass_b_zero = 5.27958

    hInvMass = ROOT.TH1F("hInvMass","Invariant Mass",300,0,6)

    # Генерація частинок
    for i in tqdm.tqdm(range(nParticles)):
        tracks_k0 = generate_and_decay(mass_k_zero, mass_pi_ch, mass_pi_ch)
        tracks_d0 = generate_and_decay(mass_d_zero, mass_pi_ch, mass_pi_ch)
        tracks_b0 = generate_and_decay(mass_b_zero, mass_pi_ch, mass_pi_ch)
        all_tracks = tracks_k0 + tracks_d0 + tracks_b0

        for itr1 in range(len(all_tracks)):
            for itr2 in range(itr1+1, len(all_tracks)):
                hInvMass.Fill((all_tracks[itr1]+all_tracks[itr2]).M())

    # -------------------------------
    # Повний фіт (K0 + D0 + B0 + поліном)
    # -------------------------------
    fitFull = ROOT.TF1("fitFull",
                       "[0]/(sqrt(2*TMath::Pi())*[2])*TMath::Exp(-0.5*((x-[1])/[2])^2)"   # K0
                       " + [3]/(sqrt(2*TMath::Pi())*[5])*TMath::Exp(-0.5*((x-[4])/[5])^2)" # D0
                       " + [6]/(sqrt(2*TMath::Pi())*[8])*TMath::Exp(-0.5*((x-[7])/[8])^2)" # B0
                       " + [9] + [10]*x + [11]*x^2 + [12]*x^3 + [13]*x^4 + [14]*x^5",
                       0,6)

    params_full = array('d', [
        200., 0.4976, 0.02,  # K0
        500., 1.865, 0.03,   # D0
        300., 5.2796, 0.05,  # B0
        0.,0.,0.,0.,0.,0.    # Поліном
    ])
    fitFull.SetParameters(params_full)
    fitFull.SetParLimits(1,0.48,0.52)
    fitFull.SetParLimits(2,0.,0.1)
    fitFull.SetParLimits(4,1.80,1.90)
    fitFull.SetParLimits(5,0.,0.2)
    fitFull.SetParLimits(7,5.20,5.35)
    fitFull.SetParLimits(8,0.,0.5)
    hInvMass.Fit(fitFull,"R")

    # -------------------------------
    # Зелений фіт (K0 + D0 + поліном, без B0)
    # -------------------------------
    fitKD = ROOT.TF1("fitKD",
                     "[0]/(sqrt(2*TMath::Pi())*[2])*TMath::Exp(-0.5*((x-[1])/[2])^2)"  # K0
                     " + [3]/(sqrt(2*TMath::Pi())*[5])*TMath::Exp(-0.5*((x-[4])/[5])^2)" # D0
                     " + [6] + [7]*x + [8]*x^2 + [9]*x^3 + [10]*x^4 + [11]*x^5",
                     0,6)

    params_kd = array('d', [
        fitFull.GetParameter(0),  # K0 amp
        fitFull.GetParameter(1),  # K0 mean
        fitFull.GetParameter(2),  # K0 sigma
        fitFull.GetParameter(3),  # D0 amp
        fitFull.GetParameter(4),  # D0 mean
        fitFull.GetParameter(5),  # D0 sigma
        fitFull.GetParameter(9),  # Poly0
        fitFull.GetParameter(10), # Poly1
        fitFull.GetParameter(11), # Poly2
        fitFull.GetParameter(12), # Poly3
        fitFull.GetParameter(13), # Poly4
        fitFull.GetParameter(14), # Poly5
    ])
    fitKD.SetParameters(params_kd)
    fitKD.SetLineColor(ROOT.kGreen+2)
    fitKD.SetLineWidth(3)

    # -------------------------------
    # Малюємо все на одному канвасі
    # -------------------------------
    canvas = ROOT.TCanvas("canvas","Invariant Mass",600,600)
    hInvMass.GetXaxis().SetTitle("M(π⁺π⁻) [GeV]")
    hInvMass.GetYaxis().SetTitle("Events")
    hInvMass.Draw()
    fitFull.SetLineColor(ROOT.kRed)
    fitFull.SetLineWidth(3)
    fitFull.Draw("same")
    fitKD.Draw("same")

    canvas.SaveAs("invmass_full_vs_KD.pdf")
    canvas.SaveAs("invmass_full_vs_KD.png")

