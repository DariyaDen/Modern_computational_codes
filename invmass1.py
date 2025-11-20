import ROOT
import math
import random
import tqdm
import sys
from array import array

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

def smear_momentum(particle):
    fluct = random.normalvariate(1., 0.05)
    return ROOT.Math.PxPyPzMVector(particle.X() * fluct,
                                   particle.Y() * fluct,
                                   particle.Z() * fluct,
                                   particle.M())

def two_body_decay(parent, m1, m2):
    m0 = parent.M()
    # protect kinematics
    if m0 < (m1 + m2):
        return (ROOT.Math.PxPyPzMVector(0,0,0,m1), ROOT.Math.PxPyPzMVector(0,0,0,m2))

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
        b0_daughters = generate_and_decay(mass_b_zero, mass_pi_ch, mass_pi_ch)
        tracks += b0_daughters

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


    low_b0 = 4.7
    high_b0 = 5.35
    fitB0 = ROOT.TF1("fitB0", "gaus", low_b0, high_b0)

    
    ibin_low = hInvMass.FindBin(low_b0)
    ibin_high = hInvMass.FindBin(high_b0)
    max_height = 0.0
    max_binx = (low_b0 + high_b0) / 2.0
    for b in range(ibin_low, ibin_high+1):
        val = hInvMass.GetBinContent(b)
        if val > max_height:
            max_height = val
            max_binx = hInvMass.GetBinCenter(b)

    fitB0.SetParameters(max_height if max_height>0 else 100., 5.28, 0.1)

    hInvMass.Fit(fitB0, "RQ")

    b0_height = fitB0.GetParameter(0)
    b0_mean = fitB0.GetParameter(1)
    b0_sigma = abs(fitB0.GetParameter(2))

    # конвертація висоти->площі для parametrization area/(sqrt(2π)*σ)*exp(...)
    b0_area = b0_height * b0_sigma * math.sqrt(2.0 * math.pi)


    # Формула: area/(sqrt(2π)*σ) * exp(-0.5*((x-mean)/σ)^2) для кожного піку
    fit_expr = (
        "[0]/(sqrt(2*TMath::Pi())*[2])*TMath::Exp(-0.5*((x-[1])/[2])^2)"   # K0 (area,[1]=mean,[2]=sigma)
        " + [3]/(sqrt(2*TMath::Pi())*[5])*TMath::Exp(-0.5*((x-[4])/[5])^2)" # D0
        " + [6]/(sqrt(2*TMath::Pi())*[8])*TMath::Exp(-0.5*((x-[7])/[8])^2)" # B0
        " + [9] + [10]*x + [11]*x*x + [12]*x*x*x"                           # poly3
    )
    fitFunc = ROOT.TF1("fitFunc", fit_expr, 0., 6.)

    params = array('d', [
        200., 0.4976, 0.02,    # [0],[1],[2]  K0: area, mean, sigma
        500., 1.865, 0.03,     # [3],[4],[5]  D0
        max(2000., b0_area), b0_mean, max(0.02, b0_sigma),  # [6],[7],[8] B0: area(from prev fit), mean, sigma
        0., 0., 0., 0.         # [9..12] poly3 coeffs
    ])
    fitFunc.SetParameters(params)

    fitFunc.SetParLimits(1, 0.485, 0.51)   # K0 mean
    fitFunc.SetParLimits(2, 0.001, 0.1)    # K0 sigma
    fitFunc.SetParLimits(4, 1.80, 1.90)    # D0 mean
    fitFunc.SetParLimits(5, 0.001, 0.3)    # D0 sigma
    fitFunc.SetParLimits(7, 5.20, 5.35)    # B0 mean
    fitFunc.SetParLimits(8, 0.02, 0.8)     # B0 sigma (тепер до 0.8)

    # Виконуємо основний фіт кілька разів (for стабільності)
    hInvMass.Fit(fitFunc, "R")
    hInvMass.Fit(fitFunc, "R")
    hInvMass.Fit(fitFunc, "R")


    bin_width = hInvMass.GetBinWidth(1)
    area_K0 = fitFunc.GetParameter(0)
    area_D0 = fitFunc.GetParameter(3)
    area_B0 = fitFunc.GetParameter(6)

    err_area_K0 = fitFunc.GetParError(0)
    err_area_D0 = fitFunc.GetParError(3)
    err_area_B0 = fitFunc.GetParError(6)

    print(f'K0: area = {area_K0:.1f} ± {err_area_K0:.1f}  (events)')
    print(f'D0: area = {area_D0:.1f} ± {err_area_D0:.1f}  (events)')
    print(f'B0: area = {area_B0:.1f} ± {err_area_B0:.1f}  (events)')

    canvas = ROOT.TCanvas("canvas", "Invariant Mass", 700, 600)
    hInvMass.GetXaxis().SetTitle('M(π⁺π⁻) [GeV]')
    hInvMass.GetYaxis().SetTitle('Events')
    hInvMass.Draw()
    fitFunc.SetLineColor(ROOT.kRed)
    fitFunc.Draw("same")

    canvas.SaveAs("invmass_1_improved.pdf")
    canvas.SaveAs("invmass_1_improved.png")

    # Додатково: виведемо параметри B0
    print(f'B0 fit (from combined): mean = {fitFunc.GetParameter(7):.4f} GeV, sigma = {fitFunc.GetParameter(8):.4f} GeV')

