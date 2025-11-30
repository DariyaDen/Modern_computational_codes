import ROOT
import math
import random
import tqdm

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

random.seed(42)
nParticles = 1000

mass_pi_ch = 0.13957
mass_k_zero = 0.497611
mass_d_zero = 1.86484
mass_b_zero = 5.27958

mass_list = []
for _ in tqdm.tqdm(range(nParticles)):
    tracks_k0 = generate_and_decay(mass_k_zero, mass_pi_ch, mass_pi_ch)
    tracks_d0 = generate_and_decay(mass_d_zero, mass_pi_ch, mass_pi_ch)
    tracks_b0 = generate_and_decay(mass_b_zero, mass_pi_ch, mass_pi_ch)
    mass_list.append((tracks_k0[0]+tracks_k0[1]).M())
    mass_list.append((tracks_d0[0]+tracks_d0[1]).M())
    mass_list.append((tracks_b0[0]+tracks_b0[1]).M())

ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.WARNING)
x = ROOT.RooRealVar("x","Invariant Mass [GeV]",0,6)

data = ROOT.RooDataSet("data","Invariant Mass Data",ROOT.RooArgSet(x))
for m in mass_list:
    x.setVal(m)
    data.add(ROOT.RooArgSet(x))

meanK0 = ROOT.RooRealVar("meanK0","K0 mean", mass_k_zero, 0.48, 0.52)
sigmaK0 = ROOT.RooRealVar("sigmaK0","K0 sigma", 0.03, 0.01, 0.1)
gaussK0 = ROOT.RooGaussian("gaussK0","K0 gaussian",x,meanK0,sigmaK0)

meanD0 = ROOT.RooRealVar("meanD0","D0 mean", mass_d_zero, 1.8, 1.88)
sigmaD0 = ROOT.RooRealVar("sigmaD0","D0 sigma",0.05,0.01,0.1)
gaussD0 = ROOT.RooGaussian("gaussD0","D0 gaussian",x,meanD0,sigmaD0)

meanB0 = ROOT.RooRealVar("meanB0","B0 mean", mass_b_zero,5.2,5.35)
sigmaB0 = ROOT.RooRealVar("sigmaB0","B0 sigma",0.05,0.01,0.1)
gaussB0 = ROOT.RooGaussian("gaussB0","B0 gaussian",x,meanB0,sigmaB0)

nK0 = ROOT.RooRealVar("nK0","N K0",1000,0,5000)
nD0 = ROOT.RooRealVar("nD0","N D0",1000,0,5000)
nB0 = ROOT.RooRealVar("nB0","N B0",1000,0,5000)

slope = ROOT.RooRealVar("slope","background slope",-1, -10, 10)
bkg = ROOT.RooExponential("bkg","background",x,slope)
nBkg = ROOT.RooRealVar("nBkg","N background",1000,0,20000)

model = ROOT.RooAddPdf("model","Extended model",
                       ROOT.RooArgList(gaussK0,gaussD0,gaussB0,bkg),
                       ROOT.RooArgList(nK0,nD0,nB0,nBkg))

fit_result = model.fitTo(data, ROOT.RooFit.Extended(), ROOT.RooFit.Save())

print("\nFit results:")
for name, var in [("K0",nK0),("D0",nD0),("B0",nB0),("Background",nBkg)]:
    print(f"{name}: N = {var.getVal():.1f} ± {var.getError():.1f}")

c = ROOT.TCanvas("c","Invariant Mass Fit",900,700)

hist_data = ROOT.TH1F("hist_data","Invariant Mass Data", 100, 0, 6)
for m in mass_list:
    hist_data.Fill(m)
hist_data.SetLineColor(ROOT.kBlue)
hist_data.SetLineWidth(2)
hist_data.Draw("HIST")

frame = x.frame(ROOT.RooFit.Title(""))
data.plotOn(frame)  # точки даних для масштабу

model.plotOn(frame, ROOT.RooFit.LineColor(ROOT.kRed), ROOT.RooFit.LineWidth(2))  # модель червоною
model.plotOn(frame, ROOT.RooFit.Components("bkg"),
             ROOT.RooFit.LineStyle(ROOT.kDashed),
             ROOT.RooFit.LineColor(ROOT.kGreen),
             ROOT.RooFit.LineWidth(2))  # фон зеленим

frame.Draw("SAME")

c.Update()
c.SaveAs("invmass++.png")
c.SaveAs("invmass++.pdf")
print("\nSaved plots: invmass++.png / invmass++.pdf")

