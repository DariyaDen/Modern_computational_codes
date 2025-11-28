// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
// >>>>>>>>>>>>>>>>>>> ZSolutionKinRecoDilepton (updated) >>>>>>>>>>>>
// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
struct ZSolutionKinRecoDilepton
{
  ZSolutionKinRecoDilepton(): zWeight(-1.0), zBetaT(-1.0), zBetaTbar(-1.0) {;}
  TLorentzVector zT, zTbar;  // top and antitop four-momenta
  int zBTag;                  // number of b-tagged jets
  double zWeight;             // weight of this solution
  double zBetaT;              // |p|/E for top
  double zBetaTbar;           // |p|/E for antitop
};

// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
ZSolutionKinRecoDilepton* SolveKinRecoDilepton(const TLorentzVector& lm, const TLorentzVector& lp, 
  const TLorentzVector& b, const TLorentzVector& bbar, const double metX, const double metY, 
  TH1D* hInacc = NULL, int* ambiguity = NULL)
{
  // ... existing calculations for nu/nubar ...

  // store and return best solution as ZSolutionKinRecoDilepton instance
  ZSolutionKinRecoDilepton* solution = new ZSolutionKinRecoDilepton;
  solution->zT = (nuBest + lp + b);
  solution->zTbar = (nubarBest + lm + bbar);
  solution->zWeight = weightBest;

  // --- NEW: calculate beta for top and antitop ---
  solution->zBetaT    = solution->zT.P() / solution->zT.E();
  solution->zBetaTbar = solution->zTbar.P() / solution->zTbar.E();

  return solution;
}

// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
int KinRecoDilepton(const TLorentzVector& lm, const TLorentzVector& mp, const std::vector<TLorentzVector>& jets, 
  const double metX, const double metY, TLorentzVector& t, TLorentzVector& tbar, 
  TH1D* hInacc = NULL, TH1D* hAmbig = NULL, TH1D* hBetaT = NULL, TH1D* hBetaTbar = NULL)
{
  // ... existing code for finding best solution ...

  // after selecting best solution:
  if(solved)
  {
    if(hBetaT)    hBetaT->Fill(t.P() / t.E());
    if(hBetaTbar) hBetaTbar->Fill(tbar.P() / tbar.E());
  }

  // ... clean up and return ...
  return solved;
}
