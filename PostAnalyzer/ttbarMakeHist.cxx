// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
// This code processes ROOT ntuples for ttbar analysis (see 
// Analyzer/src/Analyzer.cc) and produces histograms, which are 
// further used to make final plots (see ttbarMakePlots.cxx).
// Run: ./ttbarMakeHist
// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
//
// additional files from this analysis (look there for description) 
#include "eventReco.h"
#include "settings.h"
//
// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
// >>>>>>>>>>>>>>>>>>>>> Main function >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
int main(int argc, char** argv)
{
  //
  // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  // >>>>>>>>>>>>>>>>>>>>> Settings >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  //
  // set directories to data and MC ntuples
  TString dataDir = gDataDir;
  TString mcDir = gMcDir;
  //
  // flags what to run
  bool flagData    = 1; // if 1, data will be processed
  bool flagMCsig   = 1; // if 1, signal MC (dileptonic decay channel) will be processed
  bool flagMCother = 1; // if 1, signal MC 'other' decay channels will be processed to form background MC histograms
  bool flagMCstop  = 1; // if 1, MC single top (background) will be processed
  bool flagMCwjets = 1; // if 1, MC W+jets (background) will be processed
  bool flagMCdy    = 1; // if 1, MC Drell-Yan (background) will be processed
  //
  // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  //
  // common purpose variables
  std::vector<TString> nameInFile; // container to store input file names
  
  // histograms
  TH1::SetDefaultSumw2(); // keep histogram weights by default
  // ZVarHisto is a simple class which incorporates a histogram and a variable name. 
  // This class is used to store needed input settings (variable names, binning) 
  // for control plots and cross sections (as in TOP-11-013)
  std::vector<ZVarHisto> vecVH, vecVHGen; // vecVH for reconstruction level, vecVHGen for generator level
  
  // histograms and variables for control plots
  vecVHGen.push_back(ZVarHisto("ptt", new TH1D("h_ptt", "pT top", 21, 0.0, 420.0)));
  vecVHGen.push_back(ZVarHisto("ptat", new TH1D("h_ptat", "pT atop", 21, 0.0, 420.0)));
  vecVHGen.push_back(ZVarHisto("pttat", new TH1D("h_pttat", "pT tatop", 21, 0.0, 420.0)));
  vecVHGen.push_back(ZVarHisto("pttt", new TH1D("h_pttt", "pT ttbar", 30, 0.0, 300.0)));
  vecVHGen.push_back(ZVarHisto("yt", new TH1D("h_yt", "y top", 26, -2.6, 2.6)));
  vecVHGen.push_back(ZVarHisto("yat", new TH1D("h_yat", "y atop", 26, -2.6, 2.6)));
  vecVHGen.push_back(ZVarHisto("ytat", new TH1D("h_ytat", "y tatop", 26, -2.6, 2.6)));
  vecVHGen.push_back(ZVarHisto("ytt", new TH1D("h_ytt", "y ttbar", 26, -2.6, 2.6)));
  // <<< NEW >>> phi(ttbar)
  vecVHGen.push_back(ZVarHisto("phi_tt", new TH1D("h_phi_tt", "phi ttbar", 32, -3.2, 3.2)));
  
  // histograms and variables for cross sections
  {
    double bins[] = {0.,80.,130.,200.,300.,400.};
    vecVHGen.push_back(ZVarHisto("ptt", new TH1D("h_ptt_cs", "pT top", 5, bins)));
    vecVHGen.push_back(ZVarHisto("ptat", new TH1D("h_ptat_cs", "pT atop", 5, bins)));
    vecVHGen.push_back(ZVarHisto("pttat", new TH1D("h_pttat_cs", "pT tatop", 5, bins)));
  }
  {
    double bins[] = {-2.5,-1.3,-0.8,-0.4,0.0,0.4,0.8,1.3,2.5};
    vecVHGen.push_back(ZVarHisto("yt", new TH1D("h_yt_cs", "y top", 8, bins)));
    vecVHGen.push_back(ZVarHisto("yat", new TH1D("h_yat_cs", "y atop", 8, bins)));
    vecVHGen.push_back(ZVarHisto("ytat", new TH1D("h_ytat_cs", "y tatop", 8, bins)));
  }
  {
    double bins[] = {0.,20.,60.,120.,300.};
    vecVHGen.push_back(ZVarHisto("pttt", new TH1D("h_pttt_cs", "pT ttbar", 4, bins)));
  }
  {
    double bins[] = {-2.5,-1.5,-0.7,0.0,0.7,1.5,2.5};
    vecVHGen.push_back(ZVarHisto("ytt", new TH1D("h_ytt_cs", "y ttbar", 6, bins)));
  }
  {
    double bins[] = {0.,345.,400.,470.,550.,650.,800.,1100.,1600.};
    vecVHGen.push_back(ZVarHisto("mtt", new TH1D("h_mtt_cs", "M ttbar", 8, bins)));
  }
  // <<< NEW >>> phi(ttbar) for cross sections
  vecVHGen.push_back(ZVarHisto("phi_tt", new TH1D("h_phi_tt_cs", "phi ttbar", 32, -3.2, 3.2)));
  
  // for reconstruction level the same binning is needed
  vecVH = vecVHGen;
  // add lepton pT histogram at reconstruction level
  vecVH.push_back(ZVarHisto("ptl", new TH1D("h_ptl", "pT leptons", 23, 30.0, 260.0)));
  
  // loop over decay channels (ch = 1 ee, ch = 2 mumu, ch = 3 emu)
  for(int ch = 1; ch <= 3; ch++)
  {
    // *****************************************
    // **************** DATA *******************
    // *****************************************
    if(flagData)
    {
      ZEventRecoInput in;
      in.Name = "data";
      in.Type = 1;
      in.Channel = ch;
      in.VecVarHisto = vecVH;
      if(ch == 1) in.AddToChain(dataDir + "/DoubleElectron/*.root");
      else if(ch == 2) in.AddToChain(dataDir + "/DoubleMu/*.root");
      else if(ch == 3) in.AddToChain(dataDir + "/MuEG/*.root");
      eventreco(in);
    }
    
    // *****************************************
    // ************** MC signal ****************
    // *****************************************
    if(flagMCsig)
    {
      ZEventRecoInput in;
      in.Weight = 0.007529;
      in.Name = "mcSigReco";
      in.Type = 2;
      in.Channel = ch;
      in.VecVarHisto = vecVH;
      in.AddToChain(mcDir + "/TTJets_TuneZ2_7TeV-madgraph-tauola/00001/*.root");
      in.AddToChain(mcDir + "/TTJets_TuneZ2_7TeV-madgraph-tauola/010000/*.root");
      in.AddToChain(mcDir + "/TTJets_TuneZ2_7TeV-madgraph-tauola/010003/*.root");
      in.AddToChain(mcDir + "/TTJets_TuneZ2_7TeV-madgraph-tauola/010002/*.root");
      in.AddToChain(mcDir + "/TTJets_TuneZ2_7TeV-madgraph-tauola/010001/*.root");
      in.AddToChain(mcDir + "/TTJets_TuneZ2_7TeV-madgraph-tauola/00000/*.root");
      eventreco(in);
      
      in.Name = "mcSigOtherReco";
      in.Type = 3;
      eventreco(in);
      
      in.Name = "mcSigGen";
      in.Type = 2;
      in.VecVarHisto = vecVHGen;
      in.Gen = true;
      eventreco(in);
    }
    
    // *****************************************
    // ************ MC single top **************
    // *****************************************
    if(flagMCstop)
    {
      ZEventRecoInput in;
      in.Name = "mcSingleTopReco";
      in.Type = 4;
      in.Weight = 0.02544;
      in.Channel = ch;
      in.VecVarHisto = vecVH;
      in.AddToChain(mcDir + "/Tbar_TuneZ2_tW-channel-DR_7TeV-powheg-tauola/*.root");
      in.AddToChain(mcDir + "/T_TuneZ2_tW-channel-DR_7TeV-powheg-tauola/*.root");
      eventreco(in);
    }
    
    // *****************************************
    // ************** MC W+jets ****************
    // *****************************************
    if(flagMCwjets)
    {
      ZEventRecoInput in;
      in.Name = "mcWjetsReco";
      in.Weight = 0.3197;
      in.Type = 4;
      in.Channel = ch;
      in.VecVarHisto = vecVH;
      in.AddToChain(mcDir + "/WJetsToLNu_TuneZ2_7TeV-madgraph-tauola/*.root");
      eventreco(in);
    }
    
    // *****************************************
    // **************** MC DY ******************
    // *****************************************
    if(flagMCdy)
    {
      ZEventRecoInput in;
      in.Type = 4;
      in.Channel = ch;
      in.VecVarHisto = vecVH;
      in.Name = "mcDYlmReco";
      in.Weight = 0.07459;
      in.AddToChain(mcDir + "/DYJetsToLL_M-10To50_TuneZ2_7TeV-pythia6/*.root");
      eventreco(in);
      
      in.Name = "mcDYhmReco";
      in.Weight = 0.2093;
      in.ClearChain();
      in.AddToChain(mcDir + "/DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola/*.root");
      eventreco(in);
    }
  }

  return 0;
}
// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
