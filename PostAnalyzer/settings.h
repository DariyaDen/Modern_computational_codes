// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
// >>>>>>>>>>>>>> General settings (directories) >>>>>>>>>>>>>>>>>>>>>>>
// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#ifndef TTBAR_SETTINGS_H
#define TTBAR_SETTINGS_H

#include <TString.h>

// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
//
// main directory (by default current directory)
TString gBaseDir  = "./"; 
//
// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
//
// these are "normal" settings 
TString gDataDir  = gBaseDir + "./ntuples-data"; // directory with data ntuples
TString gMcDir    = gBaseDir + "./ntuples-mc"; // directory with MC ntuples
TString gHistDir  = gBaseDir + "./hist"; // directory with histograms
TString gPlotsDir = gBaseDir + "./plots"; // directory with final plots
//
// For exercises, you could use existing "reference" histograms 
// (they are provided at git) to produce final plots, or even existing 
// ROOT ntuples to produce histograms, for this uncomment proper lines 
// below (and comment out corresponding lines above, of course)
//TString gDataDir  = gBaseDir + "./ntuples-data-REF";
//TString gMcDir    = gBaseDir + "./ntuples-mc-REF";
//TString gHistDir  = gBaseDir + "./hist-REF";
//
// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#endif
