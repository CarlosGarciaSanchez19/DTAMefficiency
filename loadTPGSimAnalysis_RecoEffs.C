#include "TROOT.h"

void loadTPGSimAnalysis_RecoEffs()
{
  gROOT->ProcessLine(".L DTNtupleBaseAnalyzer.C++");
  gROOT->ProcessLine(".L DTNtupleTPGSimAnalyzer_RecoEff.C++");
}
