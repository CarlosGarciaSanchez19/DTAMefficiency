#include "TROOT.h"

void loadTPGSimAnalysis_AMEffs()
{
  gROOT->ProcessLine(".L DTNtupleBaseAnalyzer.C++");
  gROOT->ProcessLine(".L DTNtupleTPGSimAnalyzer_AMEff.C++");
}
