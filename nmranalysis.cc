// Standard C++ libraries
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <iostream>

// Custom libraries
#include "NMRAnalysis.hh"
// #include "DataFile.hh"
// #include "FontColor.hh"

// Boost libraries

// Root libraries
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TApplication.h"
#include "TAxis.h"
#include "TF1.h"
// #include "TFitResultPtr.h"
// #include "TVirtualFitter.h"
// #include "TMath.h"
// #include "Math/MinimizerOptions.h"

int main(int argc, char *argv[])
{

  NMRAnalysis *nmr = new NMRAnalysis();
  nmr->GetOptions(argv);
  if(nmr->bGraphicsShow) nmr->InitGraphicsEngine(argc, argv); 


  nmr->OpenFiles();
  nmr->ReadFiles();

  if(nmr->bGraphicsShow){
    nmr->RunGraphicsEngine(); 
  }

  return(0);
}
