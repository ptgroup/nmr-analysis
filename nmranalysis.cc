// Standard C++ libraries
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <iostream>

// Custom libraries
#include "NMRAnalysis.hh"

// Boost libraries

// Root libraries
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TApplication.h"
#include "TAxis.h"
#include "TF1.h"
#include "TH1.h"

int main(int argc, char *argv[])
{
  
  NMRAnalysis *nmr = new NMRAnalysis();

  
  nmr->GetOptions(argv);
  if(nmr->bGraphicsShow) nmr->InitGraphicsEngine(argc, argv); 

  nmr->OpenFiles();
  nmr->ReadFiles();
    
  int event_number = nmr->GetValue<int>("EventNum", 0);


  int steps = nmr->GetValue<int>("ScanSteps", 0);

  float central_frequency = nmr->GetValue<double>("RFFreq", 0);
  
  TH1F *hist = new TH1F("hist", "hist", steps, central_frequency - 0.4, central_frequency + 0.4);

  for(unsigned int i = 0; i < (unsigned int)(nmr->entry.front()->data.size()); i++)
    {
      // only looking at the first entry for the moment ===> .front() is used
      nmr->entry.front()->frequency.push_back( (central_frequency - 0.4) + i*(0.8/steps) );  // f_lower = f_central - 0.4 and f_upper = f_central + 0.4
      hist->Fill(nmr->entry.front()->frequency.back(), nmr->entry.front()->data[i]);
    }
  
  TCanvas *canvas = new TCanvas("canvas", "canvas", 5);
  canvas->cd();
  
  hist->Draw("hist");
  std::cout << "Area: " << hist->Integral() << std::endl;

  // ****************************************************************************************
  
  if(nmr->bGraphicsShow){
    nmr->RunGraphicsEngine(); 
  }

  return(0);
}
