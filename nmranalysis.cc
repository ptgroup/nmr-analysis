// Standard C++ libraries
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <iostream>

// Custom libraries
#include "NMRAnalysis.hh"

// Root libraries
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TApplication.h"
#include "TAxis.h"
#include "TF1.h"
#include "TH1.h"
#include "TSpectrum.h"
#include "TStyle.h"

int main(int argc, char *argv[])
{
  
  NMRAnalysis *nmr = new NMRAnalysis();

  
  nmr->GetOptions(argv);
  if(nmr->bGraphicsShow) nmr->InitGraphicsEngine(argc, argv); 

  std::vector <boost::filesystem::path> data_files = nmr->GetFileList("data/nmr_data/data/", ".csv");
  std::vector <boost::filesystem::path> settings_files = nmr->GetFileList("data/nmr_data/settings/", ".csv");

  nmr->OpenDataFile(data_files.at(0).string().c_str());
  nmr->OpenSettingsFile(settings_files.at(0).string().c_str());
  nmr->ReadConfigurationMap("config_map.cfg");
  nmr->ReadNMRFiles();
    
  int event_number;
  int steps;
  
  double  modulation;
  float central_frequency;

  TCanvas *canvas = new TCanvas("canvas", "canvas", 5);
  canvas->Divide(1,2);
  canvas->cd(1);

  gStyle->SetOptStat(0);
  
  std::fstream out_file;
  out_file.open(Form("output/%s", (nmr->fOutputFile.c_str())), std::fstream::in | std::fstream::app | std::fstream::out);
  if( !(out_file.good()) ){
    std::cerr << __FUNCTION__ << " >> Error opening output file: " << std::endl;
    exit(1);
  }

  for(unsigned int i = 0; i < (unsigned int)(nmr->entry.size()); i++)
    {
      event_number = nmr->GetValue<int>("EventNum", i);
      steps = nmr->GetValue<int>("ScanSteps", 0);
      modulation = (1e-3)*(nmr->GetValue<double>("RFMod", 0));          // The units from the configuration file are in MHz
      central_frequency = nmr->GetValue<double>("RFFreq", 0);
      
      TH1F *hist = new TH1F("hist", "hist", steps, central_frequency - modulation, central_frequency + modulation);
      
      for(unsigned int j = 0; j < (unsigned int)(nmr->entry.at(i)->data.size()); j++){
	nmr->entry.at(i)->frequency.push_back( (central_frequency - modulation) + j*((2*modulation)/steps) );  // f_lower = f_central - 0.4 and f_upper = f_central + 0.4
	hist->Fill(nmr->entry.at(i)->frequency.back(), (nmr->kScaleFactor)*(nmr->entry.at(i)->data[j]));
      }
      hist->Smooth(1);
      hist->SetLineWidth(2);
      hist->SetTitle(Form("Average Proton NMR Signal: Event# %d", event_number));
      hist->Draw("hist");

      TH1F *background = (TH1F *)hist->ShowBackground(200, "same");
      TH1F *difference = (TH1F *)hist->Clone();
      difference->SetName("diff");
      difference->Add(background, -1.0);
  
      canvas->cd(2);
      difference->SetLineColor(kBlue+2);
      difference->SetFillColor(kBlue-8);
      difference->Draw("hist");
      
      std::cout << ">>>>> Area: " << difference->Integral() << std::endl;
      out_file << event_number << "," << difference->Integral() << std::endl;
      canvas->SaveAs(Form("output/nmr_hist_%i.C", event_number));
      
      delete hist;
      delete background;
      delete difference;
    } 

  if(nmr->bGraphicsShow){
    nmr->RunGraphicsEngine(); 
  }

  return(0);
}
