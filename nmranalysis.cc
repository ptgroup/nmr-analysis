// Standard C++ libraries
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <iostream>

// Custom libraries
#include "NMRAnalysis.hh"

// Root libraries
#include "TGraph.h"
#include "TMultiGraph.h"
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
  
  std::vector <boost::filesystem::path> data_files = nmr->GetFileList("data/nmr_data/data/", ".csv");
  std::vector <boost::filesystem::path> settings_files = nmr->GetFileList("data/nmr_data/settings/", ".csv");

  // ***********************************************************
  // Sort is needed due to the unfortunate fact that the order
  // of the file lists is not the same. It might be better to
  // have a custom file list configuration file or just build
  // the data file extensions from the settings files list.
  // The former is probably best but for now this sorts the
  // data file list so that it matches the setting file list
  // order. Of course the files my be coresponding.
  // ***********************************************************

  nmr->Sort(settings_files, data_files);   
  
  int event_number = 0;
  int steps = 500;
  
  float modulation;
  float central_frequency;
  
  std::vector <double> average;

    
  TCanvas *canvas = new TCanvas("canvas", "canvas", 5);

  std::fstream nmr_file;

  nmr_file.open(Form("output/nmr_%s", (nmr->fOutputFile.c_str())), std::fstream::in | std::fstream::app | std::fstream::out);
  if( !(nmr_file.good()) ){
    std::cerr << __FUNCTION__ << " >> Error opening output file: " << nmr->fOutputFile << std::endl;
    exit(1);
  }
  
  for(unsigned int file = 0; file < (unsigned int)data_files.size(); file++){

    nmr->OpenDataFile((std::string("data/nmr_data/data/") + data_files.at(file).string()).c_str());
    nmr->OpenSettingsFile((std::string("data/nmr_data/settings/") + settings_files.at(file).string()).c_str());
    nmr->ReadConfigurationMap("config_map.cfg");
    nmr->ReadNMRFiles();
    
    // canvas->Divide(1,2);
    // canvas->cd(1);
 
    
    gStyle->SetOptStat(0);

    // for(unsigned int i = 0; i < (unsigned int)(nmr->entry.size()); i++)
      for(unsigned int i = 0; i < 1; i++)
      {
	event_number = nmr->GetValue<int>("EventNum", 0);
	steps = nmr->GetValue<int>("ScanSteps", 0);
	modulation = (1e-3)*(nmr->GetValue<double>("RFMod", 0));          // The units from the configuration file are in MHz
	central_frequency = nmr->GetValue<double>("RFFreq", 0);

	for(unsigned int j = 0; j < (unsigned int)(nmr->entry.at(i)->data.size()); j++){
	  nmr->entry.at(i)->frequency.push_back( (central_frequency - modulation) + (j+1)*((2*modulation)/steps) );  // f_lower = f_central - 0.4 and f_upper = f_central + 0.4
	}

	canvas->cd();	
	TGraph *pdp  = new TGraph(steps, nmr->entry.at(i)->frequency.data(), nmr->entry.at(i)->data.data());
	pdp->SetMarkerStyle(kFullDotMedium);
	pdp->SetMarkerColor(kBlue+2);
	pdp->Draw("ac");

	canvas->SaveAs(Form("histogram%d.C", i));
	canvas->Clear();
	
	delete pdp;
      }
  }
  average = nmr->ComputeSignalAverage(nmr->entry);
  TGraph *bck = nmr->ComputeBackgroundSignal(nmr->entry.at(0)->data,
					     nmr->entry.at(0)->frequency,
					     0.0, 100.0, 400.0, 500.0);
  bck->SetMarkerColor(kRed);
  bck->SetMarkerStyle(kFullDotMedium);

  // ************************************************

  std::vector <double> x;
  std::vector <double> y;
  
  for(int i = -10; i < 11; i++){
    x.push_back(i);
    y.push_back(200.0 + 2.0*x.back() + 2.0*std::pow(x.back(), 2) );

    std::cout << ">>>> " << y.back() << " " << x.back() << std::endl;
  }

  nmr->PolynomialRegression(x, y);

  

  // ************************************************
  
  TMultiGraph *mgraph = new TMultiGraph();
  TGraph *avg  = new TGraph(steps, nmr->entry.at(0)->frequency.data(), nmr->entry.at(0)->data.data());
  avg->SetMarkerStyle(kFullDotMedium);
  avg->SetMarkerColor(kBlue+2);

  canvas->cd();
  mgraph->Add(avg);
  mgraph->Add(bck);
  mgraph->Draw("ap");
  canvas->SaveAs("comp.C");
  
  exit(1);
  /*  
  average /= total_events;
  
  if(nmr->bExternNMRAverageSet)
    average = nmr->kExternNMRAverage;
  
  // Write the data to an output file
  
  for(std::vector<double>::iterator it = samples.begin(); it != samples.end(); it++) std_dev += std::pow(*it - average, 2);

  std_dev /= std::pow(total_events,2);
  std_dev = std::sqrt(std_dev);
  
  if(nmr->bExternNMRAverageSet)
    std_dev = nmr->kExternNMRError;
  
  std::cout << "PDP Average: " << average << " "
	    << std_dev << std::endl;
  
  count = 0;
  if(nmr->bWriteOutputFile){    
    for(std::vector<double>::iterator it = samples.begin(); it != samples.end(); it++){
      nmr_file << (nmr->kDataSet) + 0.001*count  << "\t"
	       << *it/average << "\t"
	       << (*it/average)*std::sqrt( std::pow(std_dev/average, 2)) << std::endl;
      count++;
    }
  }

  count = 0;
  average = 0;
  std_dev = 0;
  total_events = 0;
  samples.clear();
  
  // Begin analysis of lanl data
  
  event_number = 0;
  steps = 0;
  modulation = 0;
  central_frequency = 0;
  
  NMRAnalysis *lanl = new NMRAnalysis();
  
  lanl->GetOptions(argv);

  std::vector <boost::filesystem::path> lanl_data_files = lanl->GetFileList("data/lanl_data/data/", ".csv");
  std::vector <boost::filesystem::path> lanl_settings_files = lanl->GetFileList("data/lanl_data/settings/", ".csv");

  lanl->Sort(lanl_settings_files, lanl_data_files);   


  std::fstream lanl_file;
  lanl_file.open(Form("output/lanl_%s", (lanl->fOutputFile.c_str())), std::fstream::in | std::fstream::app | std::fstream::out);
  if( !(lanl_file.good()) ){
    std::cerr << __FUNCTION__ << " >> Error opening output file: " << std::endl;
    exit(1);
  }

  for(unsigned int file = 0; file < (unsigned int)lanl_data_files.size(); file++){
  
    lanl->OpenDataFile((std::string("data/lanl_data/data/")+lanl_data_files.at(file).string()).c_str());
    lanl->OpenSettingsFile((std::string("data/lanl_data/settings/") + lanl_settings_files.at(file).string()).c_str());
    lanl->ReadConfigurationMap("lanl_map.cfg");
    lanl->ReadNMRFiles();
    
    // for(unsigned int i = 0; i < (unsigned int)(lanl->entry.size()); i++)
    for(unsigned int i = 0; i < 1; i++)
      {
	event_number = lanl->GetValue<int>("sweep_number", i);
	steps = lanl->GetValue<int>("sweeps", 0);
	central_frequency = lanl->GetValue<double>("central_frequency", 0);
	modulation = 0.5*steps*(lanl->GetValue<double>("step_size", 0));          
	
	TH1D *hist = new TH1D("hist", "hist", steps, central_frequency - modulation, central_frequency + modulation);
	
	for(unsigned int j = 0; j < (unsigned int)(lanl->entry.at(i)->data.size()); j++){
	  lanl->entry.at(i)->frequency.push_back( (central_frequency - modulation) + j*((2*modulation)/steps) );  // f_lower = f_central - 0.4 and f_upper = f_central + 0.4
	  hist->Fill(lanl->entry.at(i)->frequency.back(), (lanl->kScaleFactor)*(lanl->entry.at(i)->data[j])/0.0579478);
	}
	hist->Smooth(1);
	hist->SetLineWidth(2);
	hist->SetTitle(Form("Average Proton NMR Signal: Event# %d", event_number));
	hist->Draw("hist");

	hist->GetXaxis()->SetRange(20,500);
	
	TH1D *background = (TH1D *)hist->ShowBackground(200, "same");
	TH1D *difference = (TH1D *)hist->Clone();
	difference->SetName("diff");
	difference->Add(background, -1.0);
	
	canvas->cd(2);
	difference->SetLineColor(kBlue+2);
	difference->SetFillColor(kBlue-8);
	difference->Draw("hist");

	lal = new TH1D("lanl", "lanl", steps, central_frequency - modulation, central_frequency + modulation);
	lal = (TH1D *)difference->Clone();
	
	spectrum->Search(difference, 2, "goff", 0.01);
	
	std::cout << "Peaks found LANL: " << spectrum->GetNPeaks() << std::endl;
	for(int j = 0; j < 4; j++){
	  std::cout <<  spectrum->GetPositionX()[j] << " " << spectrum->GetPositionY()[j] << std::endl;
	}

	canvas->SaveAs(Form("histogram_lanl%d.C", i));
	
	average += difference->Integral();
	
	samples.push_back(difference->Integral());
	
	delete hist;
	delete background;
	delete difference;
      }
    total_events += lanl->entry.size();    
    lanl->Clear();
  }

  average /= total_events;
    
  if(lanl->bExternLANLAverageSet)
    average = lanl->kExternLANLAverage;
    
  for(std::vector<double>::iterator it = samples.begin(); it != samples.end(); it++) std_dev += std::pow(*it - average, 2);
  std_dev /= std::pow(total_events,2);
  std_dev = std::sqrt(std_dev);
  
  if(lanl->bExternLANLAverageSet)
    std_dev = lanl->kExternLANLError;
  
  std::cout << "LANL Average: " << average << " "
	    << std_dev << std::endl;
  
  count = 0;
  
  if(lanl->bWriteOutputFile){      
    for(std::vector<double>::iterator it = samples.begin(); it != samples.end(); it++){
      lanl_file << (nmr->kDataSet) + 0.001*count  << "\t"
		<< *it/average << "\t"
		<< (*it/average)*std::sqrt( std::pow(std_dev/average, 2)) << std::endl;
      count++;
    }
  }

  canvas->Divide(1,1);
  canvas->cd();

  pdp->Draw("hist");
  pdp->SetLineColor(kBlue+2);
  pdp->SetFillColor(kBlue-8);

  lal->Draw("hist same");
  pdp->SetLineColor(kRed+2);
  pdp->SetFillColor(kRed-8);
  
  canvas->SaveAs("comp.C");
      
  average = 0;
  std_dev = 0;
  samples.clear();
  */  
  return(0);
}
