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
  
  std::vector <double> pdp_average;
  std::vector <double> lanl_average;

    
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
	  nmr->entry.at(i)->data.at(j) *= -1.0;
	}

	canvas->cd();
	
	TGraph *pdp  = new TGraph(steps, nmr->entry.at(i)->frequency.data(), nmr->entry.at(i)->data.data());
	pdp->SetMarkerStyle(kFullDotMedium);
	pdp->SetMarkerColor(kBlue+2);
	pdp->Draw("ac");

	canvas->SaveAs(Form("pdp_histogram%d.C", i));
	canvas->Clear();
	
	delete pdp;
      }
  }

  pdp_average = nmr->ComputeSignalAverage(nmr->entry);
  std::vector <double> fit = nmr->ComputeBackgroundSignal(nmr->entry.at(0)->data,
    							  nmr->entry.at(0)->frequency,
    							  0.0, 100.0, 400.0, 500.0);
  
  // std::vector <double> subtracted;
  
  // // fit->SetMarkerColor(kRed);
  // // fit->SetMarkerStyle(kFullDotMedium);
  // int index = 0;
  // for(std::vector <double>::iterator it = fit.begin(); it != fit.end(); it++){
  //   index = std::distance(fit.begin(), it);
  //   subtracted.push_back(nmr->entry.at(0)->data.at(index) - (*it));
  // }

  TMultiGraph *mgraph = new TMultiGraph();
  TGraph *pdp_avg  = new TGraph(steps, nmr->entry.at(0)->frequency.data(), pdp_average.data());
  pdp_avg->SetMarkerStyle(kFullDotMedium);
  pdp_avg->SetMarkerColor(kBlue+2);

  // TGraph *comp  = new TGraph(steps, nmr->entry.at(0)->frequency.data(), subtracted.data());
  // comp->SetMarkerStyle(kFullDotMedium);
  // comp->SetMarkerColor(kBlue+2);

  // canvas->cd();
  // comp->Draw("ap");
  // // mgraph->Add(avg);
  // // mgraph->Add(bck);
  // // mgraph->Draw("ap");
  // canvas->SaveAs("comp.C");
  // exit(1);

 
  // average /= total_events;
  
  // if(nmr->bExternNMRAverageSet)
  //   average = nmr->kExternNMRAverage;
  
  // // Write the data to an output file
  
  // for(std::vector<double>::iterator it = samples.begin(); it != samples.end(); it++) std_dev += std::pow(*it - average, 2);

  // std_dev /= std::pow(total_events,2);
  // std_dev = std::sqrt(std_dev);
  
  // if(nmr->bExternNMRAverageSet)
  //   std_dev = nmr->kExternNMRError;
  
  // std::cout << "PDP Average: " << average << " "
  // 	    << std_dev << std::endl;
  
  // count = 0;
  // if(nmr->bWriteOutputFile){    
  //   for(std::vector<double>::iterator it = samples.begin(); it != samples.end(); it++){
  //     nmr_file << (nmr->kDataSet) + 0.001*count  << "\t"
  // 	       << *it/average << "\t"
  // 	       << (*it/average)*std::sqrt( std::pow(std_dev/average, 2)) << std::endl;
  //     count++;
  //   }
  // }

  // count = 0;
  // average = 0;
  // std_dev = 0;
  // total_events = 0;
  // samples.clear();
  
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
	
	for(unsigned int j = 0; j < (unsigned int)(lanl->entry.at(i)->data.size()); j++){
	  lanl->entry.at(i)->frequency.push_back( (central_frequency - modulation) + j*((2*modulation)/steps) );  // f_lower = f_central - 0.4 and f_upper = f_central + 0.4
	}

	TGraph *lal  = new TGraph(steps, lanl->entry.at(i)->frequency.data(), lanl->entry.at(i)->data.data());
	lal->SetMarkerStyle(kFullDotMedium);
	lal->SetMarkerColor(kRed+2);
	lal->Draw("ac");

	canvas->SaveAs(Form("lanl_histogram%d.C", i));
	canvas->Clear();
	
	delete lal;
      }
  }

  lanl_average = lanl->ComputeSignalAverage(lanl->entry);
  TGraph *lanl_avg  = new TGraph(steps, lanl->entry.at(0)->frequency.data(), lanl_average.data());
  lanl_avg->SetMarkerStyle(kFullDotMedium);
  lanl_avg->SetMarkerColor(kRed+2);
 
  canvas->cd();
  
  mgraph->Add(pdp_avg);
  mgraph->Add(lanl_avg);
  mgraph->Draw("ap");
  canvas->SaveAs("comp.C");

  
  return(0);
}
