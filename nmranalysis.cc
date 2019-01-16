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
#include "TLegend.h"

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

    for(unsigned int i = 0; i < (unsigned int)(nmr->entry.size()); i++)
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

	//	canvas->SaveAs(Form("output/pdp_histogram%d.C", i));
	canvas->Clear();
	
	delete pdp;
      }
  }

  pdp_average = nmr->ComputeSignalAverage(nmr->entry, "pdp");
  
  std::vector <double> fit = nmr->ComputeBackgroundSignal(pdp_average,
    							  nmr->entry.at(0)->frequency,
    							  0.0, 100.0, 400.0, 500.0, "pdp");
  std::vector <double> subtracted;
  
  int index = 0;
  double offset = 0;
  
  for(std::vector <double>::iterator it = fit.begin(); it != fit.end(); it++){
    index = std::distance(fit.begin(), it);
    subtracted.push_back(pdp_average.at(index) - (*it));
  }

  if(std::abs(subtracted.front()) > 1e-6){
    offset = subtracted.front();
    for(std::vector <double>::iterator it = fit.begin(); it != fit.end(); it++){
      index = std::distance(fit.begin(), it);
      subtracted.at(index) = (pdp_average.at(index) - (*it) - offset);
    }
  }

  pdp_average = nmr->ScaleData(-1.0, pdp_average);
  subtracted = nmr->ScaleData(-1.0, subtracted);

  // *******************************************************************************************
  //
  // Calculate the residuals from the background subtraction and use them to estimate a lower 
  // bound for the area calculation.
  //
  // *******************************************************************************************

  double residual = 0;
  double area_error = 0;
  
  // Get residuals from subtracted only in the region we used to fit the background. These values should all be zero
  // if our fit was perfect.
  
  for(int i = 0; i < 100; i++){
    residual += std::pow(subtracted.at(i), 2);
  }
  for(int i = 399; i < 500; i++){
    residual += std::pow(subtracted.at(i), 2);
  }

  // Calculte the standard error for the 200 points used to estimate the background
  
  residual = std::sqrt(residual)/200;

  // ****************************************************************************************************************************
  // We can't calculate the residual of the data in the signal area since we don't have and actual background so we will scale
  // the error assuming is fairly consistent over the full signal spectrum.
  //
  //      Spectrum 0:500 bins
  //      Background from 0:200 bins
  //      Scale the background by (500 - 200)/200 = 1.5
  //
  //  ****************************************************************************************************************************

  residual *= 1.5;
  std::cout << "Standard Error PDP: " << residual << std::endl;
						      
  // ****************************************************************************************************************************
  //
  // dA = x*y SQRT( (dx/x)^2 + (dy/y)^2 ) = x*y SQRT( 0 + (dy/y)^2 ) = x*y*(dy/y) = x*dy
  //
  // The total area is the sum over the frequency spacing times the dy.  A_total = N*delta_x*dy. Twice the modulation gives the
  // total frequency spectrum, N*delta_x, so dA = 2 x modulation x residual
  //
  // ****************************************************************************************************************************

  area_error = 2*modulation*residual;  

  canvas->cd();
  TMultiGraph *mgraph_pdp = new TMultiGraph();

  TGraph *pdp_avg  = new TGraph(steps, nmr->entry.at(0)->frequency.data(), pdp_average.data());
  pdp_avg->SetMarkerStyle(kFullDotMedium);
  pdp_avg->SetLineColor(kGreen+2);
  pdp_avg->SetMarkerColor(kGreen+2);
  pdp_avg->SetLineWidth(2);
  pdp_avg->Draw();
  
  TGraph *pdp_sub  = new TGraph(steps, nmr->entry.at(0)->frequency.data(), subtracted.data());
  pdp_sub->SetMarkerStyle(kFullDotMedium);
  pdp_sub->SetLineColor(kMagenta+1);
  pdp_sub->SetMarkerColor(kMagenta+1);
  pdp_sub->SetLineWidth(2);
  pdp_sub->Draw("A");
  
  mgraph_pdp->GetXaxis()->SetTitleSize(0.025);
  mgraph_pdp->GetYaxis()->SetTitleSize(0.025);
  mgraph_pdp->GetXaxis()->SetLabelSize(0.025);
  mgraph_pdp->GetYaxis()->SetLabelSize(0.025);
  mgraph_pdp->GetXaxis()->SetTitle("Frequency (MHz)");
  mgraph_pdp->GetYaxis()->SetTitle("Signal Amplitude");

  mgraph_pdp->Add(pdp_avg);
  mgraph_pdp->Add(pdp_sub);
  mgraph_pdp->Draw("ac");
 
  canvas->SaveAs(Form("output/pdp_comparison_data_set_%d.C", nmr->kDataSet));

  std::cout << "PDP Computed integral: " << pdp_sub->Integral() << " +- " << area_error << std::endl;

  nmr->PeakFinder(subtracted, "pdp");

  for(int k = 180; k < 320; k++) std::cout << "@@@@ " <<  k << " " << subtracted.at(k) << std::endl; 
  
  // Begin analysis of lanl data

  fit.clear();
  subtracted.clear();

  event_number = 0;
  steps = 0;
  modulation = 0;
  central_frequency = 0;
  area_error = 0;
  
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
    
    for(unsigned int i = 0; i < (unsigned int)(lanl->entry.size()); i++)
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

	//	canvas->SaveAs(Form("lanl_histogram%d.C", i));
	canvas->Clear();
	
	delete lal;
      }
  }

  lanl_average = lanl->ComputeSignalAverage(lanl->entry, "lanl");

  fit = lanl->ComputeBackgroundSignal(lanl_average,
				      lanl->entry.at(0)->frequency,
				      0.0, 5.0, 450.0, 500.0, "lanl");
  
  index = 0;
  offset = 0;
  
  for(std::vector <double>::iterator it = fit.begin(); it != fit.end(); it++){
    index = std::distance(fit.begin(), it);
    subtracted.push_back(lanl->entry.at(0)->data.at(index) - (*it) - offset);
  }

  if(std::abs(subtracted.front()) > 1e-6){
    offset = subtracted.front();
    for(std::vector <double>::iterator it = fit.begin(); it != fit.end(); it++){
      index = std::distance(fit.begin(), it);
      subtracted.at(index) = (lanl->entry.at(0)->data.at(index) - (*it) - offset);
    }
  }

  canvas->cd();
  TMultiGraph *mgraph_lanl = new TMultiGraph();

  lanl_average = lanl->ComputeSignalAverage(lanl->entry, "lanl");
  TGraph *lanl_avg  = new TGraph(steps, lanl->entry.at(0)->frequency.data(), lanl->entry.at(0)->data.data());
  lanl_avg->SetMarkerStyle(kFullDotMedium);
  lanl_avg->SetMarkerColor(kRed+2);
  lanl_avg->SetLineColor(kRed+2);
  lanl_avg->SetLineWidth(2);
  lanl_avg->Draw();

  
  TGraph *lanl_sub  = new TGraph(steps, lanl->entry.at(0)->frequency.data(), subtracted.data());
  lanl_sub->SetMarkerStyle(kFullDotMedium);
  lanl_sub->SetMarkerColor(kBlue+2);
  lanl_sub->SetLineColor(kBlue+2);
  lanl_sub->SetLineWidth(2);
  lanl_sub->Draw();
  
  mgraph_lanl->Add(lanl_avg);
  mgraph_lanl->Add(lanl_sub);
  mgraph_lanl->Draw("ac");
  mgraph_lanl->GetXaxis()->SetTitleSize(0.025);
  mgraph_lanl->GetYaxis()->SetTitleSize(0.025);
  mgraph_lanl->GetXaxis()->SetLabelSize(0.025);
  mgraph_lanl->GetYaxis()->SetLabelSize(0.025);
  mgraph_lanl->GetXaxis()->SetTitle("Frequency (MHz)");
  mgraph_lanl->GetYaxis()->SetTitle("Signal Amplitude");
  
  canvas->SaveAs(Form("output/lanl_comparison_data_set_%d.C", nmr->kDataSet));

  residual = 0;

  // Get residuals from subtracted only in the region we used to fit the background. These values should all be zero
  // if our fit was perfect.
  
  for(int i = 0; i < 100; i++){
    residual += std::pow(subtracted.at(i), 2);
  }
  for(int i = 399; i < 500; i++){
    residual += std::pow(subtracted.at(i), 2);
  }

  // Calculte the standard error for the 200 points used to estimate the background
  
  residual = std::sqrt(residual)/200;

  // ****************************************************************************************************************************
  // We can't calculate the residual of the data in the signal area since we don't have and actual background so we will scale
  // the error assuming is fairly consistent over the full signal spectrum.
  //
  //      Spectrum 0:500 bins
  //      Background from 0:200 bins
  //      Scale the background by (500 - 200)/200 = 1.5
  //
  //  ****************************************************************************************************************************

  residual *= 1.5;
  std::cout << "Standard Error LANL: " << residual << std::endl;
    
  area_error = 2*modulation*residual;  
  
  std::cout << "LANL Computed integral: " << lanl_sub->Integral() << " +- " << area_error << std::endl;

  lanl->PeakFinder(subtracted, "lanl");

  for(int k = 180; k < 320; k++) std::cout << "@@@@ " <<  k << " " << subtracted.at(k) << std::endl; 
  
  canvas->Clear();

  canvas->cd();

  TMultiGraph *mgraph = new TMultiGraph();
  
  mgraph->Add(pdp_sub);
  pdp_sub->SetMarkerStyle(kFullDotMedium);
  pdp_sub->SetMarkerColor(kBlue+2);
  pdp_sub->SetLineColor(kBlue+2);
  pdp_sub->SetLineWidth(2);
  
  mgraph->Add(lanl_sub);
  lanl_sub->SetMarkerStyle(kFullDotMedium);
  lanl_sub->SetMarkerColor(kRed+2);
  lanl_sub->SetLineColor(kRed+2);
  lanl_sub->SetLineWidth(2);
  
  mgraph->GetXaxis()->SetTitleSize(0.025);
  mgraph->GetYaxis()->SetTitleSize(0.025);
  mgraph->GetXaxis()->SetLabelSize(0.025);
  mgraph->GetYaxis()->SetLabelSize(0.025);
  mgraph->GetXaxis()->SetTitle("Frequency (MHz)");
  mgraph->GetYaxis()->SetTitle("Signal Amplitude");
  mgraph->Draw("ac");

  canvas->SaveAs(Form("output/system_comparison_data_set_%d.C", nmr->kDataSet));
  
  return(0);
}
