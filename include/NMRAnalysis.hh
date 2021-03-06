#ifndef nmr_analysis_hh
#define nmr_analysis_hh

// Standard c++ libraries

#include <cmath>
#include <fstream>
#include <ostream>
#include <iostream>
#include <cstdlib>
#include <vector>
#include <string>
#include <cstring>
#include <map>
#include <cctype>
#include <typeinfo>
#include <iomanip>

// ROOT libraries

#include "TString.h"
#include "TApplication.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TMultiGraph.h"
#include "TMatrixD.h"
#include "TH1.h"
#include "TSpectrum.h"

// BOOST c++ libraries

#include <boost/any.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>

// Custom libraries

class NMRAnalysis : private std::fstream {

private:

  struct run {
    int event;

    std::vector <double> data;          // Amplitude of the NMR signal at a given frequncy
    std::vector <double> frequency;     // Frequencies making up the total frequency range
    
  };

  struct  data_type {
    const char *setting;
    const char *value;     // Frequencies making up the total frequency range
    
  };

  struct config_data {
    const char * setting;          
    std::vector <boost::any> values;     // Frequencies making up the total frequency range
  };

  TApplication *app;

  std::vector <double> momentum;
  
public:

  std::fstream config_file;
  std::fstream data_file;
  std::fstream background_file;

  bool bConfigFileLoaded;
  bool bDataFileLoaded;
  bool bBackgroundFileLoaded;
  bool bFilePrefixSet;
  bool bGraphicsShow;
  bool bLoadNMRFile;
  bool bExternNMRAverageSet;
  bool bExternLANLAverageSet;
  bool bWriteOutputFile;

  int kDataSet;
  
  double kScaleFactor;
  double kExternLANLAverage;
  double kExternNMRAverage;
  double kExternLANLError;
  double kExternNMRError;
  
  std::string fFilePrefix;
  std::string fMapFile;
  std::string fOutputFile;
  
  std::vector <run *> entry;
  std::vector <run *> background;
  std::vector <data_type *> config_dict;
  std::vector <config_data *> configuration;
  std::vector <std::string> SplitVec;
  std::vector <std::string> SplitVecData;
  std::vector <double> param;
  
  NMRAnalysis();
  ~NMRAnalysis();
  
  int OpenSettingsFile();
  int OpenSettingsFile(const char *);
  int OpenDataFile();
  int OpenDataFile(const char *);
  int OpenBackgroundFile(const char *);
  
  void ReadDataFile();
  void ReadBackgroundFile();
  void ReadNMRFiles();
  void GetOptions(char **);
  void ReadConfigurationMap();
  void ReadConfigurationMap(const char *);
  void is_string(int, const char *);
  void InitGraphicsEngine(int, char** );
  void RunGraphicsEngine();
  void PrintData();
  void Clear();
  void Sort(std::vector <boost::filesystem::path> &, std::vector <boost::filesystem::path> &);

  std::vector <double> ScaleData(double, std::vector <double>);
  std::vector <double> PolynomialRegression(std::vector <double>, std::vector <double>, const char *);
  std::vector <double> ComputeBackgroundSignal(std::vector <double>, std::vector <double>, double, double, double, double, const char *);
  std::vector <double> ComputeSignalAverage(std::vector <run *>, const char *);
  std::vector <double> ComputeSignalError(std::vector <run *>, const char *);
  std::vector <double> GradientDescent(std::vector <double>, std::vector <double>, std::vector <double>, double);
  std::vector <double> PeakFinder(std::vector <double>, std::string type);
  
  std::vector <boost::filesystem::path> GetFileList(const boost::filesystem::path&, const std::string);
  
  std::pair<bool, const char *> FindSetting(const char *);
  
  template <typename T> T GetValue(const char *, int);
  template <typename T> std::vector<T> Get(const char *);
  
};

template <typename T> std::vector<T> NMRAnalysis::Get(const char *key){
  std::vector<T> v;

  for(std::vector <config_data *>::iterator it = configuration.begin(); it != configuration.end(); ++it){
    if(strcmp((*it)->setting, key) == 0){
      for(std::vector <boost::any>::iterator itv = (*it)->values.begin(); itv != (*it)->values.end(); ++itv){
	v.push_back(boost::any_cast<T>(*itv));
      }
    }
  }
  std::cout << __FUNCTION__ << ":: Failed to find setting in function. Returning empty vector." << std::endl;
  return v;
}
    
template <typename T> T NMRAnalysis::GetValue(const char *key, int index){
  for(std::vector <config_data *>::iterator it = configuration.begin(); it != configuration.end(); ++it){
    if(strcmp((*it)->setting, key) == 0){
      return boost::any_cast<T>((*it)->values.at(index));
    }
  }
  std::cout << __FUNCTION__ << ":: Failed to find value at index:: " << key << ": " << index << ". Exiting." << std::endl;
  exit(1);
}
#endif
