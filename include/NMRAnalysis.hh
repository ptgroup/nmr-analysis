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

// ROOT libraries

#include "TString.h"
#include "TApplication.h"

// BOOST c++ libraries

#include <boost/any.hpp>
#include <boost/algorithm/string.hpp>

// Custom libraries

class NMRAnalysis : private std::fstream {

private:

  struct run {
    int event;

    std::vector <float> data;          // Amplitude of the NMR signal at a given frequncy
    std::vector <float> frequency;     // Frequencies making up the total frequency range
    
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
  
public:

  std::fstream config_file;
  std::fstream data_file;

  bool bConfigFileLoaded;
  bool bDataFileLoaded;
  bool bFilePrefixSet;
  bool bGraphicsShow;

  std::string fFilePrefix;
  
  std::vector <run *> entry;
  std::vector <data_type *> config_dict;
  std::vector <config_data *> configuration;
  std::vector< std::string > SplitVec;
  std::vector< std::string > SplitVecData;

  
  NMRAnalysis();
  ~NMRAnalysis();
  
  int OpenFiles();

  void ReadFiles();
  void GetOptions(char **);
  void ReadConfigurationMap();
  void is_string(int, const char *);
  void InitGraphicsEngine(int, char** );
  void RunGraphicsEngine();

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
  std::cout << __FUNCTION__ << ":: Failed to find value at index:: " << key << ":" << index << " Exiting." << std::endl;
  exit(1);
}
#endif
