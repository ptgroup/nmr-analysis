#ifndef nmr_analysis_cc
#define nmr_analysis_cc

#include "NMRAnalysis.hh"

NMRAnalysis::NMRAnalysis(){
  bConfigFileLoaded = false;
  bDataFileLoaded = false;
  bFilePrefixSet = false;
}

NMRAnalysis::~NMRAnalysis(){
  config_file.close();
  data_file.close();
}

void NMRAnalysis::ReadConfigurationMap()
{
  char *token;
  
  std::string line;
  std::string name;
  std::string type;
  
  std::fstream map_config;

  map_config.open("config/config_map.cfg", std::fstream::in);
  
  if( !(map_config.good()) ){
    std::cerr << __FUNCTION__ << " >> Error opening file: " << std::endl;
    exit(1);
  }

  while(!map_config.eof()){
    std::getline(map_config, line);
    token = new char[line.size() + 1];
    strcpy(token, line.c_str());
    config_dict.push_back(new data_type);

    token = strtok(token, " ,\t");
    (config_dict.back())->setting = token;
    token = strtok(NULL, " ,\t");
    (config_dict.back())->value = token;
    token = strtok(NULL, " ,\t");
  }
  return;
}

std::pair<bool, const char *> NMRAnalysis::FindSetting(const char *param){
  
  for(std::vector <data_type *>::iterator it = config_dict.begin(); it != config_dict.end(); ++it){
    if(strcmp((*it)->setting, param) == 0){
      return (std::pair<bool, const char *>(true, (*it)->value));
    }
  }
  return std::pair<bool, const char *>(false, "Failure");  
}

int NMRAnalysis::OpenFiles()
{

  std::string filename;

  if(!bFilePrefixSet){
    std::cerr << "No file prefix given. Please use --file-prefix <prefix> as a command-line argument." << std::endl;
    exit(1);
  }
  
  std::cout << "Reading files with prefix: " << fFilePrefix << std::endl;
 
  filename = Form("data/%s.csv", fFilePrefix.c_str());
  std::cout << "Reading ..... " << filename << std::endl;
  
  config_file.open(filename, std::fstream::in);
 
  if( !(config_file.good()) ){
    std::cerr << __FUNCTION__ << " >> Error opening file: " << std::endl;
    bConfigFileLoaded = false;
    exit(1);
  }
  else
    bConfigFileLoaded = true;
  

  filename = Form("data/%s-PolySignal.csv", fFilePrefix.c_str());
  data_file.open(filename, std::fstream::in);

  if( !(data_file.good()) ){
    std::cerr << __FUNCTION__ << " >> Error opening file: " << std::endl;
    bDataFileLoaded = false;
    exit(1);
  }
  else
    bDataFileLoaded = true;

  return(0);
}

void NMRAnalysis::ReadFiles(){
  this->ReadConfigurationMap();

  // Read in data file

  std::string line;

  char *token;

  while(!data_file.eof()){
    std::getline(data_file, line);
    token = new char[line.size() + 1];
    strcpy(token, line.c_str());
    token = strtok(token, " ,\t");

    entry.push_back(new run); 
    entry.back()->event = atoi(token);

    token = strtok(NULL, " ,\t");
    while(token){
      if(!token){
	entry.back()->data.push_back(0);
	break;
      }
      entry.back()->data.push_back(atof(token));
      token = strtok(NULL, " ,\t");
    }
  }

  // Read in cofiguration file
  
  bool bFirstLineRead = false;
  std::pair <bool, const char *> out;  

  while(!config_file.eof()){
    std::getline(config_file, line);
    token = new char[line.size() + 1];
    strcpy(token, line.c_str());
    token = strtok(token, ",");
    if(!bFirstLineRead){
      while(token){
	configuration.push_back(new config_data);
	configuration.back()->setting = token;
	token = strtok(NULL, ",");
      }
    }
    else{ 
      // This assumes you have a csv file that follows the structure given in the first line of the configuration file.
      for(std::vector<config_data *>::iterator it = configuration.begin(); it != configuration.end(); ++it){
	if(!token){
	  (*it)->values.push_back(new char *);
	}
	else{
	  out = FindSetting((*it)->setting);	  
	  if(strcmp(out.second, "int") == 0){
	    (*it)->values.push_back(atoi(token));
	  }
	  else if(strcmp(out.second, "float") == 0){
	    (*it)->values.push_back(atof(token));
	  }
	  else if(strcmp(out.second, "string") == 0){
	    (*it)->values.push_back(std::string(token));
	  }
	  else{
	    (*it)->values.push_back(token);
	  }
	}
	token = strtok(NULL, ",");
      }
    }
    bFirstLineRead = true;
  }
  
  return;
  
}

void NMRAnalysis::is_string(int value, const char *key)
{

  return;
  
}

void NMRAnalysis::InitGraphicsEngine(int Argc, char **Argv)
{
  std::cout << "<<<< Initialize Graphics Engine." << std::endl;
  app = new TApplication("App", &Argc, Argv);
}

void NMRAnalysis::RunGraphicsEngine()
{
  std::cout << "<<<< Running Graphics Engine." << std::endl;
  app->Run();
}

void NMRAnalysis::GetOptions(char **options){
  int i = 0;
  
  std::string flag;
  
  while(options[i] != NULL){
    flag = options[i];
    
    if(flag.compare("--file-prefix") == 0){
      std::string opt(options[i+1]);
      fFilePrefix = opt;
      flag.clear();
      std::cout << "Setting file prefix: \t" 
		<< fFilePrefix 
		<< std::endl;
      bFilePrefixSet = true;
    }
    if(flag.compare("--graphics") == 0){
      flag.clear();
      bGraphicsShow = true;
    }
    
    if(flag.compare("--help") == 0){
      printf("\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n");
      printf("Usage: ./nmranalysis <options>\n");
      printf("\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n");
      exit(0);
    }
    i++;
  }
  return;
}

#endif
