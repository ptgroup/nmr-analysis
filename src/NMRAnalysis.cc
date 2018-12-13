#ifndef nmr_analysis_cc
#define nmr_analysis_cc

#include "NMRAnalysis.hh"

NMRAnalysis::NMRAnalysis(){
  bConfigFileLoaded = false;
  bDataFileLoaded = false;
  bFilePrefixSet = false;
  bExternLANLAverageSet = false;
  bExternNMRAverageSet = false;
  bWriteOutputFile = false;
  
  fMapFile = "config_map.cfg";
  fOutputFile = "output.dat";
  
  kScaleFactor = 1.0;               // Default scaling of data. This is used primarily to flip the sign.
}

NMRAnalysis::~NMRAnalysis(){
  config_file.close();
  data_file.close();
}

std::vector<boost::filesystem::path> NMRAnalysis::GetFileList(const boost::filesystem::path& root, const std::string ext)
{

  std::vector <boost::filesystem::path> file_list;
  
  if(!boost::filesystem::exists(root) || !boost::filesystem::is_directory(root)){
    std::cerr << "Filesystems failed. " << std::endl;
    exit(1);
  }

  boost::filesystem::recursive_directory_iterator it(root);
  boost::filesystem::recursive_directory_iterator endit;

  while(it != endit)
    {
      if(boost::filesystem::is_regular_file(*it) && it->path().extension() == ext) file_list.push_back(it->path().filename());
      ++it;

    }
  return file_list;
}

void NMRAnalysis::Sort(std::vector <boost::filesystem::path> &settings, std::vector <boost::filesystem::path> &data)
{
  boost::filesystem::path temp;

  std::string substr;
  std::size_t npos;

  unsigned int i = 0;
  unsigned int j = 0;
  
  for(std::vector <boost::filesystem::path>::iterator it = settings.begin(); it != settings.end(); it++){
   npos = (*it).string().find(".csv");
    substr = (*it).string().substr(0, npos);
    i = std::distance(settings.begin(), it);
    
    for(std::vector <boost::filesystem::path>::iterator jt = data.begin(); jt != data.end(); jt++){
      j = std::distance(data.begin(), jt);
      if((*jt).string().compare(0, npos, substr) == 0){
	temp = data.at(i);
	data.at(i) = (*jt);
	data.at(j) = temp;
      }
    }
  }
}

void NMRAnalysis::Clear()
{

  entry.clear();
  config_dict.clear();
  configuration.clear();
  SplitVec.clear();
  SplitVecData.clear();
}

void NMRAnalysis::ReadConfigurationMap()
{
  char *token;
  
  std::string line;
  std::string name;
  std::string type;
  
  std::fstream map_config;

  map_config.open(Form("config/%s", fMapFile.c_str()), std::fstream::in);
  
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

std::vector <double> NMRAnalysis::ComputeSignalAverage(std::vector<run *> signal)
{
  int i = 0;
  int j = 0;
  
  std::vector <double> average;

  average.resize(signal.front()->data.size());
  
  for(std::vector <run *>::iterator it = signal.begin(); it != signal.end(); it++){
    i = std::distance(signal.begin(), it);
    for(std::vector <double>::iterator jt = signal.at(i)->data.begin(); jt != signal.at(i)->data.end(); jt++){
      j = std::distance(signal.at(i)->data.begin(), jt);
      average.at(j) += *jt;
    }
  }

  for(std::vector <double>::iterator it = average.begin(); it != average.end(); it++){
    i = std::distance(average.begin(), it);
    average.at(i) /= signal.size();
  }

  return average;
}

TGraph *NMRAnalysis::ComputeBackgroundSignal(std::vector <double> signal, std::vector <double> frequency, double x0, double x1, double x2, double x3){
  std::vector <double> y;
  std::vector <double> x;
  
  // *****************************************************************************
  //
  // Create a vector containing only the the parts of the signal you want to fit
  //
  // *****************************************************************************
  
  for(int i = x0; i < x1; i++){
    y.push_back(signal.at(i));
    x.push_back(frequency.at(i));
  }
  for(int i = x2; i < x3; i++){
    y.push_back(signal.at(i));
    x.push_back(frequency.at(i));
  }

  TCanvas *c = new TCanvas("c", "c", 5);
  c->cd();
  
  TGraph *background = new TGraph(x.size(), x.data(), y.data());
  background->Draw("ap");
  c->SaveAs("background.C");

  return background;
  
}

std::vector <double> NMRAnalysis::GradientDescent(std::vector <double> x, std::vector <double> y, std::vector <double> _coeff, double gamma)
{
  int j = 0;
  
  const double N = y.size();
  
  std::vector <double> coeff(3);
  std::vector <double> gradient(3);

  for(int i = 0; i < 3; i++){
    for(std::vector <double>::iterator it = y.begin(); it != y.end(); it++){
      j = std::distance(y.begin(), it);
      gradient.at(i) += (-2/N)*std::pow(x.at(j),i)*(y.at(j)
						    - _coeff.at(0)*std::pow(x.at(j), 0)
						    - _coeff.at(1)*std::pow(x.at(j), 1)
						    - _coeff.at(2)*std::pow(x.at(j), 2));
    }
    coeff.at(i) = _coeff.at(i) - gamma*gradient.at(i) - 0.9*gamma*momentum.at(i);
    momentum.at(i) = coeff.at(i);
  }

  return coeff;
}

double NMRAnalysis::PolynomialRegression(std::vector <double> x, std::vector <double> y)
{
  int count = 0;
  
  double gamma;
  double training = 10;

  momentum.resize(3);
  momentum.assign(3, 0);
  
  std::vector <double> gd;
  gd.assign(3, 0);
  
  for(int i = 0; i < 1e8; i++){ 
    gamma = 0.002/(1 + i/training);
    gd = this->GradientDescent(x, y, gd, gamma);
    if((i % (int)1e5) == 0){
      std::cout << "<<<<< Parameters: " << i << " "
		<< gd.at(0) << " "
		<< gd.at(1) << " "
		<< gd.at(2) << std::endl;
    }
  }
  
  return 0;
}

void NMRAnalysis::ReadConfigurationMap(const char *mapfile)
{
  char *token;
  
  std::string line;
  std::string name;
  std::string type;
  
  std::fstream map_config;

  map_config.open(Form("config/%s", mapfile), std::fstream::in);
  
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
  map_config.close();
  
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

int NMRAnalysis::OpenSettingsFile()
{
  std::string filename;

  if(!bFilePrefixSet){
    std::cerr << "No file prefix given. Please use --file-prefix <prefix> as a command-line argument." << std::endl;
    exit(1);
  }
  
  std::cout << "opening files with prefix: " << fFilePrefix << std::endl;
 
  filename = Form("data/%s.csv", fFilePrefix.c_str());
  std::cout << "Opening ..... " << filename << std::endl;
  
  config_file.open(filename, std::fstream::in);
 
  if( !(config_file.good()) ){
    std::cerr << __FUNCTION__ << " >> Error opening file: " << std::endl;
    bConfigFileLoaded = false;
    exit(1);
  }
  else
    bConfigFileLoaded = true;

  std::cout << "Finished opening settings file." << std::endl;
  
  return(0);
}

int NMRAnalysis::OpenSettingsFile(const char *filename)
{
  bFilePrefixSet = true;

  std::cout << "opening settings file: " << filename << std::endl;
  
  config_file.open(Form("%s", filename), std::fstream::in);
 
  if( !(config_file.good()) ){
    std::cerr << __FUNCTION__ << " >> Error opening settings file." << std::endl;
    bConfigFileLoaded = false;
    exit(1);
  }
  else
    bConfigFileLoaded = true;

  std::cout << "Finished opening settings file." << std::endl;
  
  return(0);
}

int NMRAnalysis::OpenDataFile()
{
  std::string filename;

  if(!bFilePrefixSet){
    std::cerr << "No file prefix given. Please use --file-prefix <prefix> as a command-line argument." << std::endl;
    exit(1);
  }
  
  std::cout << "opening files with prefix: " << fFilePrefix << std::endl;
 
  filename = Form("data/%s-RawSignal.csv", fFilePrefix.c_str());
  data_file.open(filename, std::fstream::in);
  std::cout << "Opening ..... " << filename << std::endl;
  
  if( !(data_file.good()) ){
    std::cerr << __FUNCTION__ << " >> Error opening file: " << std::endl;
    bDataFileLoaded = false;
    exit(1);
  }
  else
    bDataFileLoaded = true;

  std::cout << "Finished opening data file." << std::endl;
  
  return(0);
}

int NMRAnalysis::OpenDataFile(const char *filename)
{
  bFilePrefixSet = true;
  
  std::cout << "opening data file: " << filename << std::endl;
  data_file.open(Form("%s", filename), std::fstream::in);
  
  if( !(data_file.good()) ){
    std::cerr << __FUNCTION__ << " >> Error opening data file: " << std::endl;
    bDataFileLoaded = false;
    exit(1);
  }
  else
    bDataFileLoaded = true;

  std::cout << "Finished opening data file." << std::endl;
  
  return(0);
}

int NMRAnalysis::OpenBackgroundFile(const char *filename)
{
  bFilePrefixSet = true;
  
  std::cout << "opening background file: " << filename << std::endl;
  background_file.open(Form("%s", filename), std::fstream::in);
  
  if( !(background_file.good()) ){
    std::cerr << __FUNCTION__ << " >> Error opening background file: " << std::endl;
    bBackgroundFileLoaded = false;
    exit(1);
  }
  else
    bBackgroundFileLoaded = true;

  std::cout << "Finished opening background file." << std::endl;
  
  return(0);
}

void NMRAnalysis::PrintData()
{
  for(unsigned int i = 0; i < entry.size(); i++){
    for(unsigned int j = 0; j < (entry.at(i)->data.size()); j++){
      std::cout << __FUNCTION__ << ":: [" << i <<", " << j << "] " << entry.at(i)->data.at(j) << std::endl;
    }
  }
  return;  
}

void NMRAnalysis::ReadDataFile(){
  this->ReadConfigurationMap();
  entry.clear();
  
  // Read in data file
  std::cout << "Reading data file." << std::endl;
  std::string line;

  while(!data_file.eof()){
    std::getline(data_file, line);
    if((unsigned)strlen(line.c_str()) == 0) continue;
    
    boost::split(SplitVec, line, boost::is_any_of(",\t"));

    entry.push_back(new run);
    entry.back()->event = atoi(SplitVec.front().c_str());
    
    for(std::vector<std::string>::iterator it = SplitVec.begin()+1; it != SplitVec.end(); ++it){
      entry.back()->data.push_back(atof((*it).c_str()));
    }
  }
  data_file.close();
  
  return;
}

void NMRAnalysis::ReadNMRFiles(){

  entry.clear();
  configuration.clear();
  
  // Read in data file
  std::cout << "Reading data file." << std::endl;
  std::string line;

  while(!data_file.eof()){
    std::getline(data_file, line);
    if((unsigned)strlen(line.c_str()) == 0) continue;
    
    boost::split(SplitVec, line, boost::is_any_of(","));

    entry.push_back(new run);
    entry.back()->event = atoi(SplitVec.front().c_str());
    
    for(std::vector<std::string>::iterator it = SplitVec.begin()+1; it != SplitVec.end(); ++it){
      entry.back()->data.push_back(atof((*it).c_str()));
      // std::cout << "data >>> " << entry.back()->data.back() << std::endl; 
    }
  }

  data_file.close();
  
  // Read in configuration file

  std::cout << "Reading configuration file." << std::endl;
  
  bool bFirstLineRead = false;
  std::pair <bool, const char *> out;  
  unsigned int iter = 0;

  
  while(!config_file.eof()){
    std::getline(config_file, line);
    if((unsigned)strlen(line.c_str()) == 0) continue;

    SplitVecData.clear();
    
    if(bFirstLineRead)
      boost::split(SplitVecData, line, boost::is_any_of(","));
    
    if(!bFirstLineRead){
      boost::split(SplitVec, line, boost::is_any_of(","));
      for(std::vector<std::string>::iterator it = SplitVec.begin(); it != SplitVec.end(); ++it){
	configuration.push_back(new config_data);
	configuration.back()->setting = (*it).c_str();
      }
      bFirstLineRead = true;
    }
    else{ 
      //      This assumes you have a csv file that follows the structure given in the first line of the configuration file.
	for(std::vector<config_data *>::iterator it = configuration.begin(); it != configuration.end(); ++it){
	  
	  if((unsigned)strlen(SplitVecData.at(iter).c_str()) == 0){
	   (*it)->values.push_back("NA");
	   }
	   else{
	     out = FindSetting((*it)->setting);
	   
	     if(strcmp(out.second, "int") == 0){
	       (*it)->values.push_back(atoi(SplitVecData.at(iter).c_str()));
	     }
	     else if(strcmp(out.second, "float") == 0){
	       (*it)->values.push_back(atof(SplitVecData.at(iter).c_str() ));
	     }
	     else if(strcmp(out.second, "string") == 0){
	       (*it)->values.push_back(SplitVecData.at(iter));
	     }
	     else{
	       (*it)->values.push_back(SplitVecData.at(iter));
	       continue;
	     }
	   }
	  ++iter;
	}
	iter = 0;
    }
  }
  config_file.close();
  
  return;
}

void NMRAnalysis::ReadBackgroundFile(){
  background.clear();
    
  // Read in background file
  std::cout << "Reading background file." << std::endl;
  std::string line;

  while(!background_file.eof()){
    std::getline(background_file, line);
    if((unsigned)strlen(line.c_str()) == 0) continue;

    SplitVec.clear();    
    boost::split(SplitVec, line, boost::is_any_of(","));

    background.push_back(new run);
    background.back()->event = atoi(SplitVec.front().c_str());
    
    for(std::vector<std::string>::iterator it = SplitVec.begin()+1; it != SplitVec.end(); ++it){
      background.back()->data.push_back(atof((*it).c_str()));
    }
  }
  background_file.close();
  
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
    if(flag.compare("--output-file") == 0){
      std::string opt(options[i+1]);
      fOutputFile = opt;
      flag.clear();
      std::cout << "Setting output file: \t" 
		<< fOutputFile 
		<< std::endl;
      bWriteOutputFile = true;
    }
    if(flag.compare("--scale-factor") == 0){
      std::string opt(options[i+1]);
      kScaleFactor = atof(opt.c_str());
      flag.clear();

    }
    if(flag.compare("--set-lanl-average") == 0){
      std::string opt(options[i+1]);
      std::vector <std::string> splitvec;

      boost::split(splitvec, opt, boost::is_any_of(":"));
		   
      kExternLANLAverage = atof(splitvec.front().c_str());
      kExternLANLError = atof(splitvec.back().c_str());
      bExternLANLAverageSet = true;
      flag.clear();
    }
    if(flag.compare("--set-nmr-average") == 0){
      std::string opt(options[i+1]);
      std::vector <std::string> splitvec;

      boost::split(splitvec, opt, boost::is_any_of(":"));
      
      kExternNMRAverage = atof(splitvec.front().c_str());
      kExternNMRError = atof(splitvec.back().c_str());
      bExternNMRAverageSet = true;
      flag.clear();
    }
    if(flag.compare("--data-set") == 0){
      std::string opt(options[i+1]);
      kDataSet = atof(opt.c_str());
      flag.clear();
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
