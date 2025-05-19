/*****************************/
/*                           */
/*  Quick optical simulator  */
/*                           */
/*  Author: Tatsuya Kikawa   */
/*                           */
/*****************************/

#include <iostream>
#include <iomanip>
#include <vector>
#include <sstream>
#include <random>
#include <getopt.h>
#include <sys/time.h>
#include <TApplication.h>
#include "Global.hh"
#include "Geometry.hh"
#include "Material.hh"
#include "Source.hh"
#include "Triangle.hh"
#include "Simulator.hh"

void PrintUsage(){
    std::cout << "Usage:" << std::endl;
    std::cout << " OptSim [input card file] [output root file]" << std::endl;
    std::cout << "Options:" << std::endl;
    std::cout << " -s : Set the seed of random number generator" << std::endl;
    std::cout << " -d : Display mode" << std::endl;
    exit(1);
}

int main(int argc, char *argv[])
{
  if(argc<3) {
    PrintUsage();
  }
  TApplication* app = new TApplication("app", 0, 0);

  struct timeval  start_time, end_time;
  gettimeofday(&start_time, NULL);

  int seed=-1;
  bool dmode=false;
  int c;
  while ((c = getopt(argc, argv, "s:d")) != -1) {
    switch(c){
    case 's':
      seed = atoi(optarg);
      break;
    case 'd':
      dmode = true;
      break;
    default:
      PrintUsage();
    }
  }

  std::string input;
  std::string output;
  int op=0;
  while (optind < argc){
    switch(op){
    case 0:
      input = argv[optind++];
      break;
    case 1:
      output = argv[optind++];
      break;
    default:
      std::cerr << "Too many arguments." << std::endl;
      PrintUsage();
    }
    op++;
  }
  if(input.length()==0 || output.length()==0)PrintUsage();
  
  std::cout<<"Reading input card file."<<std::endl;
  Config config;
  ReadInFile(input.c_str(), config);

  std::cout<<"Setting pseudorandom number generator."<<std::endl;
  std::mt19937 mt;//Mersenne twister generator
  if(seed!=-1){
    mt.seed(seed);
  }
  else{
    std::random_device rnd;//Non-deterministic random number generator
    mt.seed(rnd());
  }

  std::cout<<"Loading source data."<<std::endl;
  Source *src = new Source(mt,config);

  std::cout<<"Loading material data."<<std::endl;
  std::vector<Material*> materials;
  std::string matfile, type;
  int id;
  double index, attlen, scatlen;
  for(std::map<std::string, std::string>::iterator i = config["Material"].begin(); i != config["Material"].end(); i++){
    std::istringstream(i->first) >> id;
    std::istringstream ss(i->second);
    ss >> matfile >> type;
    if(type=="medium" || type=="converter"){
      ss >> index >> attlen >> scatlen;
    }
    else{
      index=1;
      attlen=0;
      scatlen=0;
    }
    Material *mat = new Material(mt,id, matfile, type, index, attlen, scatlen);
    materials.push_back(mat);
  }

  std::cout<<"Setting up simulator."<<std::endl;
  Simulator *sim = new Simulator(mt,config,output,src,materials);
  
  std::cout<<"Start simulation."<<std::endl;
  if(!dmode)sim->Run();//Normal run mode
  else sim->Display();//Display mode
  
  gettimeofday(&end_time, NULL);
  double time_diff = end_time.tv_sec - start_time.tv_sec + (double)(end_time.tv_usec - start_time.tv_usec)/1e6;
  std::cout << "CPU time: "<< std::setprecision(4) << time_diff <<" sec."<< std::endl;
  std::cout << "All finished. \\(^o^)/" << std::endl;

  app->Terminate();
  return 0;
}
