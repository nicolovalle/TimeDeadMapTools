#include <string>
#include <fstream>
#include <iostream>
#include <time.h>
#include <algorithm>

#include <map>
#include <vector>

#include <TBufferJSON.h>
#include <TH1.h>
#include <TH2.h>
#include <THnSparse.h>
#include <TFile.h>
#include <TTree.h>
#include <TSystem.h>
#include <TColor.h>
#include <TCanvas.h>
#include <TVector2.h>
#include <TH2Poly.h>
#include <TStyle.h>
#include <TGraph.h>
#include <TMath.h>
#include <nlohmann/json.hpp>

#pragma cling add_include_path(".")
#include "Logger.h"


#include "DataFormatsITSMFT/TimeDeadMap.h"
#include "DataFormatsITSMFT/NoiseMap.h"

using namespace TMath;


/// ________________________________________________________________________________________________________
/// settings
TString InputFile = "dmap.root";  // can be changed as argument of the macro
TString logfilename = "DeadMapTREE.log";
TString treefile = "DeadMapTREE.root";
std::string ccdbHost = "http://alice-ccdb.cern.ch"; // for RCT and CTP time stamps
/// _____________________________________________________


// global variables handled by the functions
std::map<Long64_t, std::vector<UShort_t>> MAP;  // MAP[orbit] = vector of dead chip IDs
std::vector<ULong_t> MAPKeys;
Int_t nkeys;
std::vector<Int_t> MAPNwords;
std::vector<UShort_t> SMAP;
std::map<TString, TString> QAcheck;
Long64_t runstart = -1, runstop = -1;
std::string mapver = "0";
int isdefault = false;
Logger QALOG;

float LHCOrbitNS = 88924.6; // o2::constants::lhc::LHCOrbitNS


std::vector<uint16_t> expandvector(std::vector<uint16_t> words, int version);
bool fillmap(TString fname, int map_sampling);


  
//////////////// _ MAIN _ /////////////////////
void DeadMapTREE(TString FILENAME = InputFile, int runnumber = -1, TString outdir="./", int map_sampling=1){

  QALOG.open(outdir+logfilename);
  QAcheck.clear();

  QALOG<<"Reading file "<<FILENAME<<". Run "<<runnumber<<"\n";


  //// Preparing trees
  TFile *f = TFile::Open(treefile, "RECREATE");

  if (!f || f->IsZombie()) {
    QALOG << "ERROR: cannot create file " << treefile << "\n";
    return; // or throw
  }
  else{
    QALOG << "File "<<treefile<<" will be written \n";
  }
  
  f->cd();  // <-- ensure gDirectory is the file

  
  TTree *td = new TTree("t_dynamic","dead map");
  TTree *ts = new TTree("t_static","info and static map");
  td->SetDirectory(f);
  ts->SetDirectory(f);

  std::vector<std::string> fatals = {};

  ts->Branch("version", &mapver);
  ts->Branch("run", &runnumber);
  ts->Branch("rctstart", &runstart);
  ts->Branch("rctstop", &runstop);
  ts->Branch("isdefault", &isdefault);
  ts->Branch("static", &SMAP);
  ts->Branch("nkeys", &nkeys);
  ts->Branch("fatal", &fatals);

  std::vector<UShort_t> chips;
  Long64_t orbit;
  Int_t nwords;

  td->Branch("nwords", &nwords);
  td->Branch("key", &orbit);
  td->Branch("deadchips", &chips);

  //////////////////////////////////////////////////

  bool mapdecoded = fillmap(FILENAME, map_sampling); // fill both MAP and SMAP. False if vector is not filled.
  nkeys = MAPKeys.size();
  
  if (!mapdecoded){
    QAcheck["Map decoded"] = "FATAL";
  }
 
  // Reading RCT
  if (runnumber <= 0) {
    QALOG<<"WARNING empty runnumber provided. Checks on run duration will not be effective\n";
  }
  else{
    QALOG<<"Reading RCT info for run "<<runnumber<<" from "<<ccdbHost<<"\n";
    auto& cm = o2::ccdb::BasicCCDBManager::instance();
    cm.setURL(ccdbHost);
    auto lims = cm.getRunDuration(runnumber);
    if (lims.first == 0 || lims.second == 0) { QALOG<<"ERROR, failed to fetch run info from RCT";}
    else { runstart = (Long64_t)lims.first;  runstop = (Long64_t)lims.second;}
  }
    
  QALOG<<"Writing decoded map to "<<treefile<<" ...\n";

 
  Long64_t llast=-1;
  int keycounter = 0;
  for (auto it = MAP.begin(); it!=MAP.end(); ++it) {
    orbit = it->first;
    chips = it->second;
    nwords = MAPNwords[keycounter];
    if ((long)orbit < llast){
      QALOG<<"FATAL PROBLEM in key orderding: "<<orbit<<" < "<<llast<<"\n";
      QAcheck["Key ordering"] = "FATAL";
      break;
    }
    llast = (Long64_t)orbit;
    td->Fill();
    if (keycounter % 1000 == 0) std::cout<<"written key "<<keycounter<<std::endl;
    keycounter++;
  }
  
  for (const auto& [key, val] : QAcheck){
    if (val == "FATAL") fatals.push_back(std::string(key.Data()));
  }

  // Just one entry in this tree
  ts->Fill();

  f->cd();
  ts->Write();
  td->Write();
  f->Close();

  
  QALOG<<"...done\n";


  exit(0);
}




  
std::vector<uint16_t> expandvector(std::vector<uint16_t> words, std::string version){

  std::vector<uint16_t> elementlist{};
  if (version == "2" || version == "3" || version == "4"){

    uint16_t firstel = 9999, lastel = 9999;

    for (long unsigned int i=0; i<words.size(); i++){

      uint16_t w = words[i];

      if (w & 0x8000){
	firstel = w & (0x7FFF);
	lastel = words[i+1];
	i++;
      }
      else {
	firstel = lastel = w;
      }

      
      for (uint16_t ic = firstel; ic <= lastel; ic++){
	elementlist.push_back(ic);
      }
    } // end of loop over words
  } // end of map version condition

  else{
    QALOG<<"FATAL: map version not recognized, returning empty vector.\n";
    QAcheck["MAP version"] = "FATAL";
  }
   
  return elementlist;
}
     

      


bool fillmap(TString fname, int map_sampling){

  MAP.clear();
  MAPKeys.clear();
  SMAP.clear();
  MAPNwords.clear();

  

  TFile *f = new TFile(fname);
  o2::itsmft::TimeDeadMap* obj = nullptr;
  
  f->GetObject("ccdb_object",obj);
  f->Close();

  if (!obj){
    QAcheck["ROOT file corrupted"] = "FATAL";
    return false;
  }
  

  QALOG<<"Reading map info\n";
  mapver = obj->getMapVersion();
  QALOG<<"Map version = "<<mapver<<"\n";
  isdefault = obj->isDefault();
  QALOG<<"Is default object = "<<isdefault<<"\n";
  QALOG<<"Evolving map size = "<<obj->getEvolvingMapSize()<<"\n";

  std::vector<uint16_t> StaticMap;
  obj->getStaticMap(StaticMap);

  QALOG<<"Static map imported\n";
  QALOG<<"Static map size = "<<StaticMap.size()<<"\n";

  if ((bool)isdefault){
    QAcheck["Default object"] = "FATAL";
    return false;
  }


  SMAP = expandvector(StaticMap,mapver);


  MAPKeys = obj->getEvolvingMapKeys();

  //if (MAPKeys.size() < 1){
  //  return false;
  //}
 
  for (int i=0; i<MAPKeys.size(); i+=map_sampling){

    if (i%1000==0){
      cout<<"step "<<i/1000<<"k"<<endl;
    }
       
    Long64_t OO = MAPKeys[i];
    std::vector<uint16_t> MapAtOrbit;
    obj->getMapAtOrbit(OO, MapAtOrbit);
    MAPNwords.push_back(MapAtOrbit.size());
    MAP[OO] = expandvector(MapAtOrbit,mapver);    
  }

  QALOG<<"... done importing maps for every orbit.\n";

  return true;
}



 
