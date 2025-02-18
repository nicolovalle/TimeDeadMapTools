////
//// This macro checks the existence of ITS dead maps and performs basic sanity
//// checks on their content.
////
//// Usage:
////    .x SimpleDeadMapCheck.C("653421")
//// or:
////    .x SimpleDeadMapCheck.C("list.txt")
////
//// Errors and warnings out of the following checks are printed at the end:
////   1) By querying the start-of-run and end-of-run timestamps the same object must be retrieved
////   2) The map must not be empty or the default object
////   3) The first orbit of the map must be meaningful (not zero or negative)
////   4) The total orbit span of the map should be similar to the run duration
////   5) The maximum gap between consecutive map elements must be small
////
//// For suggestions, or to report errors: nicolo.valle@cern.ch
////

#include <string>
#include <fstream>
#include <iostream>
#include <algorithm>

#include <map>
#include <vector>

#include <TMath.h>

#include "DataFormatsITSMFT/TimeDeadMap.h"

#define myLOG std::cout<<"_________"


float LHCOrbitNS = 88924.6;
std::string ccdbHost = "http://alice-ccdb.cern.ch";
std::string detector = "ITS";

std::map<int, TString> PrintInfo;

void DownloadAndCheck(int run);

void SimpleDeadMapCheck(const TString& input){

  PrintInfo = std::map<int, TString>{};

  int run = -1;
  // if input is a run number...
  if (input.IsDigit()){
    run = input.Atoi();
    DownloadAndCheck(run);
  }


  // if input is a text file...
  else{
    std::ifstream file(input.Data());

    if (!file.is_open()){
      std::cerr<<"Error: Unable to open file "<<input<<std::endl;
      return;
    }
    std::string line;
    while (std::getline(file, line)){
      std::stringstream ss(line);
      if (ss >> run){
	DownloadAndCheck(run);
      }
    }
    file.close();
  }

  myLOG<<"=================\n";
  myLOG<<"==== SUMMARY ====\n";
  myLOG<<"=================\n";
  for (auto i : PrintInfo){
    myLOG<<"Run "<<i.first<<i.second<<"\n";
  }

  return;
}


void DownloadAndCheck(int run){

  PrintInfo[run] = TString();

  long runstart = -1, runstop = -1;
  bool queryError = false;
  
  if (run <= 0){
    myLOG<<"FATAL - invalid run number provided\n";
    queryError = true;
    return;
  }

  myLOG<<"INFO - Reading start/stop timestamps for run "<<run<<"\n";

  auto& cm = o2::ccdb::BasicCCDBManager::instance();
  cm.setURL(ccdbHost);
  auto lims = cm.getRunDuration(run);
  if (lims.first == 0 || lims.second == 0) {
    myLOG<<"ERROR - failed to fetch run info from RCT for run "<<run<<"\n";
    queryError = true;
    PrintInfo[run]+=" - ERROR: failed to fetch run info from ccdb";
    return;
  }

  runstart = (long)lims.first;
  runstop = (long)lims.second;

  myLOG<<"INFO - Checking deadmap object for run "<<run<<"\n";

  auto* obj0 = cm.getForTimeStamp<o2::itsmft::TimeDeadMap>(detector+"/Calib/TimeDeadMap",runstop);
  auto* obj1 = cm.getForTimeStamp<o2::itsmft::TimeDeadMap>(detector+"/Calib/TimeDeadMap",runstart);

  if (obj0->isDefault() && obj1->isDefault()){
    myLOG<<"ERROR - Object for run "<<run<<" is missing. Only default object has been found. Report to experts.\n";
    PrintInfo[run] +=" - ERROR: default object fetched";
    queryError = true;
  }
  
  else if (obj0->isDefault() || obj1->isDefault() || obj0->getEvolvingMapKeys() != obj1->getEvolvingMapKeys()){
    myLOG<<"ERROR - Start and Stop timestamps of run "<<run<<" result in different ccdb objects. Report to experts.\n";
    PrintInfo[run] +=" - ERROR: mismatch between queries at start and stop of the run";
    queryError = true;
  }

  if (obj0->getEvolvingMapKeys().size() < 1){
    myLOG<<"ERROR - The time-evolving map is empty\n";
    PrintInfo[run] +=" - ERROR: the time-evolving map is empty";
    queryError = true;
  }

  if (queryError){
    return;
  }

  std::vector<unsigned long> mapkeys = obj0->getEvolvingMapKeys();

  myLOG<<"INFO - map version: "<<obj0->getMapVersion()<<"\n";
  myLOG<<"INFO - number of orbits in the map: "<<mapkeys.size()<<"\n";


  long firstorbit = (long)mapkeys.front();
  long lastorbit = (long)mapkeys.back();

  double runduration = (runstop - runstart)/1000.;
  double mapduration = (lastorbit - firstorbit) * (LHCOrbitNS * 1.e-9);

  int nTooLargeGaps = 0;
  for (int ikey = 1; ikey< mapkeys.size(); ikey++){
    if ( (long)mapkeys.at(ikey) - mapkeys.at(ikey-1) > 330000){
      nTooLargeGaps++;
    }
  }

  if ((nTooLargeGaps > 0 && firstorbit > 0) || nTooLargeGaps > 1){
    myLOG<<"ERROR - There are "<<nTooLargeGaps<<" orbit gaps exceeding 330k orbits\n";
    PrintInfo[run] += Form(" - ERROR: orbit gap exceeds 330k orbits %d times",nTooLargeGaps);
  }


  if (firstorbit < 1){
    PrintInfo[run] += Form(" - WARNING: first orbit saved in the map is %ld",firstorbit);
    if (mapkeys.size() > 2){
      mapduration = (lastorbit - (long)mapkeys.at(1)) * (LHCOrbitNS * 1.e-9);
    }
  }

  if (TMath::Abs(runduration - mapduration) <= 5){
    PrintInfo[run] += Form(" - run duration: %.0f sec, map duration: %.0f sec: OK",runduration,mapduration);
  }
  else if (TMath::Abs(runduration - mapduration) <= 60){
    PrintInfo[run] += Form(" - WARNING: run duration: %.0f sec, map duration: %.0f sec.",runduration,mapduration);
  }
  else{
    PrintInfo[run] += Form(" - ERROR: run duration: %.0f sec, map duration: %.0f sec.",runduration,mapduration);
  }

    
  return;
  
}
