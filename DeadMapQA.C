#include <string>
#include <fstream>
#include <iostream>
#include <time.h>
#include <algorithm>

#include <map>
#include <vector>

#include <TBufferJSON.h>
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

#include "Logger.h"


#include "DataFormatsITSMFT/TimeDeadMap.h"
#include "DataFormatsITSMFT/NoiseMap.h"

using namespace TMath;


/// ________________________________________________________________________________________________________
/// settings
TString InputFile = "dmap.root";  // can be changed as argument of the macro
TString logfilename = "DeadMapQA.log";
double SecForWorse = 15;
bool ExitWhenFinish = true;
std::string ccdbHost = "http://alice-ccdb.cern.ch"; // for RCT and CTP time stamps
Long_t NominalGap = 32000;
Long_t UnanchorableThreshold = 330000;
const std::vector<std::vector<int>> Enabled{ // not in use yet
  {0,1,2,3,4,5,6,7,8,9,10,11}, // L0
  {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14}, // L1
  {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19}, // L2
  {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23}, // L3
  {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29}, // L4
  {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41}, // L5
  {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46} // L6
};

// __________________________________________________________________________________________________________


// global variables handled by the functions
std::map<unsigned long, std::vector<uint16_t>> MAP;
std::vector<unsigned long> MAPKeys;
std::vector<uint16_t> SMAP;
std::map<TString, TString> QAcheck;
long runstart = -1, mapstart = -1, runstop = -1, mapstop = -1;
Logger QALOG;



const int NStaves[7] = { 12, 16, 20, 24, 30, 42, 48 };
const int NZElementsInHalfStave[7] =  {9,9,9, 4, 4, 7, 7};
const int NSegmentsStave[7] = {1, 1, 1, 4, 4, 4, 4};
const int N_LANES_IB = 432;
const int N_LANES_ML = 864; // L3,4
const int N_LANES = 3816;
const int N_STAVES_IB = 12+16+20;
const int N_STAVES = 192;
int LaneToLayer[N_LANES]; // filled by "getlanecoordinates" when called
int LaneToStave[N_LANES]; // filled by "getlanecoordinates" when called
int LaneToStaveInLayer[N_LANES]; // filled by "getlanecoordinates" when called

float LHCOrbitNS = 88924.6; // o2::constants::lhc::LHCOrbitNS



void interpolatestave(double r1, double r2, double phi1, double phi2, int z, int nz, int a, int na, TVector2 *TVec);
void getlanecoordinates(int laneid, double *px, double *py);

uint16_t isFirstOfLane(uint16_t chipid);
uint16_t isLastOfLane(uint16_t chipid);
uint16_t StaveToLayer(uint16_t stv);
uint16_t FirstLaneOfStave(uint16_t stv);
uint16_t LastLaneOfStave(uint16_t stv);

std::vector<uint16_t> expandvector(std::vector<uint16_t> words, int version);

uint16_t ChipToLane(uint16_t chipid);

void fillmap(TString fname);

void RemoveAxis(TH2Poly *HP);
void GetTimeStamps(int runnumber, uint32_t orbit1, uint32_t orbit2);

void PrintAndExit(TString spec=""){

  QALOG<<"\nQA SUMMARY\n";

  bool IsQAMedium=false, IsQABad=false, IsQAFatal=false;
  for (auto cc : QAcheck){
    QALOG<<"QA CHECK - "<<cc.first<<": "<<cc.second<<"\n";
    if (cc.second == "MEDIUM") IsQAMedium = true;
    if (cc.second == "BAD") IsQABad = true;
    if (cc.second == "FATAL") IsQAFatal = true;
  }

  if (IsQAFatal){
    QALOG<<"FATAL - Global quality is FATAL\n";
  }
  else if (IsQABad){
    QALOG<<"ERROR - Global quality is BAD\n";
  }
  else if (IsQAMedium){
    QALOG<<"WARNING - Global quality is MEDIUM\n";
  }
  else{
    QALOG<<"INFO - Global quality is GOOD\n";
  }
  
  QALOG<<"Exiting the macro. "<<spec<<"\n";

  QALOG.close();
  
  exit(0);
}
  

void DeadMapQA(TString FILENAME = InputFile, int runnumber = -1, TString outdir="./"){

  QALOG.open(outdir+logfilename);
  QAcheck.clear();

  QALOG<<"Checking file "<<FILENAME<<". Run "<<runnumber<<"\n";

  TH2Poly *HMAP = new TH2Poly(); // evolving map
  TH2Poly *HSMAP = new TH2Poly(); // static part

  TH2Poly *WorstOB = new TH2Poly();
  TH2Poly *WorstIB = new TH2Poly();

  TH2Poly *LastMAP = new TH2Poly(); // last snapshot

  TCanvas *c2 = new TCanvas("QAsummary2","QAsummary2", 6300, 2960);
  c2->Divide(3,1);

  TCanvas *c1 = new TCanvas("QAsummary1","QAsummary1",2400,1800);
  c1->Divide(4,3);

  TCanvas *c3 = new TCanvas("QAsummary3","QAsummary3",4000,1500);

  gStyle->SetOptStat(0);
  gStyle->SetPalette(kBlackBody);
  TColor::InvertPalette();

  
  for (int i=0; i<N_LANES; i++){
    double *px = new double[4];
    double *py = new double[4];
    getlanecoordinates(i, px, py);
    HMAP->AddBin(4,px,py);
    HSMAP->AddBin(4,px,py);
    WorstOB->AddBin(4,px,py);
    WorstIB->AddBin(4,px,py);
    LastMAP->AddBin(4,px,py);
  }
  RemoveAxis(HMAP);
  RemoveAxis(HSMAP);
  RemoveAxis(WorstOB);
  RemoveAxis(WorstIB);
  RemoveAxis(LastMAP);


  QAcheck["Chip interval"] = "GOOD";
  fillmap(FILENAME); // fill both MAP and SMAP, checking them. Exit if map is empty or default.
  
  TH1F *hOrb = new TH1F("Orbits gap","Orbits gap",MAP.size()-1,1,MAP.size());
  TH1F *hEffOB = new TH1F("OB dead fraction","OB (blue) and IB (red) dead fraction",MAP.size()-1,1,MAP.size());
  TH1F *hEffIB = new TH1F("IB dead fraction","IB dead fraction",MAP.size()-1,1,MAP.size());
  TH1F *hTimeSpan = new TH1F("Time range","Time range",2,0,2);

  const Int_t hInbin = MAP.size()-1;
  Double_t hIbins[hInbin+1];
  for (int i=0; i < hInbin+1; i++) hIbins[i] = (Double_t)MAPKeys[i];
  TH2F *hStatusTimeIB = new TH2F("Status vs time IB",Form("Run %d - Status vs time IB;Orbit;lane ID",runnumber), hInbin, hIbins, N_LANES_IB, 0, N_LANES_IB);
  TH2F *hStatusTimeML = new TH2F("Status vs time ML",Form("Run %d - Status vs time ML;Orbit;lane ID",runnumber), hInbin, hIbins, N_LANES_ML, N_LANES_IB, N_LANES_IB+N_LANES_ML);
  TH2F *hStatusTimeOL = new TH2F("Status vs time OL",Form("Run %d - Status vs time OL;Orbit;lane ID",runnumber), hInbin, hIbins, N_LANES-N_LANES_IB-N_LANES_ML, N_LANES_IB+N_LANES_ML, N_LANES);

  TH1F *hStaveDeadTime = new TH1F("Stave dead time",Form("Run %d - Stave dead time;;dead time",runnumber), N_STAVES,0,N_STAVES);
  for (int i=0; i < N_LANES; i++) hStaveDeadTime->GetXaxis()->SetBinLabel(LaneToStave[i]+1, Form("#color[%d]{L%d_%d}",1,LaneToLayer[i],LaneToStaveInLayer[i]));
  
  
  HMAP->Clear();
  HSMAP->Clear();
  HMAP->SetTitle("Lane dead time"); HMAP->SetName("Lane dead time");
  HSMAP->SetTitle("Number of fully dead chips in the lane"); HSMAP->SetName("Number of fully dead chips in the lane");
  hOrb->GetXaxis()->SetTitle("step");
  hOrb->GetYaxis()->SetTitle("Delta(orbit) from previous step");
  WorstOB->Clear();
  WorstIB->Clear();
   
  
  Long_t firstorbit = 0;
  int countstep = 0;

  Long_t currentorbit = 0;
  Long_t previousorbit = 0;
  
  Long_t maxgap = 0;
  int ngap_overnominal = 0;

  Long_t unAnchorable = 0;

  const int NSteps = MAP.size();

  const int nn = N_LANES;

  double eff[nn], effstat[nn];
  for (int i =0 ; i<nn; i++) eff[i]=effstat[i]=0;

  for (auto c : SMAP){
    //cout<<"DEBUG Incrementing static map counter chip / lane "<<c<<"/"<<ChipToLane(c);
    effstat[ChipToLane(c)]+=1;
  }

  for (int i = 0; i<nn; i++) {
    //cout<<"DEBUG bin "<<i+1<<" content "<<effstat[i];
    HSMAP->SetBinContent(i+1,effstat[i]);
  }

  int nfullydeadIB = 0;
  int nwithfullydeadOB = 0;
  for (int i=0; i<nn; i++){
    if (effstat[i] > 0 && i < N_LANES_IB) nfullydeadIB++;
    if (effstat[i] > 0 && i >= N_LANES_IB) nwithfullydeadOB++;
  }
  QALOG<<"Static map: IB dead chips: "<<nfullydeadIB<<"\n";
  QALOG<<"Static map: OB lanes with at least one fully dead chip: "<<nwithfullydeadOB<<"\n";
  QAcheck["Fully dead IB"] = (nfullydeadIB == 0) ? "GOOD" : (nfullydeadIB < 4) ? "MEDIUM" : "BAD";
  QAcheck["Fully dead OB"] = (1.*nwithfullydeadOB/N_LANES) < 0.02 ? "GOOD" : "BAD";
  


  int worstOBcount = -1, worstIBcount = -1;
  int worstOBstep = -1, worstIBstep = -1;
  Long_t worstOBorbit = 0, worstIBorbit = 0;

  double TimeStampFromStart[NSteps];
  double BarrelEfficiency[2][NSteps];

  double maprange = (double)(MAP.rbegin()->first - MAP.begin()->first) * 89.e-6; 
  QALOG<<"Orbit range "<<MAP.begin()->first<<" to "<<MAP.rbegin()->first<<" , in seconds: "<<maprange<<"\n";
  
  for (auto M : MAP){

    if (firstorbit == 0) firstorbit = currentorbit = previousorbit = M.first;
    previousorbit = currentorbit;
    currentorbit = M.first;
    if (countstep > 0)  {
      Long_t ogap = currentorbit-previousorbit;
      hOrb->SetBinContent(countstep, ogap);
      maxgap = TMath::Max(maxgap,ogap);
      if (ogap > NominalGap){
	ngap_overnominal++;
      }
      if (ogap > UnanchorableThreshold){
	unAnchorable += (ogap - UnanchorableThreshold);
      }
    }

    
    currentorbit = M.first;

    TimeStampFromStart[countstep] = (currentorbit - firstorbit) * LHCOrbitNS * 1.e-9;
    

    //if (M.second.size() > 50) cout<<"DEBUG Step "<<countstep<<"/"<<NSteps<<". Orbit = "<<M.first<<" dead lanes: "<<M.second.size();

    int OBdead = 0, IBdead = 0;
    for (int il = 0; il<2; il++) BarrelEfficiency[il][countstep] = 0.;
    
    for (uint lan : M.second) {

      if (lan < N_LANES_IB) hStatusTimeIB->Fill(currentorbit+1, lan); // "+1" to avoid edge effect in conversion from long to double
      else if (lan < N_LANES_IB + N_LANES_ML) hStatusTimeML->Fill(currentorbit+1, lan);
      else hStatusTimeOL->Fill(currentorbit+1,lan);
      
      if (lan < N_LANES_IB) IBdead++;
      else OBdead++;
      int ibarrel = (int)(LaneToLayer[lan] > 2);
      BarrelEfficiency[ibarrel][countstep] += 1./((ibarrel <1) ? N_LANES_IB : (N_LANES - N_LANES_IB));
    }

    if ( (maprange < SecForWorse) || (currentorbit - firstorbit)*89.e-6 > SecForWorse ){ // save worst cases after SecForWorse seconds if the map lasts at least SecForWorse seconds
      if (OBdead > worstOBcount){
	worstOBstep = countstep;
	worstOBorbit = currentorbit;
	worstOBcount = OBdead;
      }
      if (IBdead > worstIBcount){
	worstIBstep = countstep;
	worstIBorbit = currentorbit;
	worstIBcount = IBdead;
      }   
    }
      
    if (countstep > 0){
      hEffOB->SetBinContent(countstep,1.*OBdead/(N_LANES-N_LANES_IB));
      hEffIB->SetBinContent(countstep,1.*IBdead/N_LANES_IB);
    }

    for (uint ii : M.second) eff[ii] += (currentorbit - previousorbit);

    // step up
    countstep++;
  
  } // end of for M:MAP

  

  QALOG<<"Max orbit gap: "<<maxgap<<"\n";
  QALOG<<"Number of gaps over nominal ("<<NominalGap<<"): "<<ngap_overnominal<<"\n";

  QAcheck["Orbit gaps"] = (maxgap > 2.*UnanchorableThreshold || ngap_overnominal > 0.25*NSteps ) ? "BAD" : (maxgap > NominalGap && ngap_overnominal > 2 ) ? "MEDIUM" : "GOOD";

  double unAnchorableFrac = 1.*unAnchorable/ (currentorbit-firstorbit);

  QALOG<<"Un-anchorable number of orbits: "<<unAnchorable<<", corresponding to a fraction of the run of "<<unAnchorableFrac<<"\n";

  QAcheck["Un-anchorable fraction"] = (unAnchorableFrac < 0.02) ? "GOOD" : (unAnchorableFrac < 0.05) ? "MEDIUM" : "BAD";

  QALOG<<"Worst cases computed skipping first "<<SecForWorse<<" seconds.\n";
  QALOG<<"Worst OB case: orbit "<<worstOBorbit<<" step #"<<worstOBstep<<" dead lanes "<<worstOBcount<<"\n";
  QALOG<<"Worst IB case: orbit "<<worstIBorbit<<" step #"<<worstIBstep<<" dead lanes "<<worstIBcount<<"\n";

  if ((MAP.begin()->first != firstorbit) || (MAP.rbegin()->first != currentorbit)){
    QALOG<<"ERROR after checking the map the first and last orbit don't match\n";
  }

  for (uint chip : MAP.rbegin()->second) LastMAP->SetBinContent(chip+1,1); 

  double AvgEffIB=0, AvgEffOB=0, cib=0, cob=0; // average dead time over lanes. eff[i] is the dead time for lane i
  double AvgEffStave[N_STAVES]; double cst[N_STAVES]; for (int i=0; i<N_STAVES;i++) AvgEffStave[i]=cst[i]=0; // average dead time staves
  
  for (int i=0; i<nn; i++) {
    eff[i] = (1.*eff[i]/(currentorbit-firstorbit)); // normalizing eff[i] to time
    if (i<N_LANES_IB) {AvgEffIB += eff[i]; cib+=1;}
    else {AvgEffOB += eff[i]; cob+=1;}
    AvgEffStave[LaneToStave[i]] += eff[i];
    cst[LaneToStave[i]] += 1;
  }
  if (cib > 0) AvgEffIB /= cib; else AvgEffIB = 1.1111;
  if (cob > 0) AvgEffOB /= cob; else AvgEffOB = 1.1111;
  for (int i=0; i<N_STAVES; i++){
    if (cst[i]>0) AvgEffStave[i] /= cst[i]; else AvgEffStave[i] = 1.1111;
  }
  

  QALOG<<"Average IB dead time: "<<AvgEffIB<<"\n";
  QALOG<<"Average OB dead time: "<<AvgEffOB<<"\n";

  QAcheck["Avg dead time IB"] = (AvgEffIB < 0.03) ? "GOOD" : (AvgEffIB < 0.10) ? "MEDIUM" : "BAD";
  QAcheck["Avg dead time OB"] = (AvgEffOB < 0.05) ? "GOOD" : (AvgEffOB < 0.10) ? "MEDIUM" : "BAD";
    
  for (int i=0; i<nn; i++) HMAP->SetBinContent(i+1,eff[i]);
  for (int i=0; i<N_STAVES; i++) if (AvgEffStave[i]>0) hStaveDeadTime->SetBinContent(i+1, AvgEffStave[i]);
  if (worstOBorbit > 0) for (uint chip : MAP[worstOBorbit]) WorstOB->SetBinContent(chip+1,1);
  if (worstIBorbit > 0) for (uint chip : MAP[worstIBorbit]) WorstIB->SetBinContent(chip+1,1);
  WorstOB->SetTitle(Form("Worst OB, step %d = %d sec",worstOBstep,(int)((worstOBorbit-firstorbit)*89.e-6)));
  WorstOB->SetName(Form("Worst OB, step %d = %d sec",worstOBstep,(int)((worstOBorbit-firstorbit)*89.e-6)));
  WorstIB->SetTitle(Form("Worst IB, step %d = %d sec",worstIBstep,(int)((worstIBorbit-firstorbit)*89.e-6)));
  WorstIB->SetName(Form("Worst IB, step %d = %d sec",worstIBstep,(int)((worstIBorbit-firstorbit)*89.e-6)));

  QALOG<<"Comparing orbit range and run duration\n";

  GetTimeStamps(runnumber, firstorbit, currentorbit);

  QALOG<<"Assuming run start = 0. Map start = "<<(mapstart-runstart)/1000.<<", Run stop = "<<(runstop-runstart)/1000.<<", Map stop = "<<(mapstop-runstart)/1000.<<" (sec)\n";

  double RCTrunduration = (runstop - runstart)/1000.;
  double MAPduration = (currentorbit-firstorbit) * (LHCOrbitNS *1.e-9);

  QALOG<<"RCT run duration: "<<RCTrunduration<<". MAP duration: "<<MAPduration<<". Difference: "<<RCTrunduration - MAPduration<<" (sec)\n";

  QAcheck["Orbit range"] =
    (MAPduration > RCTrunduration) ? "MEDIUM":
    (MAPduration >= RCTrunduration - 5) ? "GOOD" :
    (MAPduration >= RCTrunduration - 30) ? "MEDIUM" :
    "BAD";

  hTimeSpan->SetBinContent(1,RCTrunduration);
  hTimeSpan->SetBinContent(2,MAPduration);

  TFile outroot(Form("%s/DeadMapQA.root",outdir.Data()),"RECREATE");
 
  c1->cd(1); // average evolving map
  HMAP->Draw("lcolz");
  gPad->SetLogz();
  HMAP->Write();
  
  c1->cd(2); // static map
  HSMAP->Draw("lcolz");
  HSMAP->SetMinimum(1);
  HSMAP->Write();

  c1->cd(3);
  hTimeSpan->GetXaxis()->SetBinLabel(1,"RCT");
  hTimeSpan->GetXaxis()->SetBinLabel(2,"MAP");
  hTimeSpan->GetYaxis()->SetTitle("sec");
  hTimeSpan->Draw("histo");
  hTimeSpan->Write();

  c1->cd(4);

  // second row
  
  c1->cd(5); // orbit offset
  hOrb->Draw("histo");
  hOrb->Write();

  c1->cd(6); // IB and OB efficiency
  hEffOB->GetXaxis()->SetTitle("step");
  hEffOB->Draw("histo");
  hEffIB->SetLineColor(2);
  hEffIB->Draw("histo same");
  gPad->SetLogy();
  hEffOB->Write();
  hEffIB->Write();

  c1->cd(7); // IB and OB efficiency vs time
  TGraph *grIB = new TGraph(NSteps,TimeStampFromStart,BarrelEfficiency[0]);
  TGraph *grOB = new TGraph(NSteps,TimeStampFromStart,BarrelEfficiency[1]);
  grIB->SetLineColor(kRed);
  grOB->SetLineColor(kBlue);
  grOB->SetMinimum(0);
  grOB->SetMaximum(1);
  grOB->SetTitle("Dead fraction vs time");
  grOB->GetXaxis()->SetTitle("time (sec)");
  grOB->Draw();
  grIB->Draw("same");
  gPad->SetLogy();
  grOB->SetName("Dead fraction vs time OB");
  grIB->SetName("Dead fraction vs time IB");
  grOB->Write();
  grIB->Write();

  c1->cd(8);

  // third row

  c1->cd(9); // worst OB
  WorstOB->Draw("lcolz");
  WorstOB->Write();
  
  c1->cd(10); // worst IB
  WorstIB->Draw("lcolz");
  WorstIB->Write();

  c1->cd(11); // last snapshot
  LastMAP->SetTitle("Last snapshot");
  LastMAP->SetName("Last snapshot");
  LastMAP->Draw("lcolz");
  LastMAP->Write();

  c1->cd(4); // text summary
  TLatex latex;
  latex.SetNDC();
  double yPos = .9;
  latex.SetTextSize(0.04);
  latex.SetTextAlign(13); // left alignment
  latex.DrawLatex(0.1,yPos,Form("File %s",FILENAME.Data())); 
  yPos -= 0.05;
  latex.DrawLatex(0.1,yPos,Form("Run %d",runnumber));
  yPos -= 0.05;
  for (auto cc : QAcheck){
    yPos -= 0.05;
    if (cc.second == "GOOD"){
      latex.SetTextColor(kGreen);
    }
    else if (cc.second == "MEDIUM"){
      latex.SetTextColor(kOrange);
    }
    else{
      latex.SetTextColor(kRed);
    }
    latex.DrawLatex(0.1,yPos,Form("%s : %s",cc.first.Data(),cc.second.Data()));
  }

  c2->cd(1);
  hStatusTimeIB->GetYaxis()->SetNdivisions(0.-(N_LANES_IB/9));
  hStatusTimeIB->Draw("col");    
  hStatusTimeIB->Write();
  gPad->SetGrid(0,1);

  c2->cd(2);
  hStatusTimeML->GetYaxis()->SetNdivisions(0.-(N_LANES_ML/16));
  hStatusTimeML->Draw("col");
  hStatusTimeML->Write();
  gPad->SetGrid(0,1);

  c2->cd(3);
  hStatusTimeOL->GetYaxis()->SetNdivisions(0.-(N_LANES-N_LANES_ML-N_LANES_IB/28));
  hStatusTimeOL->Draw("col");
  hStatusTimeOL->Write();
  //gPad->SetGrid(0,1);

  c3->cd();
  hStaveDeadTime->GetXaxis()->SetNdivisions(-64);
  hStaveDeadTime->GetXaxis()->SetLabelSize(0.02);
  hStaveDeadTime->SetMarkerStyle(20);
  hStaveDeadTime->Draw("P");
  TLatex latexl;
  //latexl.SetTextAlign(23); // Center alignment
  latexl.SetTextAngle(90); // Rotate text 90 degrees for vertical display
  latexl.SetTextSize(0.0125); 
  for (int i = 0; i < N_STAVES; i++) {
        double binContent = hStaveDeadTime->GetBinContent(i+1);
        double binCenter = hStaveDeadTime->GetBinCenter(i+1);
        if (binContent != 0) latexl.DrawLatex(binCenter, binContent*1.03, hStaveDeadTime->GetXaxis()->GetBinLabel(i+1));      
  }
  hStaveDeadTime->Write();
  gPad->SetLogy();
  gPad->SetGrid(1,1);
  c3->Update();
  

  c1->Write();
  c2->Write();
  c3->Write();
  c1->SaveAs(Form("%s/DeadMapQA1.png",outdir.Data()));
  c2->SaveAs(Form("%s/DeadMapQA2.png",outdir.Data()));
  c3->SaveAs(Form("%s/DeadMapQA3.png",outdir.Data()));

  outroot.Close();
  
  QALOG<<"Orbits: "<<firstorbit<<" to "<<currentorbit<<" corrsponding to "<<(currentorbit - firstorbit)*89.e-6 / 60.<<" minutes\n";

  
  
  if (ExitWhenFinish){
    PrintAndExit();
  }
 
}




void interpolatestave(double r1, double r2, double phi1, double phi2, int z, int nz, int a, int na, TVector2 *TVec){

  TVector2 extr[4];

  double phi1pp = phi1 + 0.05*(phi2-phi1);
  double phi2pp = phi2 - 0.05*(phi2-phi1);
  extr[0].SetMagPhi(r1,phi1pp);
  extr[1].SetMagPhi(r2,phi1pp);
  extr[2].SetMagPhi(r1,phi2pp);
  extr[3].SetMagPhi(r2,phi2pp);

  TVector2 seg[4] = {extr[0]+(1.*z/nz)*(extr[1]-extr[0]), extr[2]+(1.*z/nz)*(extr[3]-extr[2]),
		     extr[0]+(1.*(z+1)/nz)*(extr[1]-extr[0]), extr[2]+(1.*(z+1)/nz)*(extr[3]-extr[2])};

  /* seg:

  | | | |
  3 | | 2
  | | | |
  1 | | 0
  | | | |

  */

  TVec[0] = seg[0]+ (1.*a/na)*(seg[1]-seg[0]);
  TVec[1] = seg[0]+ (1.*(a+1)/na)*(seg[1]-seg[0]);

  TVec[2] = seg[2]+ (1.*(a+1)/na)*(seg[3]-seg[2]);
  TVec[3] = seg[2]+ (1.*(a)/na)*(seg[3]-seg[2]);
  
}

void RemoveAxis(TH2Poly *HP){
  HP->GetXaxis()->SetLabelSize(0);
  HP->GetXaxis()->SetTickLength(0);
  HP->GetYaxis()->SetLabelSize(0);
  HP->GetYaxis()->SetTickLength(0);
}

void getlanecoordinates(int laneid, double *px, double *py){

  int laneob = laneid - N_LANES_IB;
  int layer, staveinlayer, laneinlayer;

  if (laneid < 12*9) {laneinlayer = laneid; layer=0;}
  else if (laneid < (12+16)*9) {laneinlayer = laneid - (12*9); layer=1;}
  else if (laneid < N_LANES_IB) { laneinlayer = laneid - (12+16)*9; layer=2;}
  else if (laneob < 24*4*4) {layer = 3; laneinlayer = laneob;}
  else if (laneob < (24+30)*4*4) { layer = 4; laneinlayer = laneob - 24*4*4;}
  else if (laneob < (24+30)*4*4 + 42*7*4) { layer = 5; laneinlayer = laneob - (24+30)*4*4;}
  else { layer = 6; laneinlayer = laneob - (24+30)*4*4 - 42*7*4;}

  LaneToLayer[laneid] = layer;

  int nz = NZElementsInHalfStave[layer];
  int nseg = NSegmentsStave[layer];
  int nsegh = (nseg == 1) ? 1 : nseg/2;

  
  staveinlayer = laneinlayer / (nz*nseg);
  int laneinstave = laneinlayer % (nz*nseg);
  int stave = 0;
  for (int l=0; l<7; l++) stave += (l<layer)*NStaves[l]+(l==layer)*staveinlayer;

  LaneToStave[laneid] = stave;
  LaneToStaveInLayer[laneid] = staveinlayer;

  int halfstave = (layer < 3) ? 0 : (int)( laneinstave >= nz*nsegh); // 0 or 1 
  int laneinhalfstave = laneinstave - halfstave*nz*nsegh;

  int phi_matrix = (int)(laneinhalfstave % nsegh) + nsegh*halfstave;
  int z_matrix = (int)(laneinhalfstave/nsegh);

  double r1 = 2. + 0.4*layer + 9*layer + 2*(layer > 2);
  double r2 = r1+9;
  double phi1 = Pi()*2 *staveinlayer / NStaves[layer];
  double phi2 = phi1 + Pi()*2 / NStaves[layer];
  
  TVector2 *tvec = new TVector2[4];

  interpolatestave(r1, r2, phi1, phi2, z_matrix, NZElementsInHalfStave[layer], phi_matrix, nseg, tvec);

  for (int i=0; i<4; i++) {
    px[i] = tvec[i].X();
    py[i] = tvec[i].Y();
  }
}


uint16_t isFirstOfLane(uint16_t chipid){ // 9999 if it is not
  
  if (chipid < N_LANES_IB){
    return chipid;
  }
  else if ((chipid - N_LANES_IB) % 7 == 0){
    return N_LANES_IB + (uint16_t)( (chipid - N_LANES_IB) / 7);
  }
  return 9999;
}


uint16_t isLastOfLane(uint16_t chipid){ // 9999 if it is not

  if (chipid < N_LANES_IB){
    return chipid;
  }
  else if (chipid < N_LANES_IB + 7 - 1){
    return 9999;
  }
  else return isFirstOfLane(chipid - 6);
}

uint16_t StaveToLayer(uint16_t stv){
  int count = 0;
  for (int i=0; i<7; i++){
    count+=NStaves[i];
    if (stv < count) return (uint16_t)i;
  }
}

uint16_t FirstLaneOfStave(uint16_t stv){

  uint16_t lay = StaveToLayer(stv);
  if (lay < 3) return stv*9;
  else if (lay < 5) return N_LANES_IB + (stv - NStaves[0] - NStaves[1] - NStaves[2])*16;
  else return N_LANES_IB + N_LANES_ML + (stv - NStaves[0] - NStaves[1] - NStaves[2] - NStaves[3] - NStaves[4]) * 28;
}


uint16_t LastLaneOfStave(uint16_t stv){

  uint16_t lay = StaveToLayer(stv);
  if (lay < 3) return FirstLaneOfStave(stv)+8;
  if (lay < 5) return FirstLaneOfStave(stv)+15;
  else return FirstLaneOfStave(stv)+27;
}


  
std::vector<uint16_t> expandvector(std::vector<uint16_t> words, std::string version, TString opt = "lane"){

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

      if (opt == "chip"){
	QALOG<<"Static map: decoded interval of "<<lastel - firstel + 1<<" chips: "<<firstel<<":"<<lastel<<"\n";
	for (uint16_t ic = firstel; ic <= lastel; ic++){
	  elementlist.push_back(ic);
	}
      }

      else if (opt == "lane"){

	uint16_t firstlane = isFirstOfLane(firstel);
	uint16_t lastlane = isLastOfLane(lastel);

	if (firstlane == 9999 || lastlane == 9999){
	  QALOG<<"FATAL decoding of chip intervals, returning emtpy vector\n";
	  QAcheck["Chip interval"] = "FATAL";
	  return elementlist;
	}
	else{
	  for (uint16_t il = firstlane ; il <= lastlane; il++){
	    elementlist.push_back(il);
	  }
	}
      } // end opt == lane
    } // end loop word
  } // end map version


  else{
    QALOG<<"FATAL: map version not recognized, returning empty vector.\n";
    QAcheck["MAP version"] = "FATAL";
  }

  return elementlist;
}


uint16_t ChipToLane(uint16_t chipid){

  if (chipid < N_LANES_IB) return chipid;
  else return N_LANES_IB + (uint16_t)(( chipid - N_LANES_IB) / 7);
}

void fillmap(TString fname){

  MAP.clear();
  MAPKeys.clear();
  SMAP.clear();

  TFile *f = new TFile(fname);
  o2::itsmft::TimeDeadMap* obj = nullptr;
  
  f->GetObject("ccdb_object",obj);
  f->Close();

  if (!obj){
    QAcheck["ROOT file corrupted"] = "FATAL";
    PrintAndExit("FATAL object not found");
  }
  

  QALOG<<"Reading map info\n";
  std::string mapver = obj->getMapVersion();
  QALOG<<"Map version = "<<mapver<<"\n";
  QALOG<<"Is default object = "<<obj->isDefault()<<"\n";
  QALOG<<"Evolving map size = "<<obj->getEvolvingMapSize()<<"\n";

  std::vector<uint16_t> StaticMap;
  obj->getStaticMap(StaticMap);

  QALOG<<"Static map imported\n";
  QALOG<<"Static map size = "<<StaticMap.size()<<"\n";

  QAcheck["Map size"] = (obj->getEvolvingMapSize() > 0 && StaticMap.size() > 0) ? "GOOD" : "BAD";
  if (obj->isDefault()){
    if (obj->getEvolvingMapSize() == 0 && StaticMap.size() == 0){
      QAcheck["Map size"] = "GOOD";
    }
    QAcheck["Default object"] = "FATAL";
    PrintAndExit("Exiting because default object");
  }
  if (obj->getEvolvingMapSize() == 0){
    QAcheck["Map size"] = "FATAL";
    PrintAndExit("Exiting because evolving map is empty");
  }

  SMAP = expandvector(StaticMap,mapver,"chip");

  MAPKeys = obj->getEvolvingMapKeys();

  QALOG<<"Orbit keys imported. Importing maps...\n";
  for (auto OO : MAPKeys){
    std::vector<uint16_t> MapAtOrbit;
    obj->getMapAtOrbit(OO, MapAtOrbit);
    MAP[OO] = expandvector(MapAtOrbit,mapver,"lane");
  }

  QALOG<<"... done importing maps for every orbit.\n";
}

void GetTimeStamps(int runnumber, uint32_t orbit1, uint32_t orbit2){

  if (runnumber <= 0) {
    QALOG<<"WARNING empty runnumber provided. Checks on run duration will not be effective\n";
    return;
  }
  
  QALOG<<"Reading RCT and CTP info for run "<<runnumber<<" from "<<ccdbHost<<"\n";
  static int64_t orbitResetMUS = 0;
  auto& cm = o2::ccdb::BasicCCDBManager::instance();
  cm.setURL(ccdbHost);
  auto lims = cm.getRunDuration(runnumber);
  if (lims.first == 0 && lims.second == 0) {
    QALOG<<"ERROR, failed to fetch run info from RCT";
    return;
  }

  runstart = (long)lims.first;
  runstop = (long)lims.second;
  
  auto* orbitReset = cm.getForTimeStamp<std::vector<Long64_t>>("CTP/Calib/OrbitReset", lims.first);
  orbitResetMUS = (*orbitReset)[0];

  //cout<<"DEBUG "<<orbit1<<" "<<LHCOrbitNS<<" "<<orbitResetMUS<<" ";

  mapstart = std::ceil((orbit1 * LHCOrbitNS / 1000 + orbitResetMUS) / 1000);
  mapstop = std::ceil((orbit2 * LHCOrbitNS / 1000 + orbitResetMUS) / 1000);

  QALOG<<"run start/stop: "<<runstart<<"/"<<runstop<<" map start/stop: "<<mapstart<<"/"<<mapstop<<" (ms)\n";
}
    

