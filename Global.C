#define Global_cxx
#include "Global.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TString.h>
#include <TFile.h>
#include <iostream>

using namespace std;

int main(int argc, char** argv)
{
  Global t(argv[1]);
  if (argc > 3)
    t.Loop(TString(argv[2]).Atoi(), TString(argv[3]).Atoi());
  else
    t.Loop(TString(argv[2]).Atoi(), -1);
}

void Global::Loop(Bool_t doTemplates, Int_t maxEvents)
{

  TH1D* map = new TH1D("map", "", 100, 89., 92.);
  TH1D* mt_pseudo = new TH1D("mt_pseudo", "", 300, 50., 200.);
  TH1D* mt_templates[100];
  TH1D* pt_pseudo = new TH1D("pt_pseudo", "", 200, 0., 100.);
  TH1D* pt_templates[100];
  TH1D* met_pseudo = new TH1D("met_pseudo", "", 200, 0., 100.);
  TH1D* met_templates[100];
  for (Int_t i=0; i<100; i++) {
    mt_templates[i] = new TH1D(TString::Format("mt_templates_%d",i), "", 300, 50., 200.);
    pt_templates[i] = new TH1D(TString::Format("pt_templates_%d",i), "", 200, 0., 100.);
    met_templates[i] = new TH1D(TString::Format("met_templates_%d",i), "", 200, 0., 100.);
  }

  if (fChain == 0) return;
  
  Long64_t nentries = fChain->GetEntriesFast();
  Int_t ct = 0;
  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    if (Cut(ientry) < 0) continue;     
    mt_pseudo->Fill(mT);
    pt_pseudo->Fill(ElecpT);
    met_pseudo->Fill(MET);
    if (doTemplates) {
      for (Int_t i=0; i<100; i++) {
	mt_templates[i]->Fill(mT, Reweight(80.419+(((Double_t) i)-50)*0.02));
	pt_templates[i]->Fill(ElecpT, Reweight(80.419+(((Double_t) i)-50)*0.02));
	met_templates[i]->Fill(MET, Reweight(80.419+(((Double_t) i)-50)*0.02));
      }
    }
    ct++;
    if (ct > maxEvents && maxEvents > 0)
      break;
  }
  
  
  // Now convert to "PMCS"
  TFile* output = new TFile("output.root", "recreate");
  output->mkdir("default");
  output->cd("default");
  
  map->GetXaxis()->Set(99, 80.419-50*0.002, 80.419+49*0.002);
  map->Write("histd1map_WMassTemplates");
  
  mt_pseudo->Write("hWcandMt_CC");
  pt_pseudo->Write("hWcandElecPt_CC");
  met_pseudo->Write("hWcandMet_CC");
  if (doTemplates) {
    for (Int_t i=0; i<100; i++) {
      mt_templates[i]->Write(TString::Format("hWcandMt_CC_%d", i));
      pt_templates[i]->Write(TString::Format("hWcandElecPt_CC_%d", i));
      met_templates[i]->Write(TString::Format("hWcandMet_CC_%d", i));
    }
  }
  output->Close();

}


