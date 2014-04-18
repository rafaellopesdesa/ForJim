//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Apr 15 10:18:57 2014 by ROOT version 5.22/00a
// from TTree Global/Fake pmcs output
// found on file: /home/rclsa/2673/SAM_SAMPLES4/resbos_wen_cteq66_v2_0.root
//////////////////////////////////////////////////////////

#ifndef Global_h
#define Global_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TMath.h>
#include <TVector2.h>
#include <TLorentzVector.h>
#include <iostream>

class Global {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   TVector2 elec, nu;
   TLorentzVector WZ;
   Double_t mT;
   Double_t ElecpT;
   Double_t MET;
   Double_t bw_default;

   // Declaration of leaf types
   Int_t           pmcs_ana_nevtp;
   Int_t           pmcs_ana_npart;
   Int_t           pmcs_ana_nvtx;
   Float_t         pmcs_ana_evnum[1000];   //[nevtp]
   Float_t         pmcs_ana_evrun[1000];   //[nevtp]
   Float_t         pmcs_ana_evwt[1000];   //[nevtp]
   Float_t         pmcs_ana_evxs[1000];   //[nevtp]
   Float_t         pmcs_ana_pE[1000];   //[npart]
   Float_t         pmcs_ana_pcid[1000];   //[npart]
   Float_t         pmcs_ana_pcnum[1000];   //[npart]
   Float_t         pmcs_ana_pdvtx[1000];   //[npart]
   Float_t         pmcs_ana_peta[1000];   //[npart]
   Float_t         pmcs_ana_pidx[1000];   //[npart]
   Float_t         pmcs_ana_pistable[1000];   //[npart]
   Float_t         pmcs_ana_pphi[1000];   //[npart]
   Float_t         pmcs_ana_ppid[1000];   //[npart]
   Float_t         pmcs_ana_ppt[1000];   //[npart]
   Float_t         pmcs_ana_ppvtx[1000];   //[npart]
   Float_t         pmcs_ana_ppx[1000];   //[npart]
   Float_t         pmcs_ana_ppy[1000];   //[npart]
   Float_t         pmcs_ana_ppz[1000];   //[npart]
   Float_t         pmcs_ana_vcid[1000];   //[nvtx]
   Float_t         pmcs_ana_vcnum[1000];   //[nvtx]
   Float_t         pmcs_ana_vct[1000];   //[nvtx]
   Float_t         pmcs_ana_vidx[1000];   //[nvtx]
   Float_t         pmcs_ana_visdisp[1000];   //[nvtx]
   Float_t         pmcs_ana_vpprt[1000];   //[nvtx]
   Float_t         pmcs_ana_vx[1000];   //[nvtx]
   Float_t         pmcs_ana_vy[1000];   //[nvtx]
   Float_t         pmcs_ana_vz[1000];   //[nvtx]
   Int_t           pmcs_em_nelg;
   Int_t           pmcs_em_nels;
   Int_t           pmcs_em_nphg;
   Int_t           pmcs_em_nphs;
   Float_t         pmcs_em_elcalphis[120];   //[nels]
   Float_t         pmcs_em_eleg[120];   //[nelg]
   Float_t         pmcs_em_elelmergedEg[120];   //[nelg]
   Float_t         pmcs_em_eles[120];   //[nelg]
   Float_t         pmcs_em_eletads[120];   //[nels]
   Float_t         pmcs_em_eletag[120];   //[nelg]
   Float_t         pmcs_em_eletas[120];   //[nels]
   Float_t         pmcs_em_elfid[120];   //[nelg]
   Int_t           pmcs_em_elhastrack[120];   //[nels]
   Float_t         pmcs_em_eliso[120];   //[nels]
   Int_t           pmcs_em_elmergedg[120];   //[nelg]
   Int_t           pmcs_em_elpasshmtx[120];   //[nels]
   Int_t           pmcs_em_elpassid1011[120];   //[nels]
   Float_t         pmcs_em_elphig[120];   //[nelg]
   Float_t         pmcs_em_elphis[120];   //[nels]
   Float_t         pmcs_em_elphmergedEg[120];   //[nelg]
   Int_t           pmcs_em_elpntg[120];   //[nels]
   Int_t           pmcs_em_elpnts[120];   //[nelg]
   Float_t         pmcs_em_elptg[120];   //[nelg]
   Float_t         pmcs_em_elpts[120];   //[nels]
   Int_t           pmcs_em_elpttr[120];   //[nels]
   Float_t         pmcs_em_pheg[120];   //[nphg]
   Float_t         pmcs_em_phes[120];   //[nphs]
   Float_t         pmcs_em_phetads[120];   //[nphs]
   Float_t         pmcs_em_phetag[120];   //[nphg]
   Float_t         pmcs_em_phetas[120];   //[nphs]
   Float_t         pmcs_em_phfid[120];   //[nphg]
   Int_t           pmcs_em_phhastrack[120];   //[nels]
   Float_t         pmcs_em_phiso[120];   //[nphs]
   Int_t           pmcs_em_phpasshmtx[120];   //[nels]
   Float_t         pmcs_em_phphig[120];   //[nphg]
   Float_t         pmcs_em_phphis[120];   //[nphs]
   Int_t           pmcs_em_phpntg[120];   //[nphs]
   Int_t           pmcs_em_phpnts[120];   //[nphg]
   Float_t         pmcs_em_phptg[120];   //[nphg]
   Float_t         pmcs_em_phpts[120];   //[nphs]
   Int_t           pmcs_met_nmetg;
   Int_t           pmcs_met_nmets;
   Float_t         pmcs_met_metg[200];   //[nmetg]
   Float_t         pmcs_met_metphig[200];   //[nmetg]
   Float_t         pmcs_met_metphis[200];   //[nmets]
   Float_t         pmcs_met_mets[200];   //[nmets]
   Float_t         pmcs_met_metxg[200];   //[nmetg]
   Float_t         pmcs_met_metxs[200];   //[nmets]
   Float_t         pmcs_met_metyg[200];   //[nmetg]
   Float_t         pmcs_met_metys[200];   //[nmets]
   Float_t         pmcs_met_scalarg[200];   //[nmetg]
   Float_t         pmcs_met_scalars[200];   //[nmets]
   Int_t           pmcs_vtx_nvtxg;
   Int_t           pmcs_vtx_nvtxs;
   Float_t         pmcs_vtx_vtnds[150];   //[nvtxs]
   Float_t         pmcs_vtx_vtxxs[150];   //[nvtxs]
   Float_t         pmcs_vtx_vtxys[150];   //[nvtxs]
   Float_t         pmcs_vtx_vtxzg[150];   //[nvtxg]
   Float_t         pmcs_vtx_vtxzs[150];   //[nvtxs]

   // List of branches
   TBranch        *b_pmcs_ana;   //!
   TBranch        *b_pmcs_em;   //!
   TBranch        *b_pmcs_met;   //!
   TBranch        *b_pmcs_vtx;   //!

   Global(const char* fileName);
   virtual ~Global();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop(Bool_t doTemplates, Int_t maxEvents);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
   virtual Double_t Reweight(Double_t mass);
};

#endif

#ifdef Global_cxx
Global::Global(const char* fileName)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
  TFile* f = TFile::Open(fileName);
  TTree* tree = (TTree*) f->Get("Global");
  Init(tree);
}

Global::~Global()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t Global::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t Global::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (!fChain->InheritsFrom(TChain::Class()))  return centry;
   TChain *chain = (TChain*)fChain;
   if (chain->GetTreeNumber() != fCurrent) {
      fCurrent = chain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void Global::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("pmcs_ana", &pmcs_ana_nevtp, &b_pmcs_ana);
   fChain->SetBranchAddress("pmcs_em", &pmcs_em_nelg, &b_pmcs_em);
   fChain->SetBranchAddress("pmcs_met", &pmcs_met_nmetg, &b_pmcs_met);
   fChain->SetBranchAddress("pmcs_vtx", &pmcs_vtx_nvtxg, &b_pmcs_vtx);
   Notify();
}

Bool_t Global::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void Global::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t Global::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.

  for (Int_t ipart=0; ipart<pmcs_ana_npart; ipart++) {
    if (TMath::Abs(pmcs_ana_ppid[ipart]) == 11) { // electron
      if (TMath::Abs(pmcs_ana_peta[ipart]) > 1.05) return -1;
      if (TMath::Abs(pmcs_ana_ppt[ipart]) < 25) return -1;      
      elec.Set(pmcs_ana_ppx[ipart], pmcs_ana_ppy[ipart]);
    } else if (TMath::Abs(pmcs_ana_ppid[ipart]) == 12) { // electron
      if (TMath::Abs(pmcs_ana_ppt[ipart]) < 25) return -1;
      nu.Set(pmcs_ana_ppx[ipart], pmcs_ana_ppy[ipart]);
    } else if (TMath::Abs(pmcs_ana_ppid[ipart]) == 24) { // electron
      if (TMath::Abs(pmcs_ana_ppt[ipart]) > 15) return -1;
      WZ.SetXYZT(pmcs_ana_ppx[ipart], pmcs_ana_ppy[ipart], pmcs_ana_ppz[ipart], pmcs_ana_pE[ipart]);
    }
  }
  mT = TMath::Sqrt(2*(elec.Mod()*nu.Mod()-elec*nu));
  ElecpT = elec.Mod();
  MET = nu.Mod();
  if (mT < 50 || mT>200) return -1;
  bw_default = TMath::Power(WZ.M()*2.048/80.419,2.)/(TMath::Power(WZ.M()*WZ.M()-80.419*80.419,2.) + TMath::Power(WZ.M()*WZ.M()*2.048/80.419,2.));
  return 1;
}

Double_t Global::Reweight(Double_t mass)
{
  Double_t bw_reweight = TMath::Power(WZ.M()*2.048/mass,2.)/(TMath::Power(WZ.M()*WZ.M()-mass*mass,2.)+TMath::Power(WZ.M()*WZ.M()*2.048/mass,2.));
  return bw_reweight/bw_default;
}
  

#endif // #ifdef Global_cxx
