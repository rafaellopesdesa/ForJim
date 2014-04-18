#include <TFile.h>
#include <TString.h>
#include <TGraphErrors.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TMath.h>
#include <TMatrixD.h>
#include <TVectorD.h>
#include <TH2.h>
#include <TF2.h> 
#include <TCanvas.h>
#include <TPad.h>
#include <TGaxis.h>
#include <TLegend.h>

#include <fstream>
#include <iostream>

#undef doGaussian
#define rebin 1

using namespace std;

Double_t chiSquare(Double_t* x, Double_t* p) {

  return p[0]*TMath::Power(x[0],p[1]/2.-1.)*TMath::Exp(-x[0]/2.);

}

Double_t gaussian2d(Double_t* x, Double_t* p) {

  return p[0]*TMath::Exp(-(1./(2.*(1-p[5]))) * ((x[0]-p[1])*(x[0]-p[1])/(p[3]*p[3])+(x[1]-p[2])*(x[1]-p[2])/(p[4]*p[4])-2*p[5]*(x[0]-p[1])*(x[1]-p[2])/(p[3]*p[4])));

}

int main(int argc, char** argv)
{

  gROOT->SetStyle("Plain");
  gStyle->SetOptStat("MR");
  gStyle->SetOptFit(0111);
  gStyle->SetTitleBorderSize(0);
  gStyle->SetPalette(1, 0); 

  // The ResBos mass
  const Double_t trueMass = 80.419;

  // Define what we want
  TH1D* mt_dists[300/rebin];
  TH1D* pt_dists[200/rebin];
  TH1D* met_dists[200/rebin];

  TH1D* mt_poisson[300/rebin];
  TH1D* pt_poisson[200/rebin];
  TH1D* met_poisson[200/rebin];

  TH1D* mtUnbiasedPull = new TH1D("mtUnbiasedPull", "", 50, -5., 5.);
  TH1D* ptUnbiasedPull = new TH1D("ptUnbiasedPull", "", 50, -5., 5.);
  TH1D* metUnbiasedPull = new TH1D("metUnbiasedPull", "", 50, -5., 5.);

  TH2D* mtptUnbiasedPull = new TH2D("mtptUnbiasedPull", "", 20, -5., 5., 20, -5., 5.);
  TH2D* ptmetUnbiasedPull = new TH2D("ptmetUnbiasedPull", "", 20, -5., 5., 20, -5., 5.);
  TH2D* metmtUnbiasedPull = new TH2D("metmtUnbiasedPull", "", 20, -5., 5., 20, -5., 5.);

  TH1D* mtPull = new TH1D("mtPull", "", 50, -5., 5.);
  TH1D* ptPull = new TH1D("ptPull", "", 50, -5., 5.);
  TH1D* metPull = new TH1D("metPull", "", 50, -5., 5.);

  TH2D* mtptPull = new TH2D("mtptPull", "", 20, -5., 5., 20, -5., 5.);
  TH2D* ptmetPull = new TH2D("ptmetPull", "", 20, -5., 5., 20, -5., 5.);
  TH2D* metmtPull = new TH2D("metmtPull", "", 20, -5., 5., 20, -5., 5.);
  
  TH1D* mtJKMean = new TH1D("mtJKMean", "", 50, 80.380, 80.440);
  TH1D* ptJKMean = new TH1D("ptJKMean", "", 50, 80.380, 80.440);
  TH1D* metJKMean = new TH1D("metJKMean", "", 50, 80.380, 80.440);

  TH1D* mtJKSigma = new TH1D("mtJKSigma", "", 50, 0.005, 0.060);
  TH1D* ptJKSigma = new TH1D("ptJKSigma", "", 50, 0.005, 0.060);
  TH1D* metJKSigma = new TH1D("metJKSigma", "", 50, 0.005, 0.060);

  TH1D* mtJKUnbiasedSigma = new TH1D("mtJKUnbiasedSigma", "", 50, 0.005, 0.060);
  TH1D* ptJKUnbiasedSigma = new TH1D("ptJKUnbiasedSigma", "", 50, 0.005, 0.060);
  TH1D* metJKUnbiasedSigma = new TH1D("metJKUnbiasedSigma", "", 50, 0.005, 0.060);

  TH1D* mtptJKCorr = new TH1D("mtptJKCorr", "", 20, 0.0, 1.0);
  TH1D* metmtJKCorr = new TH1D("metmtJKCorr", "", 20, 0.0, 1.0);
  TH1D* ptmetJKCorr = new TH1D("ptmetJKCorr", "", 20, 0.0, 1.0);

  // Storage for the data
  TGraph* mtptData = new TGraph();
  TGraph* metmtData = new TGraph();
  TGraph* ptmetData = new TGraph();

  TGraphErrors* mtptScatter = new TGraphErrors();
  TGraphErrors* metmtScatter = new TGraphErrors();
  TGraphErrors* ptmetScatter = new TGraphErrors();

  // Reads the data
  cout << "Getting the data" << endl;
  TH1D* mt_sum = new TH1D("mt_sum", "", 300/rebin, 50., 200.);
  TH1D* mt_sqrtmean = new TH1D("mt_sqrtmean", "", 300/rebin, 50., 200.);
  TH1D* mt_mean = new TH1D("mt_mean", "", 300/rebin, 50., 200.);
  TH1D* mt_rms = new TH1D("mt_rms", "", 300/rebin, 50., 200.);
  TH1D* mt_rms_ratio = new TH1D("mt_rms_ratio", "", 300/rebin, 50., 200.);
  TH1D* mt_rmsfrac = new TH1D("mt_rmsfrac", "", 300/rebin, 50., 200.);
  TH2D* mt_sum2 = new TH2D("mt_sum2", "", 300/rebin, 50., 200., 300/rebin, 50., 200.);
  TH2D* mt_corr = new TH2D("mt_corr", "", 300/rebin, 50., 200., 300/rebin, 50., 200.);
  TH1D* elecpt_sum = new TH1D("elecpt_sum", "", 200/rebin, 0., 100.);
  TH1D* elecpt_mean = new TH1D("elecpt_mean", "", 200/rebin, 0., 100.);
  TH1D* elecpt_sqrtmean = new TH1D("elecpt_sqrtmean", "", 200/rebin, 0., 100.);
  TH1D* elecpt_rms = new TH1D("elecpt_rms", "", 200/rebin, 0., 100.);
  TH1D* elecpt_rms_ratio = new TH1D("elecpt_rms_ratio", "", 200/rebin, 0., 100.);
  TH1D* elecpt_rmsfrac = new TH1D("elecpt_rmsfrac", "", 200/rebin, 0., 100.);
  TH2D* elecpt_sum2 = new TH2D("elecpt_sum2", "", 200/rebin, 0., 100., 200/rebin, 0., 100.);
  TH2D* elecpt_corr = new TH2D("elecpt_corr", "", 200/rebin, 0., 100., 200/rebin, 0., 100.);
  TH1D* met_sum = new TH1D("met_sum", "", 200/rebin, 0., 100.);
  TH1D* met_mean = new TH1D("met_mean", "", 200/rebin, 0., 100.);
  TH1D* met_sqrtmean = new TH1D("met_sqrtmean", "", 200/rebin, 0., 100.);
  TH1D* met_rms = new TH1D("met_rms", "", 200/rebin, 0., 100.);
  TH1D* met_rms_ratio = new TH1D("met_rms_ratio", "", 200/rebin, 0., 100.);
  TH1D* met_rmsfrac = new TH1D("met_rmsfrac", "", 200/rebin, 0., 100.);
  TH2D* met_sum2 = new TH2D("met_sum2", "", 200/rebin, 0., 100., 200/rebin, 0., 100.);
  TH2D* met_corr = new TH2D("met_corr", "", 200/rebin, 0., 100., 200/rebin, 0., 100.);
  Double_t corr_n = 0;

  ifstream _fileList(argv[2]);
  TString _reader;
  while (!_fileList.eof()) {
    _reader.ReadLine(_fileList);
    if (_reader.Length() < 5) continue;
    TFile* tempFile = TFile::Open(_reader);
    if (!tempFile->IsOpen()) continue;
    TH1D* tempmt = (TH1D*) tempFile->Get("default/hWcandMt_CC");
    TH1D* tempelecpt = (TH1D*) tempFile->Get("default/hWcandElecPt_CC");
    TH1D* tempmet = (TH1D*) tempFile->Get("default/hWcandMet_CC");
    tempmt->Rebin(rebin);
    tempelecpt->Rebin(rebin);
    tempmet->Rebin(rebin);
    for (int ibin=1; ibin<=tempmt->GetNbinsX(); ibin++) {
      mt_sum->Fill(tempmt->GetXaxis()->GetBinCenter(ibin), tempmt->GetBinContent(ibin));
      for (int jbin=1; jbin<=tempmt->GetNbinsX(); jbin++) {
	mt_sum2->Fill(tempmt->GetXaxis()->GetBinCenter(ibin), tempmt->GetXaxis()->GetBinCenter(jbin), tempmt->GetBinContent(ibin)*tempmt->GetBinContent(jbin));
      }
    }
    for (int ibin=1; ibin<=tempelecpt->GetNbinsX(); ibin++) {
      elecpt_sum->Fill(tempelecpt->GetXaxis()->GetBinCenter(ibin), tempelecpt->GetBinContent(ibin));
      for (int jbin=1; jbin<=tempelecpt->GetNbinsX(); jbin++) {
	elecpt_sum2->Fill(tempelecpt->GetXaxis()->GetBinCenter(ibin), tempelecpt->GetXaxis()->GetBinCenter(jbin), tempelecpt->GetBinContent(ibin)*tempelecpt->GetBinContent(jbin));
      }
    }
    for (int ibin=1; ibin<=tempmet->GetNbinsX(); ibin++) {
      met_sum->Fill(tempmet->GetXaxis()->GetBinCenter(ibin), tempmet->GetBinContent(ibin));
      for (int jbin=1; jbin<=tempmet->GetNbinsX(); jbin++) {
	met_sum2->Fill(tempmet->GetXaxis()->GetBinCenter(ibin), tempmet->GetXaxis()->GetBinCenter(jbin), tempmet->GetBinContent(ibin)*tempmet->GetBinContent(jbin));
      }
    }
    corr_n++;
    tempFile->Close();
  }
  _fileList.close();
  for (int ibin=1; ibin<=mt_rms->GetNbinsX(); ibin++) {
    mt_mean->SetBinContent(ibin, mt_sum->GetBinContent(ibin)/corr_n);
    mt_sqrtmean->SetBinContent(ibin, TMath::Sqrt(mt_sum->GetBinContent(ibin)/corr_n));
    mt_rms->SetBinContent(ibin, TMath::Sqrt(mt_sum2->GetBinContent(ibin, ibin)/corr_n-mt_sum->GetBinContent(ibin)*mt_sum->GetBinContent(ibin)/(corr_n*corr_n)));
    mt_rmsfrac->SetBinContent(ibin, mt_rms->GetBinContent(ibin)/mt_mean->GetBinContent(ibin));
  }  
  mt_rms_ratio->Divide(mt_rms, mt_sqrtmean, 1., 1.);
  for (int ibin=1; ibin<=elecpt_rms->GetNbinsX(); ibin++) {
    elecpt_mean->SetBinContent(ibin, elecpt_sum->GetBinContent(ibin)/corr_n);
    elecpt_sqrtmean->SetBinContent(ibin, TMath::Sqrt(elecpt_sum->GetBinContent(ibin)/corr_n));
    elecpt_rms->SetBinContent(ibin, TMath::Sqrt(elecpt_sum2->GetBinContent(ibin, ibin)/corr_n-elecpt_sum->GetBinContent(ibin)*elecpt_sum->GetBinContent(ibin)/(corr_n*corr_n)));
    elecpt_rmsfrac->SetBinContent(ibin, elecpt_rms->GetBinContent(ibin)/elecpt_mean->GetBinContent(ibin));
  }
  elecpt_rms_ratio->Divide(elecpt_rms, elecpt_sqrtmean, 1., 1.);
  for (int ibin=1; ibin<=met_rms->GetNbinsX(); ibin++) {
    met_mean->SetBinContent(ibin, met_sum->GetBinContent(ibin)/corr_n);
    met_sqrtmean->SetBinContent(ibin, TMath::Sqrt(met_sum->GetBinContent(ibin)/corr_n));
    met_rms->SetBinContent(ibin, TMath::Sqrt(met_sum2->GetBinContent(ibin, ibin)/corr_n-met_sum->GetBinContent(ibin)*met_sum->GetBinContent(ibin)/(corr_n*corr_n)));
    met_rmsfrac->SetBinContent(ibin, met_rms->GetBinContent(ibin)/met_mean->GetBinContent(ibin));
  }
  met_rms_ratio->Divide(met_rms, met_sqrtmean, 1., 1.);

  for (int ibin=1; ibin<=mt_corr->GetNbinsX(); ibin++) 
    for (int jbin=1; jbin<=mt_corr->GetNbinsY(); jbin++) {
      if (ibin == jbin) continue;
      mt_corr->SetBinContent(ibin, jbin, (mt_sum2->GetBinContent(ibin,jbin)/corr_n - mt_sum->GetBinContent(ibin)*mt_sum->GetBinContent(jbin)/(corr_n*corr_n))/(mt_rms->GetBinContent(ibin)*mt_rms->GetBinContent(jbin)));
    }
  for (int ibin=1; ibin<=elecpt_corr->GetNbinsX(); ibin++) 
    for (int jbin=1; jbin<=elecpt_corr->GetNbinsY(); jbin++) {
      if (ibin == jbin) continue;
      elecpt_corr->SetBinContent(ibin, jbin, (elecpt_sum2->GetBinContent(ibin,jbin)/corr_n - elecpt_sum->GetBinContent(ibin)*elecpt_sum->GetBinContent(jbin)/(corr_n*corr_n))/(elecpt_rms->GetBinContent(ibin)*elecpt_rms->GetBinContent(jbin)));
    }
  for (int ibin=1; ibin<=met_corr->GetNbinsX(); ibin++) 
    for (int jbin=1; jbin<=met_corr->GetNbinsY(); jbin++) {
      if (ibin == jbin) continue;
      met_corr->SetBinContent(ibin, jbin, (met_sum2->GetBinContent(ibin,jbin)/corr_n - met_sum->GetBinContent(ibin)*met_sum->GetBinContent(jbin)/(corr_n*corr_n))/(met_rms->GetBinContent(ibin)*met_rms->GetBinContent(jbin)));
    }

  for (int i=0; i<300/rebin; i++) {
    mt_dists[i] = new TH1D(TString::Format("mt_dists_%d", i), "", 25, mt_mean->GetBinContent(i)-5*mt_rms->GetBinContent(i), mt_mean->GetBinContent(i)+5*mt_rms->GetBinContent(i));
    mt_poisson[i] = new TH1D(TString::Format("mt_poisson_%d", i), "", 25, mt_mean->GetBinContent(i)-5*mt_rms->GetBinContent(i), mt_mean->GetBinContent(i)+5*mt_rms->GetBinContent(i));
  }
  for (int i=0; i<200/rebin; i++) {
    pt_dists[i] = new TH1D(TString::Format("pt_dists_%d", i), "", 25, elecpt_mean->GetBinContent(i)-5*elecpt_rms->GetBinContent(i), elecpt_mean->GetBinContent(i)+5*elecpt_rms->GetBinContent(i));
    pt_poisson[i] = new TH1D(TString::Format("pt_poisson_%d", i), "", 25, elecpt_mean->GetBinContent(i)-5*elecpt_rms->GetBinContent(i), elecpt_mean->GetBinContent(i)+5*elecpt_rms->GetBinContent(i));
  }
  for (int i=0; i<200/rebin; i++) {
    met_dists[i] = new TH1D(TString::Format("met_dists_%d", i), "", 25, met_mean->GetBinContent(i)-5*met_rms->GetBinContent(i), met_mean->GetBinContent(i)+5*met_rms->GetBinContent(i));
    met_poisson[i] = new TH1D(TString::Format("met_poisson_%d", i), "", 25, met_mean->GetBinContent(i)-5*met_rms->GetBinContent(i), met_mean->GetBinContent(i)+5*met_rms->GetBinContent(i));
  }
  ifstream _fileList2(argv[2]);
  while (!_fileList2.eof()) {
    _reader.ReadLine(_fileList2);
    if (_reader.Length() < 5) continue;
    TFile* tempFile = TFile::Open(_reader);
    if (!tempFile->IsOpen()) continue;
    TH1D* tempmt = (TH1D*) tempFile->Get("default/hWcandMt_CC");
    TH1D* tempelecpt = (TH1D*) tempFile->Get("default/hWcandElecPt_CC");
    TH1D* tempmet = (TH1D*) tempFile->Get("default/hWcandMet_CC");
    tempmt->Rebin(rebin);
    tempelecpt->Rebin(rebin);
    tempmet->Rebin(rebin);
    for (int i=0; i<300/rebin; i++)
      mt_dists[i]->Fill(tempmt->GetBinContent(i));
    for (int i=0; i<200/rebin; i++)
      pt_dists[i]->Fill(tempelecpt->GetBinContent(i));
    for (int i=0; i<200/rebin; i++)
      met_dists[i]->Fill(tempmet->GetBinContent(i));
  }
  for (int i=0; i<300/rebin; i++) {
    Double_t mt_dists_integral = mt_dists[i]->Integral();
    if (mt_dists_integral == 0) continue;
    mt_dists[i]->Sumw2();
    mt_dists[i]->Scale(1/mt_dists_integral);
    for (int j=1; j<=mt_poisson[i]->GetNbinsX(); j++) {
#ifdef doGaussian
      mt_poisson[i]->SetBinContent(j, TMath::Gaus(mt_poisson[i]->GetXaxis()->GetBinCenter(j), mt_mean->GetBinContent(i), 1.18411*TMath::Sqrt(mt_mean->GetBinContent(i))));
#else
      mt_poisson[i]->SetBinContent(j, TMath::PoissonI(mt_poisson[i]->GetXaxis()->GetBinCenter(j), mt_mean->GetBinContent(i)));
#endif
    }
  }
  for (int i=0; i<200/rebin; i++) {
    Double_t pt_dists_integral = pt_dists[i]->Integral();
    if (pt_dists_integral == 0) continue;
    pt_dists[i]->Sumw2();
    pt_dists[i]->Scale(1/pt_dists_integral);
    for (int j=1; j<=pt_poisson[i]->GetNbinsX(); j++) {
#ifdef doGaussian
      pt_poisson[i]->SetBinContent(j, TMath::Gaus(pt_poisson[i]->GetXaxis()->GetBinCenter(j), elecpt_mean->GetBinContent(i), 1.18378*TMath::Sqrt(elecpt_mean->GetBinContent(i))));
#else
      pt_poisson[i]->SetBinContent(j, TMath::PoissonI(pt_poisson[i]->GetXaxis()->GetBinCenter(j), elecpt_mean->GetBinContent(i)));
#endif
    }
  }
  for (int i=0; i<200/rebin; i++) {
    Double_t met_dists_integral = met_dists[i]->Integral();
    if (met_dists_integral == 0) continue;
    met_dists[i]->Sumw2();
    met_dists[i]->Scale(1/met_dists_integral);
    for (int j=1; j<=met_poisson[i]->GetNbinsX(); j++) {      
#ifdef doGaussian
      met_poisson[i]->SetBinContent(j, TMath::Gaus(met_poisson[i]->GetXaxis()->GetBinCenter(j), met_mean->GetBinContent(i), 1.18663*TMath::Sqrt(met_mean->GetBinContent(i))));
#else
      met_poisson[i]->SetBinContent(j, TMath::PoissonI(met_poisson[i]->GetXaxis()->GetBinCenter(j), met_mean->GetBinContent(i)));
#endif
    }
  }

  _fileList2.close();
  
  ifstream fileList(argv[1]);
  TString reader;
  Int_t pseudoNumber = 0;
  while (!fileList.eof()) {
    reader.ReadLine(fileList);
    if (reader.Length() < 5) continue;
    TFile* tempFile = TFile::Open(reader);
    if (!tempFile->IsOpen()) continue;
    TGraphErrors* mtGraph = (TGraphErrors*) tempFile->Get("wmass_Mt");
    TGraphErrors* ptGraph = (TGraphErrors*) tempFile->Get("wmass_Pt");
    TGraphErrors* metGraph = (TGraphErrors*) tempFile->Get("wmass_MET");
    mtptData->SetPoint(pseudoNumber, mtGraph->GetY()[0], ptGraph->GetY()[0]);
    metmtData->SetPoint(pseudoNumber, metGraph->GetY()[0], mtGraph->GetY()[0]);
    ptmetData->SetPoint(pseudoNumber, ptGraph->GetY()[0], metGraph->GetY()[0]);
    mtptScatter->SetPoint(pseudoNumber, mtGraph->GetY()[0], ptGraph->GetY()[0]);
    metmtScatter->SetPoint(pseudoNumber, metGraph->GetY()[0], mtGraph->GetY()[0]);
    ptmetScatter->SetPoint(pseudoNumber, ptGraph->GetY()[0], metGraph->GetY()[0]);
    mtptScatter->SetPointError(pseudoNumber, mtGraph->GetEY()[0], ptGraph->GetEY()[0]);
    metmtScatter->SetPointError(pseudoNumber, metGraph->GetEY()[0], mtGraph->GetEY()[0]);
    ptmetScatter->SetPointError(pseudoNumber, ptGraph->GetEY()[0], metGraph->GetEY()[0]);
    mtUnbiasedPull->Fill((mtGraph->GetY()[0]-trueMass)/mtGraph->GetEY()[0]);
    ptUnbiasedPull->Fill((ptGraph->GetY()[0]-trueMass)/ptGraph->GetEY()[0]);
    metUnbiasedPull->Fill((metGraph->GetY()[0]-trueMass)/metGraph->GetEY()[0]);
    mtptUnbiasedPull->Fill((mtGraph->GetY()[0]-trueMass)/mtGraph->GetEY()[0], (ptGraph->GetY()[0]-trueMass)/ptGraph->GetEY()[0]);
    ptmetUnbiasedPull->Fill((ptGraph->GetY()[0]-trueMass)/ptGraph->GetEY()[0], (metGraph->GetY()[0]-trueMass)/metGraph->GetEY()[0]);
    metmtUnbiasedPull->Fill((metGraph->GetY()[0]-trueMass)/metGraph->GetEY()[0], (mtGraph->GetY()[0]-trueMass)/mtGraph->GetEY()[0]);
    tempFile->Close();
    pseudoNumber++;
  }
  fileList.close();

  ifstream fileList2(argv[1]);
  while (!fileList2.eof()) {
    reader.ReadLine(fileList2);
    if (reader.Length() < 5) continue;
    TFile* tempFile = TFile::Open(reader);
    if (!tempFile->IsOpen()) continue;
    TGraphErrors* mtGraph = (TGraphErrors*) tempFile->Get("wmass_Mt");
    TGraphErrors* ptGraph = (TGraphErrors*) tempFile->Get("wmass_Pt");
    TGraphErrors* metGraph = (TGraphErrors*) tempFile->Get("wmass_MET");
    mtPull->Fill((mtGraph->GetY()[0]-mtptData->GetMean())/mtGraph->GetEY()[0]);
    ptPull->Fill((ptGraph->GetY()[0]-ptmetData->GetMean())/ptGraph->GetEY()[0]);
    metPull->Fill((metGraph->GetY()[0]-metmtData->GetMean())/metGraph->GetEY()[0]);
    mtptPull->Fill((mtGraph->GetY()[0]-mtptData->GetMean())/mtGraph->GetEY()[0], (ptGraph->GetY()[0]-ptmetData->GetMean())/ptGraph->GetEY()[0]);
    ptmetPull->Fill((ptGraph->GetY()[0]-ptmetData->GetMean())/ptGraph->GetEY()[0], (metGraph->GetY()[0]-metmtData->GetMean())/metGraph->GetEY()[0]);
    metmtPull->Fill((metGraph->GetY()[0]-metmtData->GetMean())/metGraph->GetEY()[0], (mtGraph->GetY()[0]-mtptData->GetMean())/mtGraph->GetEY()[0]);
    tempFile->Close();
  }
  fileList2.close();

  cout << mtptData->GetMean() << " " << mtptData->GetRMS() << endl;
  cout << ptmetData->GetMean() << " " << ptmetData->GetRMS() << endl;
  cout << metmtData->GetMean() << " " << metmtData->GetRMS() << endl;
  cout << mtptData->GetCorrelationFactor() << " " << ptmetData->GetCorrelationFactor() << " " << metmtData->GetCorrelationFactor() << endl;
  cout << "Data saved! " << pseudoNumber << endl;

  mtUnbiasedPull->Fit("gaus");
  ptUnbiasedPull->Fit("gaus");
  metUnbiasedPull->Fit("gaus");

  mtPull->Fit("gaus");
  ptPull->Fit("gaus");
  metPull->Fit("gaus");

  TF2* mtptUnbiasedGaussian2d = new TF2("mtptUnbiasedGaussian2d", gaussian2d, -5., 5., -5., 5., 6);
  TF2* ptmetUnbiasedGaussian2d = new TF2("ptmetUnbiasedGaussian2d", gaussian2d, -5., 5., -5., 5., 6);
  TF2* metmtUnbiasedGaussian2d = new TF2("metmtUnbiasedGaussian2d", gaussian2d, -5., 5., -5., 5., 6);

  mtptUnbiasedGaussian2d->SetParameter(0, 10.);
  mtptUnbiasedGaussian2d->SetParameter(3, 1.);
  mtptUnbiasedGaussian2d->SetParameter(4, 1.);
  mtptUnbiasedGaussian2d->SetParameter(5, 0.5);
  
  ptmetUnbiasedGaussian2d->SetParameter(0, 10.);
  ptmetUnbiasedGaussian2d->SetParameter(3, 1.);
  ptmetUnbiasedGaussian2d->SetParameter(4, 1.);
  ptmetUnbiasedGaussian2d->SetParameter(5, 0.5);

  metmtUnbiasedGaussian2d->SetParameter(0, 10.);
  metmtUnbiasedGaussian2d->SetParameter(3, 1.);
  metmtUnbiasedGaussian2d->SetParameter(4, 1.);
  metmtUnbiasedGaussian2d->SetParameter(5, 0.5);

  mtptUnbiasedPull->Fit(mtptUnbiasedGaussian2d, "R");
  ptmetUnbiasedPull->Fit(ptmetUnbiasedGaussian2d, "R");
  metmtUnbiasedPull->Fit(metmtUnbiasedGaussian2d, "R");

  ((TF2*) mtptUnbiasedPull->GetListOfFunctions()->At(0))->SetContour(4);
  ((TF2*) mtptUnbiasedPull->GetListOfFunctions()->At(0))->SetLineColor(kRed);
  ((TF2*) mtptUnbiasedPull->GetListOfFunctions()->At(0))->SetNpx(100);
  ((TF2*) mtptUnbiasedPull->GetListOfFunctions()->At(0))->SetNpy(100);

  ((TF2*) ptmetUnbiasedPull->GetListOfFunctions()->At(0))->SetContour(4);
  ((TF2*) ptmetUnbiasedPull->GetListOfFunctions()->At(0))->SetLineColor(kRed);
  ((TF2*) ptmetUnbiasedPull->GetListOfFunctions()->At(0))->SetNpx(100);
  ((TF2*) ptmetUnbiasedPull->GetListOfFunctions()->At(0))->SetNpy(100);

  ((TF2*) metmtUnbiasedPull->GetListOfFunctions()->At(0))->SetContour(4);
  ((TF2*) metmtUnbiasedPull->GetListOfFunctions()->At(0))->SetLineColor(kRed);
  ((TF2*) metmtUnbiasedPull->GetListOfFunctions()->At(0))->SetNpx(100);
  ((TF2*) metmtUnbiasedPull->GetListOfFunctions()->At(0))->SetNpy(100);

  mtptUnbiasedPull->SetLineWidth(2);
  ptmetUnbiasedPull->SetLineWidth(2);
  metmtUnbiasedPull->SetLineWidth(2);

  mtptUnbiasedPull->SetLineColor(kBlue);
  ptmetUnbiasedPull->SetLineColor(kBlue);
  metmtUnbiasedPull->SetLineColor(kBlue);

  ((TF1*) mtUnbiasedPull->GetListOfFunctions()->At(0))->SetLineColor(kRed);
  ((TF1*) ptUnbiasedPull->GetListOfFunctions()->At(0))->SetLineColor(kRed);
  ((TF1*) metUnbiasedPull->GetListOfFunctions()->At(0))->SetLineColor(kRed);

  TF2* mtptGaussian2d = new TF2("mtptGaussian2d", gaussian2d, -5., 5., -5., 5., 6);
  TF2* ptmetGaussian2d = new TF2("ptmetGaussian2d", gaussian2d, -5., 5., -5., 5., 6);
  TF2* metmtGaussian2d = new TF2("metmtGaussian2d", gaussian2d, -5., 5., -5., 5., 6);

  mtptGaussian2d->SetParameter(0, 10.);
  mtptGaussian2d->SetParameter(3, 1.);
  mtptGaussian2d->SetParameter(4, 1.);
  mtptGaussian2d->SetParameter(5, 0.5);
  
  ptmetGaussian2d->SetParameter(0, 10.);
  ptmetGaussian2d->SetParameter(3, 1.);
  ptmetGaussian2d->SetParameter(4, 1.);
  ptmetGaussian2d->SetParameter(5, 0.5);

  metmtGaussian2d->SetParameter(0, 10.);
  metmtGaussian2d->SetParameter(3, 1.);
  metmtGaussian2d->SetParameter(4, 1.);
  metmtGaussian2d->SetParameter(5, 0.5);

  mtptPull->Fit(mtptGaussian2d, "R");
  ptmetPull->Fit(ptmetGaussian2d, "R");
  metmtPull->Fit(metmtGaussian2d, "R");

  ((TF2*) mtptPull->GetListOfFunctions()->At(0))->SetContour(4);
  ((TF2*) mtptPull->GetListOfFunctions()->At(0))->SetLineColor(kRed);
  ((TF2*) mtptPull->GetListOfFunctions()->At(0))->SetNpx(100);
  ((TF2*) mtptPull->GetListOfFunctions()->At(0))->SetNpy(100);

  ((TF2*) ptmetPull->GetListOfFunctions()->At(0))->SetContour(4);
  ((TF2*) ptmetPull->GetListOfFunctions()->At(0))->SetLineColor(kRed);
  ((TF2*) ptmetPull->GetListOfFunctions()->At(0))->SetNpx(100);
  ((TF2*) ptmetPull->GetListOfFunctions()->At(0))->SetNpy(100);

  ((TF2*) metmtPull->GetListOfFunctions()->At(0))->SetContour(4);
  ((TF2*) metmtPull->GetListOfFunctions()->At(0))->SetLineColor(kRed);
  ((TF2*) metmtPull->GetListOfFunctions()->At(0))->SetNpx(100);
  ((TF2*) metmtPull->GetListOfFunctions()->At(0))->SetNpy(100);

  mtptPull->SetLineWidth(2);
  ptmetPull->SetLineWidth(2);
  metmtPull->SetLineWidth(2);

  mtptPull->SetLineColor(kBlue);
  ptmetPull->SetLineColor(kBlue);
  metmtPull->SetLineColor(kBlue);

  ((TF1*) mtPull->GetListOfFunctions()->At(0))->SetLineColor(kRed);
  ((TF1*) ptPull->GetListOfFunctions()->At(0))->SetLineColor(kRed);
  ((TF1*) metPull->GetListOfFunctions()->At(0))->SetLineColor(kRed);

  Double_t biasFactor = TMath::Sqrt(pseudoNumber/(pseudoNumber-1));
  Double_t biasFactor2 = TMath::Sqrt((pseudoNumber-1)/(pseudoNumber-2));

  // Do the JackKnife analysis
  cout << "Doing the jackknife analysis" << endl;
  for (Int_t i=0; i<pseudoNumber; i++) {
    TGraph* mtptTemp = (TGraph*) mtptData->Clone();
    TGraph* metmtTemp = (TGraph*) metmtData->Clone();
    TGraph* ptmetTemp = (TGraph*) ptmetData->Clone();
    mtptTemp->RemovePoint(i);
    metmtTemp->RemovePoint(i);
    ptmetTemp->RemovePoint(i);
    mtJKMean->Fill(mtptTemp->GetMean());
    ptJKMean->Fill(ptmetTemp->GetMean());
    metJKMean->Fill(metmtTemp->GetMean());
    mtJKSigma->Fill(biasFactor*pseudoNumber*mtptData->GetRMS()-biasFactor2*(pseudoNumber-1)*mtptTemp->GetRMS());
    ptJKSigma->Fill(biasFactor*pseudoNumber*ptmetData->GetRMS()-biasFactor2*(pseudoNumber-1)*ptmetTemp->GetRMS());
    metJKSigma->Fill(biasFactor*pseudoNumber*metmtData->GetRMS()-biasFactor2*(pseudoNumber-1)*metmtTemp->GetRMS());

    mtJKUnbiasedSigma->Fill(pseudoNumber*TMath::Sqrt(mtptData->GetRMS()*mtptData->GetRMS()+mtptData->GetMean()*mtptData->GetMean()-2*trueMass*mtptData->GetMean()+trueMass*trueMass)
			    -(pseudoNumber-1)*TMath::Sqrt(mtptTemp->GetRMS()*mtptTemp->GetRMS()+mtptTemp->GetMean()*mtptTemp->GetMean()-2*trueMass*mtptTemp->GetMean()+trueMass*trueMass));
    ptJKUnbiasedSigma->Fill(pseudoNumber*TMath::Sqrt(ptmetData->GetRMS()*ptmetData->GetRMS()+ptmetData->GetMean()*ptmetData->GetMean()-2*trueMass*ptmetData->GetMean()+trueMass*trueMass)
			    -(pseudoNumber-1)*TMath::Sqrt(ptmetTemp->GetRMS()*ptmetTemp->GetRMS()+ptmetTemp->GetMean()*ptmetTemp->GetMean()-2*trueMass*ptmetTemp->GetMean()+trueMass*trueMass));
    metJKUnbiasedSigma->Fill(pseudoNumber*TMath::Sqrt(metmtData->GetRMS()*metmtData->GetRMS()+metmtData->GetMean()*metmtData->GetMean()-2*trueMass*metmtData->GetMean()+trueMass*trueMass)
			    -(pseudoNumber-1)*TMath::Sqrt(metmtTemp->GetRMS()*metmtTemp->GetRMS()+metmtTemp->GetMean()*metmtTemp->GetMean()-2*trueMass*metmtTemp->GetMean()+trueMass*trueMass));

    mtptJKCorr->Fill(pseudoNumber*mtptData->GetCorrelationFactor()-(pseudoNumber-1)*mtptTemp->GetCorrelationFactor());
    metmtJKCorr->Fill(pseudoNumber*metmtData->GetCorrelationFactor()-(pseudoNumber-1)*metmtTemp->GetCorrelationFactor());
    ptmetJKCorr->Fill(pseudoNumber*ptmetData->GetCorrelationFactor()-(pseudoNumber-1)*ptmetTemp->GetCorrelationFactor());
  }

  // The values we are measuring
  TMatrixD correlationValue(3,3);
  TMatrixD correlationError(3,3);
  TVectorD uncertaintyValue(3);
  TVectorD uncertaintyError(3);
  TVectorD unbiasedUncertaintyValue(3);
  TVectorD unbiasedUncertaintyError(3);

  correlationValue[0][1] = mtptJKCorr->GetMean(); correlationError[0][1] = mtptJKCorr->GetRMS()/TMath::Sqrt(pseudoNumber-1.);
  correlationValue[0][2] = metmtJKCorr->GetMean(); correlationError[0][2] = metmtJKCorr->GetRMS()/TMath::Sqrt(pseudoNumber-1.);
  correlationValue[1][2] = ptmetJKCorr->GetMean(); correlationError[1][2] = ptmetJKCorr->GetRMS()/TMath::Sqrt(pseudoNumber-1.);
  
  correlationValue[1][0] = correlationValue[0][1]; correlationError[1][0] = correlationError[0][1];
  correlationValue[2][0] = correlationValue[0][2]; correlationError[2][0] = correlationError[0][2];
  correlationValue[2][1] = correlationValue[1][2]; correlationError[2][1] = correlationError[1][2];
    
  correlationValue[0][0] = 1.0; correlationError[0][0] = 0.0;
  correlationValue[1][1] = 1.0; correlationError[1][1] = 0.0;
  correlationValue[2][2] = 1.0; correlationError[2][2] = 0.0;

  uncertaintyValue[0] = mtJKSigma->GetMean(); uncertaintyError[0] = mtJKSigma->GetRMS()/TMath::Sqrt(pseudoNumber-1.);
  uncertaintyValue[1] = ptJKSigma->GetMean(); uncertaintyError[1] = ptJKSigma->GetRMS()/TMath::Sqrt(pseudoNumber-1.);
  uncertaintyValue[2] = metJKSigma->GetMean(); uncertaintyError[2] = metJKSigma->GetRMS()/TMath::Sqrt(pseudoNumber-1.);

  unbiasedUncertaintyValue[0] = mtJKUnbiasedSigma->GetMean(); unbiasedUncertaintyError[0] = mtJKUnbiasedSigma->GetRMS()/TMath::Sqrt(pseudoNumber-1.);
  unbiasedUncertaintyValue[1] = ptJKUnbiasedSigma->GetMean(); unbiasedUncertaintyError[1] = ptJKUnbiasedSigma->GetRMS()/TMath::Sqrt(pseudoNumber-1.);
  unbiasedUncertaintyValue[2] = metJKUnbiasedSigma->GetMean(); unbiasedUncertaintyError[2] = metJKUnbiasedSigma->GetRMS()/TMath::Sqrt(pseudoNumber-1.);

  cout << "The order is mT, pT, MET" << endl;

  cout << "The correlation matrix" << endl;
  correlationValue.Print();

  cout << "The uncertainty in the correlation matrix" << endl;
  correlationError.Print();

  cout << "The uncertainty values" << endl;
  uncertaintyValue.Print();

  cout << "The uncertainty in the uncertainty values" << endl;
  uncertaintyError.Print();

  cout << "The unbiased uncertainty values" << endl;
  unbiasedUncertaintyValue.Print();

  cout << "The uncertainty in the unbiased uncertainty values" << endl;
  unbiasedUncertaintyError.Print();

  cout << "Jackknife bias" << endl;
  cout << (pseudoNumber-1)*(mtJKMean->GetMean()-mtptData->GetMean()) << endl;
  cout << (pseudoNumber-1)*(ptJKMean->GetMean()-ptmetData->GetMean()) << endl;
  cout << (pseudoNumber-1)*(metJKMean->GetMean()-metmtData->GetMean()) << endl << endl;

  cout << "Unbiased Pull information" << endl;

  cout << " mT) Mean: " << mtUnbiasedPull->GetMean() << " +/- " << mtUnbiasedPull->GetMeanError() << "   RMS: " << mtUnbiasedPull->GetRMS() << " +/- " << mtUnbiasedPull->GetRMSError() << endl;
  cout << " pT) Mean: " << ptUnbiasedPull->GetMean() << " +/- " << ptUnbiasedPull->GetMeanError() << "   RMS: " << ptUnbiasedPull->GetRMS() << " +/- " << ptUnbiasedPull->GetRMSError() << endl;
  cout << "MET) Mean: " << metUnbiasedPull->GetMean() << " +/- " << metUnbiasedPull->GetMeanError() << "   RMS: " << metUnbiasedPull->GetRMS() << " +/- " << metUnbiasedPull->GetRMSError() << endl;

  cout << "Pull information" << endl;

  cout << " mT) Mean: " << mtPull->GetMean() << " +/- " << mtPull->GetMeanError() << "   RMS: " << mtPull->GetRMS() << " +/- " << mtPull->GetRMSError() << endl;
  cout << " pT) Mean: " << ptPull->GetMean() << " +/- " << ptPull->GetMeanError() << "   RMS: " << ptPull->GetRMS() << " +/- " << ptPull->GetRMSError() << endl;
  cout << "MET) Mean: " << metPull->GetMean() << " +/- " << metPull->GetMeanError() << "   RMS: " << metPull->GetRMS() << " +/- " << metPull->GetRMSError() << endl;

  TMatrixD mtptcov(2,2);
  mtptcov[0][0] = mtptGaussian2d->GetParameter(3)*mtptGaussian2d->GetParameter(3);
  mtptcov[0][1] = mtptGaussian2d->GetParameter(3)*mtptGaussian2d->GetParameter(4)*mtptGaussian2d->GetParameter(5);
  mtptcov[1][0] = mtptGaussian2d->GetParameter(3)*mtptGaussian2d->GetParameter(4)*mtptGaussian2d->GetParameter(5);  
  mtptcov[1][1] = mtptGaussian2d->GetParameter(4)*mtptGaussian2d->GetParameter(4);
  TMatrixD mtptEigenVectors(2,2);
  TVectorD mtptEigenValues(2);
  mtptEigenVectors = mtptcov.EigenVectors(mtptEigenValues);
  
  cout << "---- mT-pT ----" << endl;
  cout << "   " << mtptGaussian2d->GetParameter(1) << " +/- " << mtptGaussian2d->GetParError(1) << endl;
  cout << "   " << mtptGaussian2d->GetParameter(2) << " +/- " << mtptGaussian2d->GetParError(2) << endl;
  cout << "   " << mtptGaussian2d->GetParameter(3) << " +/- " << mtptGaussian2d->GetParError(3) << endl;
  cout << "   " << mtptGaussian2d->GetParameter(4) << " +/- " << mtptGaussian2d->GetParError(4) << endl;
  cout << "   " << mtptGaussian2d->GetParameter(5) << " +/- " << mtptGaussian2d->GetParError(5) << endl;
  mtptEigenVectors.Print();
  mtptEigenValues.Print();

  TMatrixD ptmetcov(2,2);
  ptmetcov[0][0] = ptmetGaussian2d->GetParameter(3)*ptmetGaussian2d->GetParameter(3);
  ptmetcov[0][1] = ptmetGaussian2d->GetParameter(3)*ptmetGaussian2d->GetParameter(4)*ptmetGaussian2d->GetParameter(5);
  ptmetcov[1][0] = ptmetGaussian2d->GetParameter(3)*ptmetGaussian2d->GetParameter(4)*ptmetGaussian2d->GetParameter(5);  
  ptmetcov[1][1] = ptmetGaussian2d->GetParameter(4)*ptmetGaussian2d->GetParameter(4);
  TMatrixD ptmetEigenVectors(2,2);
  TVectorD ptmetEigenValues(2);
  ptmetEigenVectors = ptmetcov.EigenVectors(ptmetEigenValues);
  
  cout << "---- pT-MET ----" << endl;
  cout << "   " << ptmetGaussian2d->GetParameter(1) << " +/- " << ptmetGaussian2d->GetParError(1) << endl;
  cout << "   " << ptmetGaussian2d->GetParameter(2) << " +/- " << ptmetGaussian2d->GetParError(2) << endl;
  cout << "   " << ptmetGaussian2d->GetParameter(3) << " +/- " << ptmetGaussian2d->GetParError(3) << endl;
  cout << "   " << ptmetGaussian2d->GetParameter(4) << " +/- " << ptmetGaussian2d->GetParError(4) << endl;
  cout << "   " << ptmetGaussian2d->GetParameter(5) << " +/- " << ptmetGaussian2d->GetParError(5) << endl;
  ptmetEigenVectors.Print();
  ptmetEigenValues.Print();

  TMatrixD metmtcov(2,2);
  metmtcov[0][0] = metmtGaussian2d->GetParameter(3)*metmtGaussian2d->GetParameter(3);
  metmtcov[0][1] = metmtGaussian2d->GetParameter(3)*metmtGaussian2d->GetParameter(4)*metmtGaussian2d->GetParameter(5);
  metmtcov[1][0] = metmtGaussian2d->GetParameter(3)*metmtGaussian2d->GetParameter(4)*metmtGaussian2d->GetParameter(5);  
  metmtcov[1][1] = metmtGaussian2d->GetParameter(4)*metmtGaussian2d->GetParameter(4);
  TMatrixD metmtEigenVectors(2,2);
  TVectorD metmtEigenValues(2);
  metmtEigenVectors = metmtcov.EigenVectors(metmtEigenValues);
  
  cout << "---- MET-mT ----" << endl;
  cout << "   " << metmtGaussian2d->GetParameter(1) << " +/- " << metmtGaussian2d->GetParError(1) << endl;
  cout << "   " << metmtGaussian2d->GetParameter(2) << " +/- " << metmtGaussian2d->GetParError(2) << endl;
  cout << "   " << metmtGaussian2d->GetParameter(3) << " +/- " << metmtGaussian2d->GetParError(3) << endl;
  cout << "   " << metmtGaussian2d->GetParameter(4) << " +/- " << metmtGaussian2d->GetParError(4) << endl;
  cout << "   " << metmtGaussian2d->GetParameter(5) << " +/- " << metmtGaussian2d->GetParError(5) << endl;
  metmtEigenVectors.Print();
  metmtEigenValues.Print();

  TFile* outputFile = TFile::Open("StatUnc.root", "recreate");
  mt_mean->Write();
  elecpt_mean->Write();
  met_mean->Write();

  mt_rms->Write();
  elecpt_rms->Write();
  met_rms->Write();

  mt_rms_ratio->Write();
  elecpt_rms_ratio->Write();
  met_rms_ratio->Write();

  mt_rmsfrac->Write();
  elecpt_rmsfrac->Write();
  met_rmsfrac->Write();

  mt_corr->Write();
  elecpt_corr->Write();
  met_corr->Write();

  mtPull->Write();
  ptPull->Write();
  metPull->Write();

  mtptPull->Write();
  ptmetPull->Write();
  metmtPull->Write();

  mtUnbiasedPull->Write();
  ptUnbiasedPull->Write();
  metUnbiasedPull->Write();

  mtptUnbiasedPull->Write();
  ptmetUnbiasedPull->Write();
  metmtUnbiasedPull->Write();
  
  mtJKMean->Write();
  ptJKMean->Write();
  metJKMean->Write();

  mtJKSigma->Write();
  ptJKSigma->Write();
  metJKSigma->Write();

  mtJKUnbiasedSigma->Write();
  ptJKUnbiasedSigma->Write();
  metJKUnbiasedSigma->Write();

  mtptJKCorr->Write();
  metmtJKCorr->Write();
  ptmetJKCorr->Write();

  for (int i=0; i<300/rebin; i++)
    mt_dists[i]->Write();
  for (int i=0; i<200/rebin; i++)
    pt_dists[i]->Write();
  for (int i=0; i<200/rebin; i++)
    met_dists[i]->Write();

  outputFile->Close();

  TCanvas* c1 = new TCanvas();
  c1->Print("StatUnc.ps[");

  {
    for (int i=1; i<=mt_mean->GetNbinsX(); i++) {
      Double_t expectedUncertainty = mt_mean->GetBinContent(i);
      expectedUncertainty = TMath::Sqrt(expectedUncertainty);
      mt_mean->SetBinContent(i, expectedUncertainty);
      mt_mean->SetBinError(i, TMath::Sqrt(expectedUncertainty));
    }
    mt_mean->SetLineColor(kBlue);
    mt_rms->SetLineColor(kRed);
    mt_rms->Sumw2();
    mt_mean->Sumw2();
    double __sum = 0;
    double __n = 0;
    for (int i=1; i<=mt_rms_ratio->GetNbinsX(); i++) {
      if (mt_rms_ratio->GetBinContent(i) > 0) {
	__sum += mt_rms_ratio->GetBinContent(i);
	__n += 1;
      }
    }
    cout << "Mt correction " << __sum/__n << endl;

    TPad *pad1 = new TPad("pad1","pad1",0,0.3,1,1);
    pad1->SetBottomMargin(0);
    pad1->Draw();
    pad1->cd();
    mt_rms->SetStats(0);
    mt_mean->SetStats(0);
    mt_rms->SetMarkerStyle(20);
    mt_mean->SetMarkerStyle(21);
    mt_rms->SetMinimum(-10);
    mt_rms->GetXaxis()->SetRangeUser(0, 120);
    mt_rms->Draw();
    mt_mean->Draw("same");

    TLegend* leg = new TLegend(0.6, 0.6, 0.8, 0.8);
    leg->AddEntry(mt_rms, "RMS", "lp");
    leg->AddEntry(mt_mean, "#sqrt{mean}", "lp");
    leg->SetFillColor(0);
    leg->SetBorderSize(0);
    leg->Draw();

    c1->cd();
    TPad *pad2 = new TPad("pad2","pad2",0,0.,1,0.35);
    pad2->SetTopMargin(0);
    pad2->SetBottomMargin(0.5);
    pad2->Draw();
    pad2->cd();
    mt_rms_ratio->SetTitle(" ;m_{T} (GeV); ");
    mt_rms_ratio->SetLineColor(kBlack);
    mt_rms_ratio->SetMarkerColor(kBlue);
    mt_rms_ratio->SetMarkerStyle(20);
    mt_rms_ratio->SetStats(0);
    mt_rms_ratio->GetXaxis()->SetLabelFont(63);
    mt_rms_ratio->GetXaxis()->SetLabelSize(16);
    mt_rms_ratio->GetYaxis()->SetLabelFont(63);
    mt_rms_ratio->GetYaxis()->SetLabelSize(12);
    mt_rms_ratio->GetYaxis()->SetRangeUser(0.45, 1.55);
    mt_rms_ratio->GetYaxis()->SetNdivisions(5);
    mt_rms_ratio->GetXaxis()->SetRangeUser(0, 120);
    mt_rms_ratio->GetXaxis()->SetTitleFont(63);
    mt_rms_ratio->GetXaxis()->SetTitleSize(16);
    mt_rms_ratio->GetXaxis()->SetTitleOffset(2.5);
    mt_rms_ratio->Draw("e0");
  }
  c1->Print("StatUnc.ps");

  c1->Clear();
  c1->cd();
  {
    for (int i=1; i<=met_mean->GetNbinsX(); i++) {
      Double_t expectedUncertainty = met_mean->GetBinContent(i);
      expectedUncertainty = TMath::Sqrt(expectedUncertainty);
      met_mean->SetBinContent(i, expectedUncertainty);
      met_mean->SetBinError(i, TMath::Sqrt(expectedUncertainty));
    }
    met_mean->SetLineColor(kBlue);
    met_rms->SetLineColor(kRed);
    met_rms->Sumw2();
    met_mean->Sumw2();
    double __sum = 0;
    double __n = 0;
    for (int i=1; i<=met_rms_ratio->GetNbinsX(); i++) {
      if (met_rms_ratio->GetBinContent(i) > 0) {
	__sum += met_rms_ratio->GetBinContent(i);
	__n += 1;
      }
    }
    cout << "Met correction " << __sum/__n << endl;

    TPad *pad1 = new TPad("pad1","pad1",0,0.3,1,1);
    pad1->SetBottomMargin(0);
    pad1->Draw();
    pad1->cd();
    met_rms->SetStats(0);
    met_mean->SetStats(0);
    met_rms->SetMarkerStyle(20);
    met_mean->SetMarkerStyle(21);
    met_rms->SetMinimum(-10);
    met_rms->GetXaxis()->SetRangeUser(25, 60);
    met_rms->Draw();
    met_mean->Draw("same");

    TLegend* leg = new TLegend(0.6, 0.6, 0.8, 0.8);
    leg->AddEntry(met_rms, "RMS", "lp");
    leg->AddEntry(met_mean, "#sqrt{mean}", "lp");
    leg->SetFillColor(0);
    leg->SetBorderSize(0);
    leg->Draw();

    c1->cd();
    TPad *pad2 = new TPad("pad2","pad2",0,0,1,0.3);
    pad2->SetTopMargin(0);
    pad2->SetBottomMargin(0.5);
    pad2->Draw();
    pad2->cd();
    met_rms_ratio->SetTitle(" ;MET (GeV); ");
    met_rms_ratio->SetLineColor(kBlack);
    met_rms_ratio->SetMarkerColor(kBlue);
    met_rms_ratio->SetMarkerStyle(20);
    met_rms_ratio->SetStats(0);
    met_rms_ratio->GetXaxis()->SetLabelFont(63);
    met_rms_ratio->GetXaxis()->SetLabelSize(16);
    met_rms_ratio->GetYaxis()->SetLabelFont(63);
    met_rms_ratio->GetYaxis()->SetLabelSize(12);
    met_rms_ratio->GetYaxis()->SetRangeUser(0.45, 1.55);
    met_rms_ratio->GetYaxis()->SetNdivisions(5);
    met_rms_ratio->GetXaxis()->SetRangeUser(25, 60);
    met_rms_ratio->GetXaxis()->SetTitleFont(63);
    met_rms_ratio->GetXaxis()->SetTitleSize(16);
    met_rms_ratio->GetXaxis()->SetTitleOffset(2.5);
    met_rms_ratio->Draw("e0");
  }
  c1->Print("StatUnc.ps");

  c1->Clear();
  c1->cd();
  {
    for (int i=1; i<=elecpt_mean->GetNbinsX(); i++) {
      Double_t expectedUncertainty = elecpt_mean->GetBinContent(i);
      expectedUncertainty = TMath::Sqrt(expectedUncertainty);
      elecpt_mean->SetBinContent(i, expectedUncertainty);
      elecpt_mean->SetBinError(i, TMath::Sqrt(expectedUncertainty));
    }
    elecpt_mean->SetLineColor(kBlue);
    elecpt_rms->SetLineColor(kRed);
    elecpt_rms->Sumw2();
    elecpt_mean->Sumw2();
    double __sum = 0;
    double __n = 0;
    for (int i=1; i<=elecpt_rms_ratio->GetNbinsX(); i++) {
      if (elecpt_rms_ratio->GetBinContent(i) > 0) {
	__sum += elecpt_rms_ratio->GetBinContent(i);
	__n += 1;
      }
    }
    cout << "Elecpt correction " << __sum/__n << endl;

    TPad *pad1 = new TPad("pad1","pad1",0,0.3,1,1);
    pad1->SetBottomMargin(0);
    pad1->Draw();
    pad1->cd();
    elecpt_rms->SetStats(0);
    elecpt_mean->SetStats(0);
    elecpt_rms->SetMarkerStyle(20);
    elecpt_mean->SetMarkerStyle(21);
    elecpt_rms->SetMinimum(-10);
    elecpt_rms->GetXaxis()->SetRangeUser(25, 60);
    elecpt_rms->Draw();
    elecpt_mean->Draw("same");

    TLegend* leg = new TLegend(0.6, 0.6, 0.8, 0.8);
    leg->AddEntry(elecpt_rms, "RMS", "lp");
    leg->AddEntry(elecpt_mean, "#sqrt{mean}", "lp");
    leg->SetFillColor(0);
    leg->SetBorderSize(0);
    leg->Draw();

    c1->cd();
    TPad *pad2 = new TPad("pad2","pad2",0,0,1,0.3);
    pad2->SetTopMargin(0);
    pad2->SetBottomMargin(0.5);
    pad2->Draw();
    pad2->cd();
    elecpt_rms_ratio->SetTitle(" ; Electron p_{T} (GeV); ");
    elecpt_rms_ratio->SetLineColor(kBlack);
    elecpt_rms_ratio->SetMarkerColor(kBlue);
    elecpt_rms_ratio->SetMarkerStyle(20);
    elecpt_rms_ratio->SetStats(0);
    elecpt_rms_ratio->GetXaxis()->SetLabelFont(63);
    elecpt_rms_ratio->GetXaxis()->SetLabelSize(16);
    elecpt_rms_ratio->GetYaxis()->SetLabelFont(63);
    elecpt_rms_ratio->GetYaxis()->SetLabelSize(12);
    elecpt_rms_ratio->GetYaxis()->SetRangeUser(0.45, 1.55);
    elecpt_rms_ratio->GetYaxis()->SetNdivisions(5);
    elecpt_rms_ratio->GetXaxis()->SetRangeUser(25, 60);
    elecpt_rms_ratio->GetXaxis()->SetTitleFont(63);
    elecpt_rms_ratio->GetXaxis()->SetTitleSize(16);
    elecpt_rms_ratio->GetXaxis()->SetTitleOffset(2.5);
    elecpt_rms_ratio->Draw("e0");
  }
  c1->Print("StatUnc.ps");

  c1->Clear();
  c1->cd();
  mt_corr->SetTitle("Bin-by-bin correlation m_{T}; m_{T} (GeV); m_{T} (GeV)");
  mt_corr->GetXaxis()->SetRangeUser(65, 90);
  mt_corr->GetYaxis()->SetRangeUser(65, 90);
  mt_corr->GetZaxis()->SetRangeUser(-0.4, 0.4);
  mt_corr->SetStats(0);
  mt_corr->Draw("colz");
  c1->Print("StatUnc.ps");

  elecpt_corr->SetTitle("Bin-by-bin correlation Electron p_{T}; Electron p_{T} (GeV); Electron p_{T} (GeV)");
  elecpt_corr->GetXaxis()->SetRangeUser(32, 48);
  elecpt_corr->GetYaxis()->SetRangeUser(32, 48);
  elecpt_corr->GetZaxis()->SetRangeUser(-0.4, 0.4);
  elecpt_corr->SetStats(0);
  elecpt_corr->Draw("colz");
  c1->Print("StatUnc.ps");  

  met_corr->SetTitle("Bin-by-bin correlation MET; MET (GeV); MET (GeV)");
  met_corr->GetXaxis()->SetRangeUser(32, 48);
  met_corr->GetYaxis()->SetRangeUser(32, 48);
  met_corr->GetZaxis()->SetRangeUser(-0.4, 0.4);
  met_corr->SetStats(0);
  met_corr->Draw("colz");
  c1->Print("StatUnc.ps");  

  mtUnbiasedPull->SetTitle("m_{T} unbiased pull; Unbiased pull; Pseudo-experiments");
  mtUnbiasedPull->Draw();
  c1->Print("StatUnc.ps");

  ptUnbiasedPull->SetTitle("p_{T} unbiased pull; Unbiased pull; Pseudo-experiments");
  ptUnbiasedPull->Draw();
  c1->Print("StatUnc.ps");

  metUnbiasedPull->SetTitle("MET unbiased pull; Unbiased pull; Pseudo-experiments");
  metUnbiasedPull->Draw();
  c1->Print("StatUnc.ps");

  mtptUnbiasedPull->SetTitle("m_{T} - p_{T} unbiased pull; m_{T} pull; p_{T} pull");
  mtptUnbiasedPull->Draw("box");
  c1->Print("StatUnc.ps");

  ptmetUnbiasedPull->SetTitle("p_{T} - MET unbiased pull; p_{T} pull; MET pull");
  ptmetUnbiasedPull->Draw("box");
  c1->Print("StatUnc.ps");

  metmtUnbiasedPull->SetTitle("MET - m_{T} unbiased pull; MET pull; m_{T} pull");
  metmtUnbiasedPull->Draw("box");
  c1->Print("StatUnc.ps");

  mtPull->SetTitle("m_{T} pull; pull; Pseudo-experiments");
  mtPull->Draw();
  c1->Print("StatUnc.ps");

  ptPull->SetTitle("p_{T} pull; pull; Pseudo-experiments");
  ptPull->Draw();
  c1->Print("StatUnc.ps");

  metPull->SetTitle("MET pull; pull; Pseudo-experiments");
  metPull->Draw();
  c1->Print("StatUnc.ps");

  mtptPull->SetTitle("m_{T} - p_{T} pull; m_{T} pull; p_{T} pull");
  mtptPull->Draw("box");
  c1->Print("StatUnc.ps");

  ptmetPull->SetTitle("p_{T} - MET pull; p_{T} pull; MET pull");
  ptmetPull->Draw("box");
  c1->Print("StatUnc.ps");

  metmtPull->SetTitle("MET - m_{T} pull; MET pull; m_{T} pull");
  metmtPull->Draw("box");
  c1->Print("StatUnc.ps");

  mtJKSigma->SetTitle("m_{T} Standard Deviation Jackknife Estimator; Psuedo-value (GeV); Pseudo-experiment");
  mtJKSigma->Draw();
  c1->Print("StatUnc.ps");

  ptJKSigma->SetTitle("p_{T} Standard Deviation Jackknife Estimator; Psuedo-value (GeV); Pseudo-experiment");
  ptJKSigma->Draw();
  c1->Print("StatUnc.ps");

  metJKSigma->SetTitle("MET Standard Deviation Jackknife Estimator; Psuedo-value (GeV); Pseudo-experiment");
  metJKSigma->Draw();
  c1->Print("StatUnc.ps");

  mtJKUnbiasedSigma->SetTitle("m_{T} Unbiased Standard Deviation Jackknife Estimator; Psuedo-value (GeV); Pseudo-experiment");
  mtJKUnbiasedSigma->Draw();
  c1->Print("StatUnc.ps");

  ptJKUnbiasedSigma->SetTitle("p_{T} Unbiased Standard Deviation Jackknife Estimator; Psuedo-value (GeV); Pseudo-experiment");
  ptJKUnbiasedSigma->Draw();
  c1->Print("StatUnc.ps");

  metJKUnbiasedSigma->SetTitle("MET Unbiased Standard Deviation Jackknife Estimator; Psuedo-value (GeV); Pseudo-experiment");
  metJKUnbiasedSigma->Draw();
  c1->Print("StatUnc.ps");

  mtptJKCorr->SetTitle("m_{T} - p_{T} Correlation Jackknife Estimator; Pseudo-value; Pseudo-experiment");
  mtptJKCorr->Draw();
  c1->Print("StatUnc.ps");

  ptmetJKCorr->SetTitle("p_{T} - MET Correlation Jackknife Estimator; Pseudo-value; Pseudo-experiment");
  ptmetJKCorr->Draw();
  c1->Print("StatUnc.ps");

  metmtJKCorr->SetTitle("MET - m_{T} Correlation Jackknife Estimator; Pseudo-value; Pseudo-experiment");
  metmtJKCorr->Draw();
  c1->Print("StatUnc.ps");

  mtptScatter->Draw("A");
  mtptScatter->SetLineColor(kRed);
  mtptScatter->SetTitle("m_{T} - p_{T} Fitted W mass; m_{T} fitted W mass (GeV); p_{T} fitted W mass (GeV)");
  mtptScatter->Draw("AP");  
  c1->Print("StatUnc.ps");

  ptmetScatter->Draw("A");
  ptmetScatter->SetLineColor(kRed);
  ptmetScatter->SetTitle("p_{T} - MET Fitted W mass; p_{T} fitted W mass (GeV); MET fitted W mass (GeV)");
  ptmetScatter->Draw("AP");  
  c1->Print("StatUnc.ps");

  metmtScatter->Draw("A");
  metmtScatter->SetLineColor(kRed);
  metmtScatter->SetTitle("MET - m_{T} Fitted W mass; MET fitted W mass (GeV); m_{T} fitted W mass (GeV)");
  metmtScatter->Draw("AP");  
  c1->Print("StatUnc.ps");

  c1->Clear();
  c1->Divide(2,2);
  Int_t __i = 1;
  for (int i=0; i<300/rebin; i++) {
    c1->cd(__i);
    Double_t mt_poisson_integral = mt_poisson[i]->Integral();
    mt_poisson[i]->Sumw2();
    mt_poisson[i]->Scale(1/mt_poisson_integral);
    mt_poisson[i]->SetLineColor(kRed);
    mt_poisson[i]->SetTitle(TString::Format("m_{T} bin [%.1f, %.1f]", mt_mean->GetXaxis()->GetBinLowEdge(i+1), mt_mean->GetXaxis()->GetBinUpEdge(i+1)));
    mt_poisson[i]->SetStats(0);    
    mt_poisson[i]->Draw("hist");
    mt_dists[i]->SetLineColor(kBlue);
    mt_dists[i]->Draw("same, hist");
    __i++;
    if (__i == 5) {
      c1->Print("StatUnc.ps");
      __i = 1;
    }
  }

    
  c1->Print("StatUnc.ps]");
}

