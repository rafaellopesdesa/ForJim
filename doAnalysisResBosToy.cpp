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
#define rebin 10

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
  TH1D* mt_poisson[300/rebin];

  TH1D* mtUnbiasedPull = new TH1D("mtUnbiasedPull", "", 50, -5., 5.);
  TH1D* mtPull = new TH1D("mtPull", "", 50, -5., 5.);
  TH1D* mtJKMean = new TH1D("mtJKMean", "", 50, 80.380, 80.440);
  TH1D* mtJKSigma = new TH1D("mtJKSigma", "", 50, 0.00, 0.01);
  TH1D* mtJKUnbiasedSigma = new TH1D("mtJKUnbiasedSigma", "", 50, 0.00, 0.01);

  // Storage for the data
  TGraph* mtptData = new TGraph();
  TGraphErrors* mtptScatter = new TGraphErrors();

  // Reads the data
  cout << "Getting the data" << endl;
  TH1D* mt_sum = new TH1D("mt_sum", "", 300/rebin, 50., 200.);
  TH1D* mt_mean = new TH1D("mt_mean", "", 300/rebin, 50., 200.);
  TH1D* mt_rms = new TH1D("mt_rms", "", 300/rebin, 50., 200.);
  TH1D* mt_rmsfrac = new TH1D("mt_rmsfrac", "", 300/rebin, 50., 200.);
  TH2D* mt_sum2 = new TH2D("mt_sum2", "", 300/rebin, 50., 200., 300/rebin, 50., 200.);
  TH2D* mt_corr = new TH2D("mt_corr", "", 300/rebin, 50., 200., 300/rebin, 50., 200.);
  Double_t corr_n = 0;

  ifstream _fileList(argv[2]);
  TString _reader;
  while (!_fileList.eof()) {
    _reader.ReadLine(_fileList);
    if (_reader.Length() < 5) continue;
    TFile* tempFile = TFile::Open(_reader);
    if (!tempFile->IsOpen()) continue;
    TH1D* tempmt = (TH1D*) tempFile->Get("default/hWcandMt_CC");
    tempmt->Rebin(rebin);
    for (int ibin=1; ibin<=tempmt->GetNbinsX(); ibin++) {
      mt_sum->Fill(tempmt->GetXaxis()->GetBinCenter(ibin), tempmt->GetBinContent(ibin));
      for (int jbin=1; jbin<=tempmt->GetNbinsX(); jbin++) {
	mt_sum2->Fill(tempmt->GetXaxis()->GetBinCenter(ibin), tempmt->GetXaxis()->GetBinCenter(jbin), tempmt->GetBinContent(ibin)*tempmt->GetBinContent(jbin));
      }
    }
    corr_n++;
    tempFile->Close();
  }
  _fileList.close();
  for (int ibin=1; ibin<=mt_rms->GetNbinsX(); ibin++) {
    mt_mean->SetBinContent(ibin, mt_sum->GetBinContent(ibin)/corr_n);
    mt_rms->SetBinContent(ibin, TMath::Sqrt(mt_sum2->GetBinContent(ibin, ibin)/corr_n-mt_sum->GetBinContent(ibin)*mt_sum->GetBinContent(ibin)/(corr_n*corr_n)));
    mt_rmsfrac->SetBinContent(ibin, mt_rms->GetBinContent(ibin)/mt_mean->GetBinContent(ibin));
  }
  for (int ibin=1; ibin<=mt_corr->GetNbinsX(); ibin++) 
    for (int jbin=1; jbin<=mt_corr->GetNbinsY(); jbin++) {
      if (ibin == jbin) continue;
      mt_corr->SetBinContent(ibin, jbin, (mt_sum2->GetBinContent(ibin,jbin)/corr_n - mt_sum->GetBinContent(ibin)*mt_sum->GetBinContent(jbin)/(corr_n*corr_n))/(mt_rms->GetBinContent(ibin)*mt_rms->GetBinContent(jbin)));
    }

  for (int i=0; i<300/rebin; i++) {
    mt_dists[i] = new TH1D(TString::Format("mt_dists_%d", i), "", 25, mt_mean->GetBinContent(i)-5*mt_rms->GetBinContent(i), mt_mean->GetBinContent(i)+5*mt_rms->GetBinContent(i));
    mt_poisson[i] = new TH1D(TString::Format("mt_poisson_%d", i), "", 25, mt_mean->GetBinContent(i)-5*mt_rms->GetBinContent(i), mt_mean->GetBinContent(i)+5*mt_rms->GetBinContent(i));
  }
  ifstream _fileList2(argv[2]);
  while (!_fileList2.eof()) {
    _reader.ReadLine(_fileList2);
    if (_reader.Length() < 5) continue;
    TFile* tempFile = TFile::Open(_reader);
    if (!tempFile->IsOpen()) continue;
    TH1D* tempmt = (TH1D*) tempFile->Get("default/hWcandMt_CC");
    tempmt->Rebin(rebin)
    for (int i=0; i<300/rebin; i++)
      mt_dists[i]->Fill(tempmt->GetBinContent(i));
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
    mtptData->SetPoint(pseudoNumber, mtGraph->GetY()[0], pseudoNumber);
    mtptScatter->SetPoint(pseudoNumber, mtGraph->GetY()[0], pseudoNumber);
    mtptScatter->SetPointError(pseudoNumber, mtGraph->GetEY()[0], 0.);
    mtUnbiasedPull->Fill((mtGraph->GetY()[0]-trueMass)/mtGraph->GetEY()[0]);
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
    mtPull->Fill((mtGraph->GetY()[0]-mtptData->GetMean())/mtGraph->GetEY()[0]);
    tempFile->Close();
  }
  fileList2.close();

  mtUnbiasedPull->Fit("gaus");
  mtPull->Fit("gaus");

  ((TF1*) mtUnbiasedPull->GetListOfFunctions()->At(0))->SetLineColor(kRed);
  ((TF1*) mtPull->GetListOfFunctions()->At(0))->SetLineColor(kRed);


  Double_t biasFactor = TMath::Sqrt(pseudoNumber/(pseudoNumber-1));
  Double_t biasFactor2 = TMath::Sqrt((pseudoNumber-1)/(pseudoNumber-2));

  // Do the JackKnife analysis
  cout << "Doing the jackknife analysis" << endl;
  for (Int_t i=0; i<pseudoNumber; i++) {
    TGraph* mtptTemp = (TGraph*) mtptData->Clone();
    mtptTemp->RemovePoint(i);
    mtJKMean->Fill(mtptTemp->GetMean());
    mtJKSigma->Fill(biasFactor*pseudoNumber*mtptData->GetRMS()-biasFactor2*(pseudoNumber-1)*mtptTemp->GetRMS());
    mtJKUnbiasedSigma->Fill(pseudoNumber*TMath::Sqrt(mtptData->GetRMS()*mtptData->GetRMS()+mtptData->GetMean()*mtptData->GetMean()-2*trueMass*mtptData->GetMean()+trueMass*trueMass)
			    -(pseudoNumber-1)*TMath::Sqrt(mtptTemp->GetRMS()*mtptTemp->GetRMS()+mtptTemp->GetMean()*mtptTemp->GetMean()-2*trueMass*mtptTemp->GetMean()+trueMass*trueMass));
  }

  TFile* outputFile = TFile::Open("StatUnc.root", "recreate");
  mt_mean->Write();
  mt_rms->Write();
  mt_rmsfrac->Write();
  mtPull->Write();
  mtUnbiasedPull->Write();
  mtJKMean->Write();
  mtJKSigma->Write();
  mtJKUnbiasedSigma->Write();

  for (int i=0; i<300/rebin; i++)
    mt_dists[i]->Write();

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
    TH1D* mt_rms_ratio = (TH1D*) mt_mean->Clone();
    mt_rms_ratio->Divide(mt_rms, mt_mean, 1., 1.);
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
  mt_corr->SetTitle("Bin-by-bin correlation m_{T}; m_{T} (GeV); m_{T} (GeV)");
  mt_corr->GetXaxis()->SetRangeUser(65, 90);
  mt_corr->GetYaxis()->SetRangeUser(65, 90);
  mt_corr->GetZaxis()->SetRangeUser(-0.2, 0.2);
  mt_corr->SetStats(0);
  mt_corr->Draw("colz");
  c1->Print("StatUnc.ps");

  mtUnbiasedPull->SetTitle("m_{T} unbiased pull; Unbiased pull; Pseudo-experiments");
  mtUnbiasedPull->Draw();
  c1->Print("StatUnc.ps");

  mtPull->SetTitle("m_{T} pull; pull; Pseudo-experiments");
  mtPull->Draw();
  c1->Print("StatUnc.ps");

  mtJKSigma->SetTitle("m_{T} Standard Deviation Jackknife Estimator; Psuedo-value (GeV); Pseudo-experiment");
  mtJKSigma->Draw();
  c1->Print("StatUnc.ps");

  mtJKUnbiasedSigma->SetTitle("m_{T} Unbiased Standard Deviation Jackknife Estimator; Psuedo-value (GeV); Pseudo-experiment");
  mtJKUnbiasedSigma->Draw();
  c1->Print("StatUnc.ps");

  mtptScatter->Draw("A");
  mtptScatter->SetLineColor(kRed);
  mtptScatter->SetTitle("m_{T} - p_{T} Fitted W mass; m_{T} fitted W mass (GeV); Pseudo-Experiment");
  mtptScatter->Draw("AP");  
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

