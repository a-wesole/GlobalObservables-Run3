#ifndef utils
#define utils
#include <TGraphAsymmErrors.h>
#include <TGraphErrors.h>
#include <TGraph.h>
#include <TF1.h>
#include <TLegend.h>
#include <TCanvas.h>
#include <TLatex.h>
#include <TLine.h>
#include <TBox.h>
#include <TCut.h>
#include <TH1F.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH2F.h>
#include <TGaxis.h>
#include <TDatime.h>
#include <TMath.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1.h>
#include <TStyle.h>
#include <TEventList.h>
#include "TROOT.h"
#include "TChain.h"
#include "TStopwatch.h"

#include <iostream>     // std::cout
#include <ctime>        // std::clock()
#include <algorithm>    // std::find()
#include <iomanip>      // std::setprecision()
#include <vector>
#include <fstream>
#include <stdlib.h>
#include <string.h>

using namespace std;

const int col[] = {1,2,3,4,6,7,28,46,41};
const int ycol[] = {8,9,28,46,41};
const int marker[] = {24,25,26,27,28,29,31,33,34};


void canvasStyle(TCanvas*& canv,
		 const Float_t leftMargin=0.17,
		 const Float_t bottomMargin=0.12,
		 const Float_t rightMargin=0.07,
		 const Float_t topMargin=0.12
		 ) {
  canv->SetTopMargin(topMargin);
  canv->SetBottomMargin(bottomMargin);
  canv->SetLeftMargin(leftMargin);
  canv->SetRightMargin(rightMargin);
}

void rebin2(TH1 *h, Int_t ngx, Int_t ngy)
{
  //Rebin 2-d histogram h, grouping ngx bins together along X
  //and ngy bins together along Y
  //NB: this macro ignores histogram errors if defined

  //make a clone of h
  TH1 *hold = (TH1*)h->Clone();
  hold->SetDirectory(0);

  Int_t  nbinsx = hold->GetXaxis()->GetNbins();
  Int_t  nbinsy = hold->GetYaxis()->GetNbins();
  Float_t xmin  = hold->GetXaxis()->GetXmin();
  Float_t xmax  = hold->GetXaxis()->GetXmax();
  Float_t ymin  = hold->GetYaxis()->GetXmin();
  Float_t ymax  = hold->GetYaxis()->GetXmax();
  Int_t nx = nbinsx/ngx;
  Int_t ny = nbinsy/ngy;
  h->SetBins (nx,xmin,xmax,ny,ymin,ymax);

  //loop on all bins to reset contents and errors
  Double_t cu;
  Float_t bx,by;
  Int_t ix,iy,ibin,bin,binx,biny;
  for (biny=1;biny<=nbinsy;biny++) {
    for (binx=1;binx<=nbinsx;binx++) {
      ibin = h->GetBin(binx,biny);
      h->SetBinContent(ibin,0);
    }
  }
  //loop on all bins and refill
  for (biny=1;biny<=nbinsy;biny++) {
    by  = hold->GetYaxis()->GetBinCenter(biny);
    iy  = h->GetYaxis()->FindBin(by);
    for (binx=1;binx<=nbinsx;binx++) {
      bx = hold->GetXaxis()->GetBinCenter(binx);
      ix  = h->GetXaxis()->FindBin(bx);
      bin = hold->GetBin(binx,biny);
      ibin= h->GetBin(ix,iy);
      cu  = hold->GetBinContent(bin);
      h->AddBinContent(ibin,cu);
    }
  }
  delete hold;
}

TH1D* xShiftHist(TH1* h, double shiftVal){
  ////////////////////////////////////////////////////////////////
  // X-shift
  TH1D* h_org = (TH1D*) h->Clone(Form("%s_org",h->GetName()));
  TH1D* h_xshifted = (TH1D*) h->Clone(Form("%s_xshifted",h->GetName()));
  int TOTNBINS = h_org->GetNbinsX();
  for(int ibin=1; ibin<TOTNBINS; ++ibin){
    double xVal = h_org->GetBinCenter(ibin);
    int shiftBin = h_org->FindBin(xVal-shiftVal);
    if(shiftBin<1) continue;
    double tempVal = h_org->GetBinContent(shiftBin);
    double tempErr = h_org->GetBinError(shiftBin);
    h_xshifted->SetBinContent(ibin,tempVal);
    h_xshifted->SetBinError(ibin,tempErr);
    if(ibin < 100) cout << "ibin = " << ibin << ", xVal = " << xVal << ", shiftBin = " << shiftBin << ", h_org_binVal = " << h_org->GetBinContent(ibin) << ", h_xshifted_val = " << tempVal << endl;
  }
  return h_xshifted;
}

double chi2(TH1* h1, TH1* h2, double rMin=-1, double rMax=-1){
  double c = 0;
  int ndof = 0;
  double tempBinMin = 1;
  double tempBinMax = h1->GetNbinsX()+1;

  if(rMin!=-1) tempBinMin = h1->GetXaxis()->FindBin(rMin);
  if(rMax!=-1) tempBinMax = h1->GetXaxis()->FindBin(rMax);
  for(int i = tempBinMin; i < tempBinMax; ++i){
    double y1 = h1->GetBinContent(i);
    double y2 = h2->GetBinContent(i);
    double e1 = h1->GetBinError(i);
    double e2 = h2->GetBinError(i);

    double dy = y1-y2;
    double de2 = e1*e1+e2*e2;
    if(de2 > 0) {
      c += dy*dy/de2;
      ndof += 1;
    }
  }

  //return c;
  return c/ndof;
}

void saveHistogramsToPicture(TH1* h, const char* fileType="pdf", const char* caption="", const char* directoryToBeSavedIn="figures", int styleIndex=0, int rebin =1){
  TCanvas* c1=new TCanvas();
  if(rebin!=1)
    {
      h->Rebin(rebin);
    }

  if(styleIndex==1)
    {
      h->Draw("E");
    }
  else
    {
      h->Draw();
      if(h->InheritsFrom("TH2"))
        {
	  h->Draw("COLZ TEXT");    // default plot style for TH2 histograms
        }
    }

  if(strcmp(directoryToBeSavedIn, "") == 0)   // save in the current directory if no directory is specified
    {
      c1->SaveAs(Form("%s_%s.%s" ,h->GetName(),caption, fileType));  // name of the file is the name of the histogram
    }
  else
    {
      c1->SaveAs(Form("%s/%s_%s.%s", directoryToBeSavedIn ,h->GetName(),caption, fileType));
    }
  c1->Close();
}

void MergeTrees(TChain& tree, std::ifstream& filelist)
{
  if (not (filelist.is_open())) {
    std::cerr << "failed in open." << std::endl;
    exit(EXIT_FAILURE);
  }

  char line[BUFSIZ];
  int counter = 0; 

  while (not (filelist.eof())) {
    filelist.getline(line, sizeof(line));

    if (strcmp(line, "") == 0) continue;

    TFile a(line);
    if (!a.IsOpen()) {
      std::cerr << "failed in open :" << line << std::endl;
      exit(EXIT_FAILURE);
    }
    if (gROOT->FindObject(tree.GetName()) == 0) {
      std::cerr << tree.GetName() << " is not found in :" << line << std::endl;
      exit(EXIT_FAILURE);
    }

    counter ++ ;
    tree.AddFile(line);
  }

  std::cout << counter << " files are loaded." << std::endl;
  return;
}

void yjStyleRoot(){
  gStyle -> SetOptStat(0);
  TH1::SetDefaultSumw2();
}
void SetyjPadStyle(){
  gStyle->SetPaperSize(20,26);
  gStyle->SetPadTopMargin(0.05);
  gStyle->SetPadRightMargin(0.05);
  gStyle->SetPadBottomMargin(0.17);
  gStyle->SetPadLeftMargin(0.17);
}
void SetHistTitleStyle(double titlesize=0.05, double labelsize = 0.04){
  gStyle->SetTextFont(42); 
  gStyle->SetTextSize(0.04); 
  gStyle->SetLabelFont( 42, "X" ); 
  gStyle->SetLabelFont( 42, "Y" ); 
  gStyle->SetLabelFont( 42, "Z" ); 
  gStyle->SetTitleSize( titlesize, "X" ); 
  gStyle->SetTitleSize( titlesize, "Y" ); 
  gStyle->SetTitleSize( titlesize, "Z" ); 
  gStyle->SetLabelSize( labelsize, "X" ); 
  gStyle->SetLabelSize( labelsize, "Y" ); 
  gStyle->SetLabelSize( labelsize, "Z" ); 
}
//void SetHistTitleStyle(double titlesize=0.06, double titleoffset=0.04, double labelsize = 0.05, double labeloffset=0.01){
//    gStyle->SetTitleSize( titlesize, "X" ); gStyle->SetTitleOffset(titleoffset, "X");
//   gStyle->SetTitleSize( titlesize, "Y" ); gStyle->SetTitleOffset(titleoffset, "Y");
//  gStyle->SetLabelSize( labelsize, "X" ); gStyle->SetLabelOffset(labeloffset, "X");
// gStyle->SetLabelSize( labelsize, "Y" ); gStyle->SetLabelOffset(labeloffset, "Y");
//}
void SetHistTitleStyle(TH1* h, double titlesize=0.055, double titleoffset=0.1, double labelsize=0.05,double labeloffset=0.01){
  h->SetTitleSize( titlesize, "X" ); h->SetTitleOffset(titleoffset, "X");
  h->SetTitleSize( titlesize, "Y" ); h->SetTitleOffset(titleoffset, "Y");
  h->SetLabelSize( labelsize, "X" ); h->SetLabelOffset(labeloffset, "X");
  h->SetLabelSize( labelsize, "Y" ); h->SetLabelOffset(labeloffset, "Y");
  h->GetXaxis()->CenterTitle();
  h->GetYaxis()->CenterTitle();
}
void thisPadStyle(){
  gPad->SetLeftMargin(0.17);
  gPad->SetRightMargin(0.08);
  gPad->SetBottomMargin(0.15);
  gPad->SetTopMargin(0.05);
}
void SetPadStyle(){
  gStyle->SetPadLeftMargin(0.14);
  gStyle->SetPadRightMargin(0.08);
  gStyle->SetPadBottomMargin(0.15);
  gStyle->SetPadTopMargin(0.10);
}
void legStyle( TLegend *a)
{
  a->SetBorderSize(0);
  a->SetFillStyle(0);
  //  a->SetHeader(head);
  //  a->SetTextFont(62);
  //  a->SetTextSize(17);
  //  a->SetLineColor(1);
  //  a->SetLineStyle(1);
  //  a->SetLineWidth(1);
  //  a->SetFillColor(0);

}
void graphStyle(TGraph *g1=0, Int_t lstyle=1, Int_t lcolor=1, Int_t mstyle=20, Int_t mcolor=1, Int_t lwidth=1, Double_t msize=1.0)
{
  g1->SetLineStyle(lstyle);
  g1->SetLineColor(lcolor);
  g1->SetLineWidth(lwidth);
  g1->SetMarkerStyle(mstyle);
  g1->SetMarkerColor(mcolor);
  g1->SetMarkerSize(msize);
}

void hLineStyle(TH1 *h1=0, Int_t lstyle=1, Int_t lcolor=1, Int_t lwidth=1, Int_t lfst=0, Int_t lfcolor=0)
{
  h1->SetLineStyle(lstyle);
  h1->SetLineColor(lcolor);
  h1->SetLineWidth(lwidth);
  h1->SetFillStyle(lfst);
  h1->SetFillColor(lfcolor);
}

void hMarkerStyle(TH1 *h1=0, Int_t mstyle=20, Int_t mcolor=1, Double_t msize=1.0)
{
  h1->SetMarkerStyle(mstyle);
  h1->SetMarkerColor(mcolor);
  h1->SetMarkerSize(msize);
}
void drawText(const char *text, float xp, float yp, bool isRightAlign=0, int textColor=kBlack, double textSize=0.04, int textFont = 42){
  TLatex *tex = new TLatex(xp,yp,text);
  tex->SetTextFont(textFont);
  //   if(bold)tex->SetTextFont(43);
  tex->SetTextSize(textSize);
  tex->SetTextColor(textColor);
  tex->SetLineWidth(1);
  tex->SetNDC();
  if(isRightAlign) tex->SetTextAlign(31);
  tex->Draw();
}

//void drawText(const char *text, float xp, float yp, int textColor=kBlack, int textSize=14){
//TLatex *tex = new TLatex(xp,yp,text);
//tex->SetTextFont(43);
////   if(bold)tex->SetTextFont(43);
//tex->SetTextSize(textSize);
//tex->SetTextColor(textColor);
//tex->SetLineWidth(1);
//tex->SetNDC();
//tex->Draw();
//}
void jumSun(Double_t x1=0,Double_t y1=0,Double_t x2=1,Double_t y2=1,Int_t color=1, Double_t width=1)
{
  TLine* t1 = new TLine(x1,y1,x2,y2);
  t1->SetLineWidth(width);
  t1->SetLineStyle(7);
  t1->SetLineColor(color);
  t1->Draw();
}

void onSun(Double_t x1=0,Double_t y1=0,Double_t x2=1,Double_t y2=1,Int_t color=1, Double_t width=1)
{
  TLine* t1 = new TLine(x1,y1,x2,y2);
  t1->SetLineWidth(width);
  t1->SetLineStyle(1);
  t1->SetLineColor(color);
  t1->Draw();
}
double findCross(TH1* h1, TH1* h2, double& frac, double& effi, double& fracErr, double& effiErr){
  Int_t nBins = h1->GetNbinsX();
  double crossVal =0;
  int binAt0 = h1->FindBin(0);
  for(Int_t ix=binAt0; ix<=nBins ;ix++){
    float yy1 = h1->GetBinContent(ix);
    float yy2 = h2->GetBinContent(ix);
    if(yy2>yy1) {
      crossVal= h1->GetBinLowEdge(ix);
      break;
    }
  }
  int crossBin = h1->FindBin(crossVal);
  frac = 1 - (h2->Integral(1,crossBin) / h1->Integral(1,crossBin) );
  effi = ( h1->Integral(1,crossBin) / h1->Integral() );
  fracErr = frac * TMath::Sqrt( (1./h2->Integral(1,crossVal)) + (1./h1->Integral(1,crossVal)) );
  effiErr = ( TMath::Sqrt(h1->Integral(1,crossVal)) / h1->Integral() ) * TMath::Sqrt(1 - (h1->Integral(1,crossVal)/h1->Integral()) );

  return crossVal;
}

void ratioPanelCanvas(TCanvas*& canv,
		      const Float_t divRatio=0.4,
		      const Float_t leftOffset=0.,
		      const Float_t bottomOffset=0.,
		      const Float_t leftMargin=0.17,
		      const Float_t bottomMargin=0.3,
		      const Float_t edge=0.05) {
  if (canv==0) {
    //Error("makeMultiPanelCanvas","Got null canvas.");
    return;
  }
  canv->Clear();


  TPad* pad1 = new TPad("pad1","",0.0,divRatio,1.0,1.0);
  canv->cd();
  pad1->SetLeftMargin(leftMargin);
  pad1->SetRightMargin(edge);
  pad1->SetTopMargin(edge);
  pad1->SetBottomMargin(edge);
  pad1->Draw();
  pad1->cd();
  pad1->SetNumber(1);

  TPad* pad2 = new TPad("pad2","",0.0,0.0,1.0,divRatio);
  canv->cd();
  pad2->SetLeftMargin(leftMargin);
  pad2->SetRightMargin(edge);
  pad2->SetTopMargin(edge);
  pad2->SetBottomMargin(bottomMargin);
  pad2->Draw();
  pad2->cd();
  pad2->SetNumber(2);
}

void makeMultiPanelCanvas(TCanvas*& canv, const Int_t columns,
			  const Int_t rows, const Float_t leftOffset=0.,
			  const Float_t bottomOffset=0.,
			  const Float_t leftMargin=0.2,
			  const Float_t bottomMargin=0.2,
			  const Float_t edge=0.05) {
  if (canv==0) {
    //Error("makeMultiPanelCanvas","Got null canvas.");
    return;
  }
  canv->Clear();

  TPad* pad[columns][rows];

  Float_t Xlow[columns];
  Float_t Xup[columns];
  Float_t Ylow[rows];
  Float_t Yup[rows];
  Float_t PadWidth =
    (1.0-leftOffset)/((1.0/(1.0-leftMargin)) +
		      (1.0/(1.0-edge))+(Float_t)columns-2.0);
  Float_t PadHeight =
    (1.0-bottomOffset)/((1.0/(1.0-bottomMargin)) +
			(1.0/(1.0-edge))+(Float_t)rows-2.0);
  Xlow[0] = leftOffset;
  Xup[0] = leftOffset + PadWidth/(1.0-leftMargin);
  Xup[columns-1] = 1;
  Xlow[columns-1] = 1.0-PadWidth/(1.0-edge);

  Yup[0] = 1;
  Ylow[0] = 1.0-PadHeight/(1.0-edge);
  Ylow[rows-1] = bottomOffset;
  Yup[rows-1] = bottomOffset + PadHeight/(1.0-bottomMargin);

  for(Int_t i=1;i<columns-1;i++) {
    Xlow[i] = Xup[0] + (i-1)*PadWidth;
    Xup[i] = Xup[0] + (i)*PadWidth;
  }
  Int_t ct = 0;
  for(Int_t i=rows-2;i>0;i--) {
    Ylow[i] = Yup[rows-1] + ct*PadHeight;
    Yup[i] = Yup[rows-1] + (ct+1)*PadHeight;
    ct++;
  }
  TString padName;
  for(Int_t i=0;i<columns;i++) {
    for(Int_t j=0;j<rows;j++) {
      canv->cd();
      padName = Form("p_%d_%d",i,j);
      pad[i][j] = new TPad(padName.Data(),padName.Data(),
			   Xlow[i],Ylow[j],Xup[i],Yup[j]);
      if(i==0) pad[i][j]->SetLeftMargin(leftMargin);
      else pad[i][j]->SetLeftMargin(0);

      if(i==(columns-1)) pad[i][j]->SetRightMargin(edge);
      else pad[i][j]->SetRightMargin(0);

      if(j==0) pad[i][j]->SetTopMargin(edge);
      else pad[i][j]->SetTopMargin(0);

      if(j==(rows-1)) pad[i][j]->SetBottomMargin(bottomMargin);
      else pad[i][j]->SetBottomMargin(0);

      pad[i][j]->Draw();
      pad[i][j]->cd();
      pad[i][j]->SetNumber(columns*j+i+1);
    }
  }
}

Double_t getDPHI( Double_t phi1, Double_t phi2) {
  Double_t dphi = phi1 - phi2;

  if ( dphi > 3.141592653589 )
    dphi = dphi - 2. * 3.141592653589;
  if ( dphi <= -3.141592653589 )
    dphi = dphi + 2. * 3.141592653589;

  if ( TMath::Abs(dphi) > 3.141592653589 ) {
    std::cout << " commonUtility::getDPHI error!!! dphi is bigger than 3.141592653589 " << std::endl;
  }

  return dphi;
}

Double_t getDETA(Double_t eta1, Double_t eta2){
  return eta1 - eta2;
}

Double_t getDR( Double_t eta1, Double_t phi1, Double_t eta2, Double_t phi2){
  Double_t theDphi = getDPHI( phi1, phi2);
  Double_t theDeta = eta1 - eta2;
  return TMath::Sqrt ( theDphi*theDphi + theDeta*theDeta);
}

Double_t cleverRange(TH1* h,Float_t fac=1.2, Float_t minY=1.e-3)
{
  Float_t maxY =  fac * h->GetBinContent(h->GetMaximumBin());
  //   cout <<" range will be set as " << minY << " ~ " << maxY << endl;
  h->SetAxisRange(minY,maxY,"Y");
  return maxY;
}


Double_t getCleverRange(TH1* h)
{
  Double_t maxY = -1000000;
  for ( Int_t ibin = 1 ; ibin <= h->GetNbinsX() ; ibin++) {
    if (maxY < h->GetBinContent(ibin) )
      maxY = h->GetBinContent(ibin);
  }
  return maxY;
}

Double_t cleverRange(TH1* h,TH1* h2, Float_t fac=1.2, Float_t minY=1.e-3)
{
  Float_t maxY1 =  fac * h->GetBinContent(h->GetMaximumBin());
  Float_t maxY2 =  fac * h2->GetBinContent(h2->GetMaximumBin());

  //   cout <<" range will be set as " << minY << " ~ " << maxY << endl;
  h->SetAxisRange(minY,max(maxY1,maxY2),"Y");
  h2->SetAxisRange(minY,max(maxY1,maxY2),"Y");
  return max(maxY1,maxY2);
}

TF1* cleverGaus(TH1* h, const char* title="h", Float_t c = 2.5, bool quietMode=true)
{
  if ( h->GetEntries() == 0 )
    {
      TF1 *fit0 = new TF1(title,"gaus",-1,1);
      fit0->SetParameters(0,0,0);
      return fit0;
    }

  Int_t peakBin  = h->GetMaximumBin();
  Double_t peak =  h->GetBinCenter(peakBin);
  Double_t sigma = h->GetRMS();

  TF1 *fit1 = new TF1(title,"gaus",peak-c*sigma,peak+c*sigma);
  fit1->SetParameter(1, 1.0);
  fit1->SetParameter(1, 0.0005);
  if (quietMode) h->Fit(fit1,"LL M O Q R");
  else    h->Fit(fit1,"LL M O Q R");
  return fit1;
}

void SetHistColor(TH1* h, Int_t color=1)
{
  h->SetMarkerColor(color);
  h->SetLineColor(color);
} 

float mean(float data[], int n)
{
  float mean=0.0;
  int i;
  for(i=0; i<n;++i)
    {
      mean+=data[i];
    }
  mean=mean/n;
  return mean;
}

float standard_deviation(float data[], int n)
{
  float mean=0.0, sum_deviation=0.0;
  int i;
  for(i=0; i<n;++i)
    {
      mean+=data[i];
    }
  mean=mean/n;
  for(i=0; i<n;++i)
    sum_deviation+=(data[i]-mean)*(data[i]-mean);
  return sqrt(sum_deviation/n);
}

void normHist(TH1* h=0, TH1* hNominal=0, double cut_i=700, double cut_f=900){
  int cutBinFrom = h->FindBin(cut_i);
  int cutBinTo = h->FindBin(cut_f);
  h->Scale(hNominal->Integral(cutBinFrom,cutBinTo)/h->Integral(cutBinFrom,cutBinTo));
}


#endif
