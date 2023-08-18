#include "utils.h"


using namespace std;
static bool doUPC = 0;
TDatime* date = new TDatime();
static int a=0;



Double_t effFunction(Double_t *x, Double_t *par)
{
    Double_t xx = x[0];
    Double_t fitVal;
    if (xx <= par[0]) fitVal = 0.5*(1+TMath::Erf((xx-par[1])/(par[2])));
    else fitVal = 0.5*(1+TMath::Erf((xx-par[1])/(xx*par[2]/par[0])));
    //else fitVal = 0.5*(1+TMath::Erf((xx-par[1])/(xx*par[2]/par[0])));
    //if (xx <= par[0]) fitVal = 0.5*(1+TMath::Erf((xx-par[1])*xx/(xx*par[2])));
    //else fitVal = 0.5*(1+TMath::Erf((xx-par[3])/(xx*par[4])));
    //if (xx <= par[0]) fitVal = 0.5*(1+TMath::Erf((xx-par[1])*xx/(xx*par[2])));
    //else fitVal = 0.5*(1+TMath::Erf((xx-par[1])/(xx*par[2])));

   // if (xx <= par[0]) return 0.5*(1+TMath::Erf((xx-par[1])*xx/(xx*par[2])));
   // else return 0.5*(1+TMath::Erf((xx-par[1])/(xx*par[2])));
    return fitVal;
}


TH1* scale(string var, TTree* nt, double s = 1, TCut cut = "", bool doDiffBins=0, int nbins = 200, double xmin=0, double xmax=6000, double normRangeMin=1000, double normRangeMax=4000, bool doMCreweight=false){
   
  Float_t hiHF;
  nt->SetBranchAddress("hiHF",&hiHF);
  // set non-constant bin width 
  int nBinWidthChange = 3;
  double binWidth[] = {1.0,2.0,20.0};
  double range[] = {15,200,xmax};
  
  int nDiffBins=0;
  float tempBins[1000];
  tempBins[0]=0;
  for(int i=0;tempBins[i]<xmax;++i){
    if(tempBins[i]<range[0]) tempBins[i+1] = tempBins[i]+binWidth[0];
    else if(tempBins[i]>=range[0] && tempBins[i]<range[1]) tempBins[i+1] = tempBins[i]+binWidth[1];
    else if(tempBins[i]>=range[1] && tempBins[i]<range[2]) tempBins[i+1] = tempBins[i]+binWidth[2];
    else cout << "[ERROR] : " <<  i << "th out of "<<nDiffBins<<"! somethings wrong";
    nDiffBins++;
  }
  double bins[nDiffBins+1];
  for(int i=0;i<nDiffBins+1;++i){
    bins[i] = tempBins[i];
  }
  
  // define histogram
  TH1D* h_temp1=0;
  TH1D* h_temp2=0;
  TH1D* h;
  
  if(!doMCreweight) {
    if(doDiffBins) h = new TH1D(Form("h%d",a),Form(";%s [GeV];event fraction",var.data()),nDiffBins,bins);
    else h = new TH1D(Form("h%d",a),Form(";%s [GeV];event fraction",var.data()),nbins,xmin,xmax);
    // fill and scale histogram with the condition of event selection
    nt->Draw(Form("%f*%s>>%s",s,var.data(),h->GetName()),cut);
  } else {
    float rangeDiv = 40;
    if(doDiffBins) h = new TH1D(Form("h%d",a),Form(";%s [GeV];event fraction",var.data()),nDiffBins,bins);
    else h = new TH1D(Form("h%d",a),Form(";%s [GeV];event fraction",var.data()),nbins,xmin,xmax);
    /////////////////////////////////////////////////
    // Event Loop for weighting
    nt->Draw(">>eventlist",cut);
    TEventList *elist = (TEventList*)gDirectory->Get("eventlist");
    Int_t nentries = elist->GetN();
    
    for (Long64_t i=0; i<nentries; ++i) {
      nt->GetEntry(elist->GetEntry(i));
      if(hiHF<rangeDiv){
	float weight = 0;
	//PV_BS_HF4, no shift, 1.06 scale (hiHF-0.0)*1.06
	//weight = 2;
	weight = (0.94038/(1+TMath::Exp(-0.10831*(hiHF-11.30114)))+3.06281/(hiHF-11.30114));
	if(hiHF<12) cout << "hiHF = " << hiHF << ", weight = " << weight << endl;
	h->Fill(hiHF,weight);
      } else{
	h->Fill(hiHF);
      }
    }
  }
  
  a++;
  h->GetXaxis()->CenterTitle();
  h->GetYaxis()->CenterTitle();
  
  // normalize within a certain range
  double rangeIntegral = 0;
  rangeIntegral = h->Integral(h->GetXaxis()->FindBin(normRangeMin),h->GetXaxis()->FindBin(normRangeMax));
  if(rangeIntegral > 0)
    h->Scale(1./rangeIntegral);
  h->Scale(1.,"width");
  
  return h;
  
}


void fit(TTree* tref=0, TTree* t=0, TTree* tupc=0, string var = "hiHF", bool doDiffBins=0, int nbins = 2000, double xmin=0, double xmax=4000, double normRangeMin=200, double normRangeMax=3500, TCut dataCut = "pPAprimaryVertexFilter && pBeamScrapingFilter && phfCoincFilter4", TCut mcCut2="", string cutname = "PV_BS_HF4", string mc="epos", double scaleMin=0.8, bool doMCevtCut=1, bool doNSDselection=0, string dataCap="XeXe 5.44 TeV", bool doMCreweight=false){
  cout << "============================================================================="<< endl; 
  cout << "DATA = " << dataCap << ", MC = " << mc << endl;
  cout << "The varible is "<< var << " on the condition of " << cutname << endl;
  cout << "doMCevtCut = " << doMCevtCut << "( doMCevtCut=1 is for sanity check for the method)" << endl;
  cout << "============================================================================="<< endl; 
  if(cutname == "") cutname = (const char*)dataCut;
  gStyle->SetOptStat(0);
  TH1::SetDefaultSumw2();
  
  TH1* hupc, *href;
  if(doUPC){
    //WARNING : No trigger for UPC
    hupc = scale(var,tupc,1.,dataCut,doDiffBins,nbins,xmin,xmax,normRangeMin,normRangeMax,0);
    hupc->SetLineColor(4);
    hupc->SetName("hupc");
  }
  
  //////////////////////////////////////////////////
  // event selection for data
  href = scale(var,tref,1.,dataCut,doDiffBins,nbins,xmin,xmax,normRangeMin,normRangeMax,0);
  cout<<"DATA entries with "<<cutname<<" from tree : "<<tref->GetEntries(dataCut)<<endl;
  cout<<"DATA entries with "<<cutname<<" from hist : "<<href->Integral("width")<<endl;
  
  //////////////////////////////////////////////////
  // scale MC to match DATA
  static const int N = 40;
  double variation = 0.2;
  TH1* h[N];
  double chi2s[N];
  double s[N];
  int ibest = 0;
  TGraphErrors* g = new TGraphErrors(N);
  
  TCut eposCut = "1";
  if(mc=="EPOS" && doNSDselection) eposCut = "ProcessID==101 || ProcessID==105";
  else if(mc=="EPOS" && !doNSDselection) eposCut = "ProcessID!=102";
  
  TCut mcCut = "1";
  if(doMCevtCut){
    if(mc=="EPOS") mcCut = eposCut && dataCut;
    //else mcCut = dataCut;
    else mcCut = mcCut2;
  }
  
  for(int i = 0; i < N; ++i){  
    s[i] = scaleMin+(variation/N)*i;
    h[i] = scale(var,t,s[i],mcCut,doDiffBins,nbins,xmin,xmax,normRangeMin,normRangeMax,doMCreweight);
    h[i]->SetLineColor(2);
    chi2s[i] = chi2(h[i], href, normRangeMin);
    g->SetPoint(i,s[i],chi2s[i]);
    if(chi2s[i]< chi2s[ibest]) ibest = i;
  }
  for(int i = 0; i < N; ++i){  
    if(i!=ibest)
      delete h[i];    
  }
  
  double dataTotalArea = href->Integral("width");
  double mcTotalArea = h[ibest]->Integral("width");
  //double mcTotalArea = h[ibest]->Integral(0,h[ibest]->GetXaxis()->FindBin(xmax),"width");
  double mcPeripheralArea = h[ibest]->Integral(0,h[ibest]->GetXaxis()->FindBin(normRangeMin),"width");
  double dataNonPeripheralArea = href->Integral(href->GetXaxis()->FindBin(normRangeMin)+1,href->GetNbinsX(),"width");
  double eff = dataTotalArea/(mcPeripheralArea + dataNonPeripheralArea); 
  
  std::cout << "############################" << std::endl;
  std::cout << dataTotalArea << " / " << " ( " << mcPeripheralArea << " + " << dataNonPeripheralArea << " ) = " << eff << std::endl;
  //////////////////////////
  // TEST! AFTER TEST PLEASE COMMENT OUT 
  //double dataPeripheralArea = href->Integral(0,href->GetXaxis()->FindBin(normRangeMin),"width");
  //cout << "TOTAL DATA AREA = " << dataTotalArea << ", DATA AREA SUM = " << dataNonPeripheralArea+dataPeripheralArea << endl;
  //////////////////////////
  
  cout<<"The best scale is : "<<s[ibest]<<" i = "<<ibest<<endl;   
  if(ibest==0 || ibest==N-1) cout << ":: WARNING :: it might not be the best chisquare. Adjust scaleMin!" << endl;
  cout<<"Efficiency+Contamination is : "<<eff<<endl;
  
  TString cap = Form("%s",var.data());
  cap += Form("_%s",cutname.data());
  cap += Form("_normRangeMin%d_normRangeMax%d",(int)normRangeMin,(int)normRangeMax);
  cap += Form("_%s",mc.data());
  if(doNSDselection) cap += "_NSDselected";
  cap += Form("_MCcut%d",(int)doMCevtCut);
  if(doDiffBins) cap += "_diffBinWidth";
  else cap += Form("_nBins%d",(int)nbins);
  //cap += Form("_%s",date->GetDate());
  
  ofstream outputf(Form("txt/centEff_fitResults_%s_%s.txt",cutname.data(),mc.data()), ios::app);
  outputf <<eff*100 << "\t" << chi2s[ibest] << endl;
  //outputf << run << "\t" << setprecision(4) <<eff*100 << "\t" << chi2s[ibest] << endl;
  outputf.close();
  
  /////////////////////////////////////
  // Chi2
  TCanvas* c1 = new TCanvas("c1","c1",600,600);
    canvasStyle(c1);
    g->SetTitle(";Scale factor;Reduced #chi^2 (normalized by ndf)");
    g->GetXaxis()->CenterTitle();
    g->GetYaxis()->CenterTitle();
    g->SetMarkerStyle(20);
    g->Draw("Ap");
    c1->Print(Form("figures/figure_%s_%s.pdf","chi2",cap.Data()));

    /////////////////////////////////////
    // MC and DATA Distribution  
    TCanvas* c2 = new TCanvas("c2","c2",600,600);
    canvasStyle(c2);
    c2->SetLogy();
    double ymin = h[ibest]->GetMinimum();
    double ymax = h[ibest]->GetMaximum();
    h[ibest]->SetMaximum(10.*ymax);
    h[ibest]->DrawCopy("hist");
    href->SetMarkerStyle(20);
    href->SetMarkerSize(0.8);
    href->DrawCopy("same");
    jumSun(normRangeMin, ymin, normRangeMin, 10*ymax);
    jumSun(normRangeMax, ymin, normRangeMax, 10*ymax);

    TLegend * leg1 = new TLegend(0.4,0.5,0.90,0.85);
    legStyle(leg1);
    leg1->SetTextSize(0.028);
    leg1->AddEntry(href,Form("%s",dataCap.data()),"");
    leg1->AddEntry(href,"DATA","p");
    leg1->AddEntry(href,Form("%s",cutname.data()),"");
    leg1->AddEntry(href,Form("norm : %.f to %.f",normRangeMin,normRangeMax),"");
    leg1->AddEntry(h[ibest],Form("%s based fit",mc.data()),"l");
    leg1->AddEntry(h[ibest],Form("Eff+Contam : %.2f %s",eff*100,"%"),"");
    leg1->AddEntry(h[ibest],Form("%s scaling : %.3f",mc.data(),s[ibest]),"");
    leg1->AddEntry(h[ibest],Form("#chi^{2}/ndf : %.3f",chi2s[ibest]),"");
    leg1->Draw();

    c2->Print(Form("figures/figure_%s_%s.pdf","fit",cap.Data()));

    /////////////////////////////////////
    // MC and DATA difference  
    TCanvas* c3 = new TCanvas("c3","c3",600,600);
    canvasStyle(c3);
    c3->SetLogx(0);
    TH1D* hdiff = (TH1D*)href->Clone("hdiff");
    hdiff->SetTitle(Form(";%s [GeV];DATA - %s",var.data(),mc.data()));
    hdiff->Add(h[ibest],-1);
	hdiff->GetXaxis()->SetRangeUser(0,1000);
    hdiff->Draw();
    jumSun(0,0,xmax,0);

    double supc = 0;
    if(doUPC){
        supc = hdiff->Integral(1,5)/hupc->Integral(1,5);
        hupc->Scale(supc);
        hupc->Draw("hist same");
    }

    //TLegend * leg2 = new TLegend(0.4,0.2,0.95,0.4);
    TLegend * leg2 = new TLegend(0.4,0.7,0.95,0.85);
    legStyle(leg2);
    leg2->SetTextSize(0.025);
    leg2->AddEntry(hdiff,Form("DATA - %s",mc.data()),"p");
    leg2->AddEntry(href,Form("%s",cutname.data()),"");
    leg2->AddEntry(href,Form("normRange : %.0f to %.0f",normRangeMin, normRangeMax),"");
    if(doUPC){
        leg2->AddEntry(hupc,"Starlight","l");
        leg2->AddEntry(hupc,Form("Vertical scaling : %f",supc),"");
    }
    leg2->Draw();

    c3->Print(Form("figures/figure_%s_%s.pdf","diff",cap.Data()));

    /////////////////////////////////////
    // Draw Effciency Curve : Distribution in zoomed x
    TCanvas* c_evtSelEff = new TCanvas("c2_2","c2",500,750);
    canvasStyle(c_evtSelEff);
    ratioPanelCanvas(c_evtSelEff);
    c_evtSelEff->cd(1);
    int drawRangeMin = 0;
    int drawRangeMax = 50;
    h[ibest]->GetXaxis()->SetRangeUser(drawRangeMin,drawRangeMax);
    h[ibest]->GetYaxis()->SetRangeUser(0,ymax*1.15);
    h[ibest]->DrawCopy("hist");
    href->DrawCopy("same");

    /////////////////////////////////////
    // Draw Effciency Curve : Ratio
    c_evtSelEff->cd(2);
    TH1D* ratio = (TH1D*) href->Clone("ratio");
    ratio->Divide(href,h[ibest]);
    ratio->GetXaxis()->SetRangeUser(drawRangeMin,drawRangeMax);
    ratio->GetYaxis()->SetRangeUser(0,1.2);
    ratio->GetYaxis()->SetTitle(Form("Eff.+Contam.(=DATA/%s)",mc.data()));
    ratio->GetXaxis()->SetTitle("#Sigma E_{HF} (GeV)");
    ratio->Draw();
    SetHistTitleStyle(ratio,0.055,1);
    jumSun(0,1,drawRangeMax,1);

    /////////////////////////////////////
    // FITTING
    int functionN = 4;
    // f_temp = new TF1("func", "[0]/(TMath::Exp((-x+[1])/[2]+1)+[3])+[4]", 3.3, 20);
    TF1 *f_temp;
    if(functionN==0) f_temp = new TF1("func", " [0]*(1+TMath::Erf((x-[1])/(TMath::Sqrt(x*x)*[2])))", 0, xmax); //nominal
    else if(functionN==1) f_temp = new TF1("func", " [0]*(1+TMath::Erf((x-[1])*x/(TMath::Sqrt(x*x)*[2])))", 0, xmax);//well-matched at high efficiency region
    else if(functionN==2) f_temp = new TF1("func", " [0]*(1+TMath::Erf((x-[1])/(TMath::Sqrt(x)*[2])))", 0, xmax);
    else if(functionN==3) f_temp = new TF1("func", " [0]*(1+TMath::Erf((x-[1])*x/(TMath::Sqrt(x)*[2])))", 0, xmax);
    else if(functionN==4) f_temp = new TF1("func", effFunction, 0, xmax,3); 
    //TF1 *f_temp = new TF1("func", " [0]*(1+TMath::Erf((x-[1])/(TMath::Sqrt(x*x)*[2])))", 0, xmax);
    //TF1 *f_temp = new TF1("func", " [0]*(1+TMath::Erf((x-[1])/(TMath::Sqrt(x*x)*[2]))*([3]*x+[4]))", 0, xmax);
    //TF1 *f_temp = new TF1("func", " [0]+[1]*(1+TMath::Erf((x-[2])/(TMath::Sqrt(x)*[3])))", 0, xmax);
    //TF1 *f_temp = new TF1("func", "[0]+[1]*(1+TMath::Erf((x-[2])/TMath::Sqrt(x)))", 0, xmax);
    //TF1 *f_temp = new TF1("func", "[0]*TMath::ATan([1]*x) +[2]", 0, xmax);
    //TF1 *f_temp = new TF1("func", " (exp([0]-x)*TMath::Erf((x-[1])/[2]))+[3]", 0, xmax);
    //TF1 *f_temp = new TF1("func", " (exp(-[0]*x)*TMath::Erf((x-[1])/[2]))+[3]", 0, xmax);
    //TF1 *f_temp = new TF1("func", "([0]*TMath::Erf((x-[1])/[2]))+[3]+[4]*x*x", 0, xmax);
    //TF1 *f_temp = new TF1("func", "(TMath::Erf([0]+[1]*x+[2]*x*x))", 0, xmax);
    //TF1 *f_temp = new TF1("func", "([0]*TMath::Erf((x-[1]/[2])))+[3]", 0, xmax);
    //TF1 *f_temp = new TF1("func", "([0]*TMath::Erf((x-[1])/[2]))+[3]", 0, xmax);
    //TF1 *f_temp = new TF1("func", "[0]/(TMath::Exp((-x+[1])/[2]+1)+[3])+[4]", 0, xmax);
    //TF1 *f_erf = new TF1("erf(x)","[0]TMath::erf(x-[1])+[2]",0,xmax);
    //TF1 *f_erf = new TF1("erf(x)","ROOT::Math::erf(x+[1])+[2]",0,xmax);
    //TF1 *f_gaussian_cdf = new TF1("gaussian_cdf(x)","ROOT::Math::gaussian_cdf(x)",0,xmax);

    f_temp->SetNpx(1000);
    if(functionN==3) f_temp->SetParameters(0.50,1.34114e+01,7.5804);
    else if(functionN==4) f_temp->SetParameters(12,1.34114e+01,7.5804,13.3955,0.2563);
    else f_temp->SetParameters(0.5,1.34114e+01,7.45804e-01);
    if(functionN!=4) f_temp->FixParameter(0, 0.5);
    if(functionN==4) f_temp->SetParLimits(0, 11,25);

    ratio->Fit(f_temp);
    ratio->Fit(f_temp);
    f_temp->Draw("same");

    drawText(Form("Normalization range [%d,%d]",(int)normRangeMin,(int)normRangeMax),0.20,0.9);
    if(functionN==0) drawText(Form("0.5*(1+Erf( #frac{(x-%.4f)}{%.4f*x} ))",f_temp->GetParameter(1),f_temp->GetParameter(2)),0.40,0.3);
    else if(functionN==1) drawText(Form("0.5*(1+Erf( #frac{x*(x-%.4f)}{%.4f*x} ))",f_temp->GetParameter(1),f_temp->GetParameter(2)),0.40,0.3);
    else if(functionN==2)drawText(Form("0.5*(1+Erf( #frac{x-%.4f}{%.4f*#sqrt{x}} ))",f_temp->GetParameter(1),f_temp->GetParameter(2)),0.40,0.3);
    else if(functionN==3)drawText(Form("0.5*(1+Erf( #frac{x*(x-%.4f)}{%.4f*#sqrt{x}} ))",f_temp->GetParameter(1),f_temp->GetParameter(2)),0.40,0.3);
    else if(functionN==4) {
        double xpos = 0.6;
        double ypos = 0.7;
        drawText(Form("if x#leq%.2f",f_temp->GetParameter(0)),xpos,ypos);
        drawText(Form("0.5*(1+Erf( #frac{x-%.4f}{%.4f} ))",f_temp->GetParameter(1),f_temp->GetParameter(2)),xpos+0.02,ypos-0.06);
        drawText(Form("if x>%.2f",f_temp->GetParameter(0)),xpos, ypos-0.15);
        drawText(Form("0.5*(1+Erf( #frac{x-%.4f}{%.4f*x} ))",f_temp->GetParameter(1),f_temp->GetParameter(2)/f_temp->GetParameter(0)),xpos+0.02,ypos-0.21);
    }

    ////////////////////////////////
    // Cross check
    TH1D* hinv = (TH1D*) href->Clone("hinv");
    for(int i=1;i<href->GetNbinsX()+1;++i){
        double xVal = href->GetBinContent(i);
        double xErr = href->GetBinError(i);
        double x = href->GetBinCenter(i);
        double scale = 1./f_temp->Eval(x);
        if(xVal==0){
            hinv->SetBinContent(i,0);
            hinv->SetBinError(i,0);
        } else{
            hinv->SetBinContent(i,xVal*scale);
            hinv->SetBinError(i,xErr);
        }
    }
    c_evtSelEff->cd(1);
    hinv->SetMarkerColor(4);
    hinv->SetMarkerStyle(32);
    hinv->Draw("p same");
    double hinvTotalArea = hinv->Integral(0,hinv->GetXaxis()->FindBin(xmax),"width");
    double sanityCheck = mcTotalArea/hinvTotalArea;
    cout << "MC : DATA/Eff. : MC/(DATA/Eff.) = " << mcTotalArea << " : " << hinvTotalArea << " : " << sanityCheck << endl; 

    double effinv = (dataTotalArea)/(dataNonPeripheralArea + hinv->Integral(0,hinv->GetXaxis()->FindBin(normRangeMin),"width"));
    double effinv_2 = href->Integral("width")/hinv->Integral("width");
    cout << "eff  = " << eff << endl;
    cout << "eff inverted = " << effinv << endl;
    cout << "eff inverted 2 = " << effinv_2 << endl;

    TLegend * leg3 = new TLegend(0.50,0.50,0.95,0.95);
    legStyle(leg3);
    //leg3->SetTextSize(0.025);
    leg3->AddEntry(href,Form("%s",dataCap.data()),"");
    leg3->AddEntry(href,"DATA","p");
    leg3->AddEntry(href,Form("%s",cutname.data()),"");
    leg3->AddEntry(h[ibest],Form("%s based fit",mc.data()),"l");
    leg3->AddEntry(h[ibest],Form("Eff+Contam : %.2f %s",eff*100,"%"),"");
    //leg3->AddEntry(h[ibest],Form("%s scaling : %.3f",mc.data(),s[ibest]),"");
    //leg3->AddEntry(h[ibest],Form("#chi^{2}/ndf : %.3f",chi2s[ibest]),"");
    leg3->AddEntry(hinv,"Data/Eff.","pl");
    leg3->AddEntry(h[ibest],Form("MC/(Data/Eff.)) : %.2f %s",sanityCheck*100,"%"),"");
    leg3->Draw();
    //drawText(Form("Eff.+Cont.(from Data/Eff.) = %.2f %s",effinv*100,"%"),0.57,0.55-0.05);
    //drawText("Evt. Sel. Eff.(from Data/Eff.)",0.6,0.6-0.05);
    //drawText(Form("= %.2f %s",effinv*100,"%"),0.5+0.25,0.6-2*0.05);
    c_evtSelEff->Print(Form("figures/figure_%s_%s.pdf","effCurve",cap.Data()));
    
    /////////////////////////////////////
    // Save output root file with histograms 
    TString fout_name = Form("output/fitting_%s.root",cap.Data());
    TFile* fout = new TFile(fout_name,"recreate");
    fout->cd();
    h[ibest]->SetName("h[ibest]");
    href->SetName("href"); 
    h[ibest]->Write(); 
    href->Write(); 
    hinv->Write(); 
    f_temp->Write(); 
    g->Write();
    fout->Close();
    cout << "output file : output/fitting_" << cap << ".root has been created" << endl;
}


void FittingMethod(){

  TFile *infData, *infhydjet, *infepos, *infampt, *infUPC;
  TTree *tref, *thydjet, *tepos, *tampt, *tupc;

  //infData = TFile::Open("/eos/cms/store/group/phys_heavyions/nsaha/GO2023/HiForest_TestRun_data/HITestRaw/HiForestAODZB_TestRun_Data_ZB1_merged_small.root");

  infData = TFile::Open("/eos/cms/store/group/phys_heavyions/anstahll/GO2023/HIForest_HIMinimumBias_HIRun2022A_20230811_SKIM.root");

  tref = (TTree*)infData->Get("hiEvtAnalyzer/HiTree");
  tref->AddFriend("hltanalysis/HltTree");
  tref->AddFriend("skimanalysis/HltTree");

  //infhydjet = TFile::Open("/eos/cms/store/group/phys_heavyions/nsaha/GO2023/HiForest_HYDJET_official_modifiedGT/MinBias_Drum5F_5p36TeV_hydjet/HiForest_TestRun2022_MC_official_GT/230809_191838/0000/HiForestMiniAOD_OfficialMC_Hyd2022_merged.root");

  infhydjet = TFile::Open("/eos/cms/store/group/phys_heavyions/anstahll/GO2023/HiForestMiniAOD_HYDJET_June7_test_SKIM.root");

  thydjet = (TTree*)infhydjet->Get("hiEvtAnalyzer/HiTree");
  thydjet->AddFriend("skimanalysis/HltTree");
  thydjet->AddFriend("hltanalysis/HltTree");


  TCut dataTotCut = "run <=362320 && HLT_HIMinimumBias_v2 && pprimaryVertexFilter && pclusterCompatibilityFilter && numMinHFTower4>=2";

  TCut mcCut2 = "HLT_HIMinimumBias_v2 && pprimaryVertexFilter && pclusterCompatibilityFilter && numMinHFTower4>=2";


  static const int Nvar = 1;
  string vars[Nvar] = {"hiHF"};

  double scaleMin_hydjet[Nvar] = {0.70};

  for(int i = 0; i < Nvar; ++i){
    fit(tref,thydjet,tupc,vars[i].data(),1,100,0,10000,1000,4000, dataTotCut, mcCut2, "MBtrig_PV_CCF_HF2Th4","HYDJET",scaleMin_hydjet[i],1,0,"PbPb 5.36 TeV", false);
    fit(tref,thydjet,tupc,vars[i].data(),0,100,0,10000,1000,4000, dataTotCut, mcCut2, "MBtrig_PV_CCF_HF2Th4","HYDJET",scaleMin_hydjet[i],1,0,"PbPb 5.36 TeV", false);
  }

}




//void fit(TTree* tref=0, TTree* t=0, TTree* tupc=0, string var = "hiHF", bool doDiffBins=0, int nbins = 2000, double xmin=0, double xmax=4000, double normRangeMin=200, double normRangeMax=3500, TCut dataCut = "pPAprimaryVertexFilter && pBeamScrapingFilter && phfCoincFilter4", TCut mcCut2="", string cutname = "PV_BS_HF4", string mc="epos", double scaleMin=0.8, bool doMCevtCut=1, bool doNSDselection=0, string dataCap="XeXe 5.44 TeV", bool doMCreweight=false)
