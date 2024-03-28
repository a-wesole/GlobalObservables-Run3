#include "TFile.h"
#include "TTree.h"
#include "TNtuple.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TString.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>

#include "DataFormats/HeavyIonEvent/interface/CentralityBins.h"

using namespace std;


void makeMCCentralityTable(
			   //const TString input_file = "forest_HYD_RECODEBUG.txt",
			   const char * tag = "CentralityTable_HFtowers200_HydjetDrum5F_v1302x04_official_MC2023",
			   const char* date = "Mar28",
			   const string label = "hiHF")
{
  
  TH1::SetDefaultSumw2();

  int nbins = 200;
  //TString inFileName = "/eos/cms/store/group/phys_heavyions/nsaha/GO2023/2023PbPbRun3/forest_Run3_HYD_RECODEBUG_09102023/GENSIM_HydjetnoPU_CMSSW_13_2_0_pre3_28Jun2023/HiForest_Run3_HYD_RECODEBUG_09102023/231009_185152/0000/Forest_Run3_HYD_RECODEBUG_merged.root";

  TString inFileName = "/eos/cms/store/group/phys_heavyions/nsaha/GO2023/2023PbPbRun3/forest_Run3_HYD_official_23032024/MinBias_Drum5F_5p36TeV_hydjet/HiForest_Run3_HYD_official_23032024/240323_071705/0000/HYD_official_GT132X_mcRun3_2023_realistic_HI_v9_out_combined.root";
  TString outFileName = Form("CentralityTable_%s200_HydjetDrum5F_official_MC2023_%s.root",label.c_str(), date);
  ofstream txtfile(Form("output_HydjetDrum5F_hiHF_official_MC2023_%s.txt", date));

  TFile *inFile = TFile::Open(inFileName);

  TTree *t = (TTree*)inFile->Get("hiEvtAnalyzer/HiTree");
  TTree *t_skim = (TTree*)inFile->Get("skimanalysis/HltTree");
  t->AddFriend(t_skim);


  /*
  // Process data                                                                                                               
  TString Str;
  ifstream fpr(Form("%s",input_file.Data()), ios::in);
  if(!fpr.is_open()){
    cout << "List of input files not found!" << endl;
    return;
  }

  std::vector<TString> file_name_vector;
  string file_chain;
  while(getline(fpr, file_chain))
    {
      file_name_vector.push_back(file_chain);
    }

  TChain *t = new TChain("hiEvtAnalyzer/HiTree");
  TChain *tskimanalysis = new TChain("skimanalysis/HltTree");
  

  for (std::vector<TString>::iterator listIterator = file_name_vector.begin(); listIterator != file_name_vector.end(); listIterator++)
    {
      TFile *file = TFile::Open(*listIterator);
      cout << "Adding file:--- " << *listIterator << "--- to the chains" << endl;
      t->Add(*listIterator);
      tskimanalysis->Add(*listIterator);
  
    }

  t->AddFriend(tskimanalysis);
  */
  
    //TNtuple * nt = new TNtuple("nt","","value:bin:b:npart:ncoll:nhard");
  
  const int runNum = 1;
  CentralityBins * bins = new CentralityBins(Form("run%d",runNum), tag, nbins);
  bins->table_.reserve(nbins);
  

  txtfile << "Input tree: " << inFileName << endl;
  
  double binboundaries[nbins+1];
  vector<float> values;

  float hiHF;
  int hiBin, pprimaryVertexFilter, pclusterCompatibilityFilter, pphfCoincFilter2Th4;
  
  t->SetBranchAddress("hiHF", &hiHF);
  t->SetBranchAddress("hiBin", &hiBin);
  t->SetBranchAddress("pprimaryVertexFilter",	&pprimaryVertexFilter);
  t->SetBranchAddress("pclusterCompatibilityFilter", &pclusterCompatibilityFilter);
  t->SetBranchAddress("pphfCoincFilter2Th4", &pphfCoincFilter2Th4);

  t->SetBranchStatus("*", 0);
  for (const auto& p : {"hiHF", "pprimaryVertexFilter", "pclusterCompatibilityFilter", "pphfCoincFilter2Th4"})
    t->SetBranchStatus(p, 1);

  
  /*bool binB = label.compare("b") == 0;
  bool binNpart = label.compare("Npart") == 0;
  bool binNcoll = label.compare("Ncoll") == 0;
  bool binNhard = label.compare("Nhard") == 0;*/
  //  bool binHF = label.compare("hiHF") == 1;
  /*bool binHFplus = label.compare("HFtowersPlus") == 0;
  bool binHFminus = label.compare("HFtowersMinus") == 0;
  bool binHFplusTrunc = label.compare("HFtowersPlusTrunc") == 0;
  bool binHFminusTrunc = label.compare("HFtowersMinusTrunc") == 0;
  bool binNpix = label.compare("PixelHits") == 0;
  bool binNpixTrks = label.compare("PixelTracks") == 0;
  bool binNtrks = label.compare("Tracks") == 0;*/

  const bool pass = (pprimaryVertexFilter>0 && pclusterCompatibilityFilter>0 && pphfCoincFilter2Th4>0);
  
  unsigned int Nevents = t->GetEntries();
  txtfile << "Number of events = " << Nevents << endl << endl;
  
  for(unsigned int iev = 0; iev < Nevents; iev++) {
    
    if(iev%100000 == 0) cout<<"Processing event: " << iev << " / " << Nevents << endl;
    t->GetEntry(iev);


      float parameter = -1;
      if (pprimaryVertexFilter>0 && pclusterCompatibilityFilter>0 && pphfCoincFilter2Th4>0){
	parameter = hiHF;
	values.push_back(parameter);
      }
  }
  
  sort(values.begin(),values.end());
  
  txtfile << "-------------------------------------" << endl;
  txtfile << label.data() << " based cuts are: " << endl;
  txtfile << "(";
  
  int size = values.size();
  binboundaries[nbins] = values[size-1];
  
  for(int i = 0; i < nbins; i++) {
    int entry = (int)(i*(size/((float)(nbins))));
    if(entry < 0 || i == 0) binboundaries[i] = 0;
    else binboundaries[i] = values[entry];
  }
  
  for(int i = 0; i < nbins; i++) {
    if(binboundaries[i] < 0) binboundaries[i] = 0;
    txtfile << binboundaries[i] << ", ";
  }
  txtfile << binboundaries[nbins] << ")" << endl << "-------------------------------------" << endl;
  
  //dir->cd();
  
  /*for(unsigned int iev = 0; iev < Nevents; iev++) {
    
    if( iev % 100000 == 0 ) cout<<"Processing event : " << iev << endl;
    t->GetEntry(iev);

    //std::cout<<"*****************"<<std::endl;

    float parameter = -1;    
    if(pass){
      parameter = hiHF;
    }
    }*/
  txtfile<<"-------------------------------------"<<endl;
  txtfile<<"# Bin  BinEdge"<<endl;
  for(int i = 0; i < nbins; i++){
    
    int ii = nbins-i;
    bins->table_[i].bin_edge = binboundaries[ii-1];
    
    txtfile << i << " " << binboundaries[ii-1] << " " <<endl;
  }
  txtfile<<"-------------------------------------"<<endl;

  TFile * outFile = new TFile(outFileName,"recreate");
  TDirectory* dir = outFile->mkdir(tag);
  dir->cd();
  outFile->cd();
  dir->cd();
  bins->Write();
  //nt->Write();  
  bins->Delete();
  outFile->Write();
  txtfile.close();
  
}
