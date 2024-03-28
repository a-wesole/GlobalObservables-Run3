#include "TFile.h"
#include "TTree.h"
#include <iostream>
using namespace std;

void makeCentrality(){

  TFile* infData = TFile::Open("/eos/cms/store/group/phys_heavyions/nsaha/GO2023/2023PbPbRun3/forest_Run3_HYD_official_23032024/MinBias_Drum5F_5p36TeV_hydjet/HiForest_Run3_HYD_official_23032024/240323_071705/0000/HYD_official_GT132X_mcRun3_2023_realistic_HI_v9_out_combined.root");

  float bounds[201] = {0, 15.0365, 16.3021, 17.3476, 18.3277, 19.3203, 20.2781, 21.2408, 22.2337, 23.2456, 24.3082, 25.4276, 26.567, 27.7402, 28.9504, 30.1923, 31.5016, 32.8665, 34.2718, 35.7318, 37.2509, 38.8256,40.4977, 42.2255, 44.0698, 46.0185, 48.0379, 50.0949, 52.2153, 54.4985, 56.8187, 59.2729, 61.8737, 64.6133, 67.4205, 70.357, 73.422, 76.5244, 79.7694, 83.1915, 86.6754, 90.3001, 94.0548, 97.9411, 102.009, 106.248, 110.552, 115.184, 119.935, 124.858, 129.798, 134.979, 140.409, 146.117, 151.898, 157.978, 164.199, 170.559, 177.1, 184.063, 191.166, 198.364, 205.969, 213.741, 221.796, 230.179, 238.81, 247.581, 256.624, 265.723, 275.253, 285.127, 295.458, 305.858, 316.574, 327.518, 338.791, 350.37, 362.262, 374.273, 386.983, 399.716, 412.856, 426.142, 440.284, 454.66, 469.073, 483.935, 499.153, 514.763, 530.465, 547.383, 564.34, 581.293, 599.188, 617.165, 635.016, 653.826, 672.91, 692.084, 711.944, 732.364, 753.297, 774.166, 795.464, 817.674, 840.086, 862.789, 886.251, 909.937, 934.059, 958.609, 984.044, 1009.93, 1036, 1062.89, 1089.63, 1116.86, 1144.52, 1172.52, 1201.75, 1231.18, 1260.42, 1291.67, 1322.75, 1354.25, 1386.73, 1419.17, 1451.9, 1485.42, 1519.31, 1555.19, 1590.89, 1626.27, 1662.77, 1699.2, 1737.53, 1776.21, 1815.6, 1855.92, 1897.57, 1939.22, 1981.14, 2024.42, 2067, 2111.15, 2155.79, 2201.62, 2247.9, 2294.02, 2342.46, 2391.8, 2441.55, 2491.13, 2542.18, 2592.97, 2645.65, 2698.57, 2753.28, 2808.11, 2864.78, 2922.48, 2981.14, 3039.09, 3098.47, 3159.66, 3221.18, 3283.42, 3348.79, 3413.59, 3479.99, 3546.98, 3617.37, 3688.87, 3760.49, 3834.05, 3908.83, 3986.92, 4063.08, 4141.67, 4220.78, 4303.65, 4387.79, 4472.19, 4559.3, 4647.69, 4738.31, 4831.4, 4928.41, 5024.6, 5123.44, 5224.7, 5328.55, 5437.2, 5543.26, 5657.01, 5774.24, 5901.36, 6049.93, 6245.47, 7387.56};



  float hiHF; 
  int hiBin, pprimaryVertexFilter, pclusterCompatibilityFilter, pphfCoincFilter2Th4;

  TTree* tref = (TTree*)infData->Get("hiEvtAnalyzer/HiTree");
  TTree *tskim = (TTree*)infData->Get("skimanalysis/HltTree");

  tref->AddFriend(tskim);
  
  tref->SetBranchAddress("hiHF", &hiHF);
  tref->SetBranchAddress("hiBin", &hiBin);
  tref->SetBranchAddress("pprimaryVertexFilter",   &pprimaryVertexFilter);
  tref->SetBranchAddress("pclusterCompatibilityFilter", &pclusterCompatibilityFilter);
  tref->SetBranchAddress("pphfCoincFilter2Th4", &pphfCoincFilter2Th4);



  TFile* outf = new TFile("compare_centralitybins_MCHYD_official_28Mar.root","recreate");

  int newbin;
  TTree* t = new TTree("anaCentrality","analysis level centrality");
  t->Branch("newBin",&newbin,"newBin/I");
  t->Branch("oldBin",&hiBin,"oldBin/I");

  int N = tref->GetEntries();
  for(int i = 0; i < N; ++i){
    tref->GetEntry(i);

    if(i % 100000 == 0) cout<<"processing event : "<<i<<endl;

    if(pprimaryVertexFilter != 1 && pclusterCompatibilityFilter != 1 && pphfCoincFilter2Th4 != 1) continue;

    newbin = 199;
    for(int b = 0; b < 200; ++b){
      if(hiHF >= bounds[199-b]){
	newbin = b;
	break;	    
      }
    }

    t->Fill();
  }

  t->Write();
  outf->Write();

}


