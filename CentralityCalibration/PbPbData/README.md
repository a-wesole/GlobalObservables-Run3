You'll need to change the following accordingly!

* L29: TString input_file = "Data_forest.txt //data forest files

* L31: const char* HLT_trg = "HLT_HIMinimumBiasHF1ANDZDC1nOR_v4 //HLT trigger

* L36: const char* CoinFilter = "pphfCoincFilterPF3Th5", // CoincFilter

* L252: TFile inputMCfile ("path/MC_forest.root", "READ") //MC forest file

* //first run the *Nominal.C
* //notes are there 

-- you also need a python file 
'makeDMFromTFile.py' 
how to run it cmsenv
cmsRun that_python_code.py 
