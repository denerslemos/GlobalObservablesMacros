#include "TFile.h"
#include "TTree.h"
#include <iostream>
#include "centralityCalibration.h"

using namespace std;

void Draw_hiBin(const TString input_file = "run374354_PhysicsHIPhysicsRawPrime0.txt",
		const char* HLT = "HLT_HIMinimumBiasHF1AND_v2",

		//const TString input_file = "run374719_HIPhysicsRawPrime0.txt",
		//const char* HLT = "HLT_HIMinimumBiasHF1ANDZDC2nOR_v3",
		
		const char* message = "Run374354"
		){

  

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
  TChain *thltanalysis = new TChain("hltanalysis/HltTree");
  
  for (std::vector<TString>::iterator listIterator = file_name_vector.begin(); listIterator != file_name_vector.end(); listIterator++)
    {
      TFile *file = TFile::Open(*listIterator);
      cout << "Adding file:--- " << *listIterator << "--- to the chains" << endl;
      t->Add(*listIterator);
      tskimanalysis->Add(*listIterator);
      thltanalysis->Add(*listIterator);
    }
  
  t->AddFriend(tskimanalysis);
  t->AddFriend(thltanalysis);
  

  float hiHF; 
  int hiBin, hiBinfinal, pphfCoincFilter2Th4;
  int phfCoincFilter3, pprimaryVertexFilter, pclusterCompatibilityFilter;
  int HLT_mb, HLT_HIMinimumBiasHF1AND_v2;

  t->SetBranchStatus("*", 0);
  for (const auto& p : {"hiHF", "hiBin", HLT, "pprimaryVertexFilter", "pclusterCompatibilityFilter", "pphfCoincFilter2Th4"}) t->SetBranchStatus(p, 1);

  t->SetBranchAddress("hiHF", &hiHF);
  t->SetBranchAddress("hiBin", &hiBin);
  t->SetBranchAddress(HLT, &HLT_mb);
  t->SetBranchAddress("pprimaryVertexFilter", &pprimaryVertexFilter);
  t->SetBranchAddress("pclusterCompatibilityFilter", &pclusterCompatibilityFilter);
  t->SetBranchAddress("pphfCoincFilter2Th4", &pphfCoincFilter2Th4);
  //t->SetBranchAddress("numMinHFTower4", &numMinHFTower4);

  TFile* outf = new TFile(Form("hiBin_Th100_MC_NOMINAL_%s.root", message),"recreate");  
  int newbin;
  TTree* t_out = new TTree("anaCentrality","analysis level centrality");
  t_out->Branch("newBin",&newbin,"newBin/I");
  t_out->Branch("oldBin",&hiBinfinal,"oldBin/I");

  int N = t->GetEntries();
  for(int i = 0; i < N; ++i){
    t->GetEntry(i);
    
    if(i % 100000 == 0) cout<<"processing event : "<<i<<endl;

    if(HLT_mb != 1) continue;
    if(pprimaryVertexFilter != 1) continue;
    if(pclusterCompatibilityFilter != 1) continue;
    if(pphfCoincFilter2Th4 != 1)continue;
    hiBinfinal = hiBin;
    int calibrated_hiBin =  GOCalibration2023PbPb_getHiBinFromhiHF(hiHF,"Run374354",i);
    newbin = calibrated_hiBin;    
    t_out->Fill();

    }
  
  outf->cd();
  t_out->Write();
  outf->Close();

}

