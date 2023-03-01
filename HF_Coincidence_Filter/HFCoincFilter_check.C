#import "call_libraries.h"


void HFCoincFilter_check(TString inputfile, TString outputfile){

	using namespace std;

	clock_t sec_start,sec_end;

	sec_start = clock();

	// Read the input file(s)
	fstream inputfile;
	inputfile.open(Form("%s",input_file.Data()), ios::in);
	if(!inputfile.is_open()){cout << "List of input files not founded!" << endl; return;}{cout << "List of input files founded! --> " << input_file.Data() << endl;}
	// Make a chain and a vector of file names
	std::vector<TString> file_name_vector;
	string file_chain;
	while(getline(inputfile, file_chain)){file_name_vector.push_back(file_chain.c_str());}
	inputfile.close();

	// Read centrality and filter tree's
	TChain *hea_tree = new TChain("hiEvtAnalyzer/HiTree");
	TChain *ski_tree = new TChain("skimanalysis/HltTree");

	for (std::vector<TString>::iterator listIterator = file_name_vector.begin(); listIterator != file_name_vector.end(); listIterator++){
		cout << "Adding file " << *listIterator << " to the chains" << endl;
		hea_tree->Add(*listIterator);
		ski_tree->Add(*listIterator);
	}

	hea_tree->AddFriend(ski_tree);	

	int hiBin; // event centrality
	float vertexz; // event z vertex
	int pprimaryVertexFilter, pclusterCompatibilityFilter; // event filters
	float hiHF_pfle; // PF leading energy
	int nCountsHFPlus_pf, nCountsHFMinus_pf; // number of pf candidates in HF+ and HF- sides
	int numMinHFTower2, numMinHFTower3, numMinHFTower4, numMinHFTower5; // minimum number of towers (towermaker) between HF+ and HF- for thresholds 2, 3, 4 and 5 GeV

	hea_tree->SetBranchStatus("*", 0); // disable all branches - this is important while reading big files
    hea_tree->SetBranchStatus("pprimaryVertexFilter", 1);
    hea_tree->SetBranchAddress("pprimaryVertexFilter", &pprimaryVertexFilter);    
    hea_tree->SetBranchStatus("pclusterCompatibilityFilter", 1);
    hea_tree->SetBranchAddress("pclusterCompatibilityFilter", &pclusterCompatibilityFilter);    
    hea_tree->SetBranchStatus("hiBin", 1);
    hea_tree->SetBranchAddress("hiBin", &hiBin);
    hea_tree->SetBranchStatus("vz", 1);
    hea_tree->SetBranchAddress("vz", &vertexz);
    hea_tree->SetBranchStatus("hiHF_pfle", 1);
    hea_tree->SetBranchStatus("nCountsHFPlus_pf", 1);
    hea_tree->SetBranchAddress("nCountsHFPlus_pf", &nCountsHFPlus_pf); 
    hea_tree->SetBranchStatus("nCountsHFMinus_pf", 1);
    hea_tree->SetBranchAddress("nCountsHFMinus_pf", &nCountsHFMinus_pf); 
    hea_tree->SetBranchStatus("numMinHFTower2", 1);
    hea_tree->SetBranchAddress("numMinHFTower2", &numMinHFTower2);  
    hea_tree->SetBranchStatus("numMinHFTower3", 1);
    hea_tree->SetBranchAddress("numMinHFTower3", &numMinHFTower3);  
    hea_tree->SetBranchStatus("numMinHFTower4", 1);
    hea_tree->SetBranchAddress("numMinHFTower4", &numMinHFTower4);  
    hea_tree->SetBranchStatus("numMinHFTower5", 1);
    hea_tree->SetBranchAddress("numMinHFTower5", &numMinHFTower5);  

	// make histograms needed (all hiBin)
	TH1D *hist_hibin_all = new TH1D("hist_hibin_all","hist_hibin_all",201,0.0,201.0);
	hist_hibin_all->Sumw2();
	TH1D *hist_hibin_PV = new TH1D("hist_hibin_PV","hist_hibin_PV",201,0.0,201.0);
	hist_hibin_PV->Sumw2();
	TH1D *hist_hibin_CC = new TH1D("hist_hibin_CC","hist_hibin_CC",201,0.0,201.0);
	hist_hibin_CC->Sumw2();
	TH1D *hist_hibin_PVCC = new TH1D("hist_hibin_PVCC","hist_hibin_PVCC",201,0.0,201.0);
	hist_hibin_PVCC->Sumw2();

	// from HF Tower
	TH1D *hist_hibin_HF1Th2 = new TH1D("hist_hibin_HF1Th2","hist_hibin_HF1Th2",201,0.0,201.0);
	hist_hibin_HF1Th2->Sumw2();
	TH1D *hist_hibin_HF2Th2 = new TH1D("hist_hibin_HF2Th2","hist_hibin_HF2Th2",201,0.0,201.0);
	hist_hibin_HF2Th2->Sumw2();
	TH1D *hist_hibin_HF3Th2 = new TH1D("hist_hibin_HF3Th2","hist_hibin_HF3Th2",201,0.0,201.0);
	hist_hibin_HF3Th2->Sumw2();
	TH1D *hist_hibin_HF4Th2 = new TH1D("hist_hibin_HF4Th2","hist_hibin_HF4Th2",201,0.0,201.0);
	hist_hibin_HF4Th2->Sumw2();

	TH1D *hist_hibin_HF1Th3 = new TH1D("hist_hibin_HF1Th3","hist_hibin_HF1Th3",201,0.0,201.0);
	hist_hibin_HF1Th3->Sumw2();
	TH1D *hist_hibin_HF2Th3 = new TH1D("hist_hibin_HF2Th3","hist_hibin_HF2Th3",201,0.0,201.0);
	hist_hibin_HF2Th3->Sumw2();
	TH1D *hist_hibin_HF3Th3 = new TH1D("hist_hibin_HF3Th3","hist_hibin_HF3Th3",201,0.0,201.0);
	hist_hibin_HF3Th3->Sumw2();
	TH1D *hist_hibin_HF4Th3 = new TH1D("hist_hibin_HF4Th3","hist_hibin_HF4Th3",201,0.0,201.0);
	hist_hibin_HF4Th3->Sumw2();

	TH1D *hist_hibin_HF1Th4 = new TH1D("hist_hibin_HF1Th4","hist_hibin_HF1Th4",201,0.0,201.0);
	hist_hibin_HF1Th4->Sumw2();
	TH1D *hist_hibin_HF2Th4 = new TH1D("hist_hibin_HF2Th4","hist_hibin_HF2Th4",201,0.0,201.0);
	hist_hibin_HF2Th4->Sumw2();
	TH1D *hist_hibin_HF3Th4 = new TH1D("hist_hibin_HF3Th4","hist_hibin_HF3Th4",201,0.0,201.0);
	hist_hibin_HF3Th4->Sumw2();
	TH1D *hist_hibin_HF4Th4 = new TH1D("hist_hibin_HF4Th4","hist_hibin_HF4Th4",201,0.0,201.0);
	hist_hibin_HF4Th4->Sumw2();

	TH1D *hist_hibin_HF1Th5 = new TH1D("hist_hibin_HF1Th5","hist_hibin_HF1Th5",201,0.0,201.0);
	hist_hibin_HF1Th5->Sumw2();
	TH1D *hist_hibin_HF2Th5 = new TH1D("hist_hibin_HF2Th5","hist_hibin_HF2Th5",201,0.0,201.0);
	hist_hibin_HF2Th5->Sumw2();
	TH1D *hist_hibin_HF3Th5 = new TH1D("hist_hibin_HF3Th5","hist_hibin_HF3Th5",201,0.0,201.0);
	hist_hibin_HF3Th5->Sumw2();
	TH1D *hist_hibin_HF4Th5 = new TH1D("hist_hibin_HF4Th5","hist_hibin_HF4Th5",201,0.0,201.0);
	hist_hibin_HF4Th5->Sumw2();

	// from PF Candidates

	TH1D *hist_hibin_PF1Th1 = new TH1D("hist_hibin_PF1Th1","hist_hibin_PF1Th1",201,0.0,201.0);
	hist_hibin_PF1Th1->Sumw2();
	TH1D *hist_hibin_PF2Th1 = new TH1D("hist_hibin_PF2Th1","hist_hibin_PF2Th1",201,0.0,201.0);
	hist_hibin_PF2Th1->Sumw2();
	TH1D *hist_hibin_PF3Th1 = new TH1D("hist_hibin_PF3Th1","hist_hibin_PF3Th1",201,0.0,201.0);
	hist_hibin_PF3Th1->Sumw2();
	TH1D *hist_hibin_PF4Th1 = new TH1D("hist_hibin_PF4Th1","hist_hibin_PF4Th1",201,0.0,201.0);
	hist_hibin_PF4Th1->Sumw2();

	TH1D *hist_hibin_PF1Th2 = new TH1D("hist_hibin_PF1Th2","hist_hibin_PF1Th2",201,0.0,201.0);
	hist_hibin_PF1Th2->Sumw2();
	TH1D *hist_hibin_PF2Th2 = new TH1D("hist_hibin_PF2Th2","hist_hibin_PF2Th2",201,0.0,201.0);
	hist_hibin_PF2Th2->Sumw2();
	TH1D *hist_hibin_PF3Th2 = new TH1D("hist_hibin_PF3Th2","hist_hibin_PF3Th2",201,0.0,201.0);
	hist_hibin_PF3Th2->Sumw2();
	TH1D *hist_hibin_PF4Th2 = new TH1D("hist_hibin_PF4Th2","hist_hibin_PF4Th2",201,0.0,201.0);
	hist_hibin_PF4Th2->Sumw2();

	TH1D *hist_hibin_PF1Th3 = new TH1D("hist_hibin_PF1Th3","hist_hibin_PF1Th3",201,0.0,201.0);
	hist_hibin_PF1Th3->Sumw2();
	TH1D *hist_hibin_PF2Th3 = new TH1D("hist_hibin_PF2Th3","hist_hibin_PF2Th3",201,0.0,201.0);
	hist_hibin_PF2Th3->Sumw2();
	TH1D *hist_hibin_PF3Th3 = new TH1D("hist_hibin_PF3Th3","hist_hibin_PF3Th3",201,0.0,201.0);
	hist_hibin_PF3Th3->Sumw2();
	TH1D *hist_hibin_PF4Th3 = new TH1D("hist_hibin_PF4Th3","hist_hibin_PF4Th3",201,0.0,201.0);
	hist_hibin_PF4Th3->Sumw2();

	TH1D *hist_hibin_PF1Th4 = new TH1D("hist_hibin_PF1Th4","hist_hibin_PF1Th4",201,0.0,201.0);
	hist_hibin_PF1Th4->Sumw2();
	TH1D *hist_hibin_PF2Th4 = new TH1D("hist_hibin_PF2Th4","hist_hibin_PF2Th4",201,0.0,201.0);
	hist_hibin_PF2Th4->Sumw2();
	TH1D *hist_hibin_PF3Th4 = new TH1D("hist_hibin_PF3Th4","hist_hibin_PF3Th4",201,0.0,201.0);
	hist_hibin_PF3Th4->Sumw2();
	TH1D *hist_hibin_PF4Th4 = new TH1D("hist_hibin_PF4Th4","hist_hibin_PF4Th4",201,0.0,201.0);
	hist_hibin_PF4Th4->Sumw2();

	TH1D *hist_hibin_PF1Th5 = new TH1D("hist_hibin_PF1Th5","hist_hibin_PF1Th5",201,0.0,201.0);
	hist_hibin_PF1Th5->Sumw2();
	TH1D *hist_hibin_PF2Th5 = new TH1D("hist_hibin_PF2Th5","hist_hibin_PF2Th5",201,0.0,201.0);
	hist_hibin_PF2Th5->Sumw2();
	TH1D *hist_hibin_PF3Th5 = new TH1D("hist_hibin_PF3Th5","hist_hibin_PF3Th5",201,0.0,201.0);
	hist_hibin_PF3Th5->Sumw2();
	TH1D *hist_hibin_PF4Th5 = new TH1D("hist_hibin_PF4Th5","hist_hibin_PF4Th5",201,0.0,201.0);
	hist_hibin_PF4Th5->Sumw2();

	// loop over events
	int nevents = hea_tree->GetEntries(); // number of events
	cout << endl;
	cout << "Total number of events in those files: "<< nevents << endl;
	cout << endl;


	// Start loop over events
	double nev = (double)nevents;
	for(int i = 0; i < nevents; i++){

		hea_tree->GetEntry(i);

		if(i % 10000 == 0){double alpha = (double)i; cout << " Running -> percentage: " << std::setprecision(3) << ((alpha / nev) * 100) << "%" << endl;}
		if(fabs(vertexz) > 15.0) continue;

		hist_hibin_all->Fill(hiBin);

		if(pprimaryVertexFilter != 1){hist_hibin_PV->Fill(hiBin);}
		if(pclusterCompatibilityFilter != 1){hist_hibin_CC->Fill(hiBin);}

		if(pprimaryVertexFilter != 1) continue;
		if(pclusterCompatibilityFilter != 1) continue;
		hist_hibin_PVCC->Fill(hiBin);

		// Tower Maker
		if(numMinHFTower2 > 0){hist_hibin_HF1Th2->Fill(hiBin);}
		if(numMinHFTower2 > 1){hist_hibin_HF2Th2->Fill(hiBin);}
		if(numMinHFTower2 > 2){hist_hibin_HF3Th2->Fill(hiBin);}
		if(numMinHFTower2 > 3){hist_hibin_HF4Th2->Fill(hiBin);}

		if(numMinHFTower3 > 0){hist_hibin_HF1Th3->Fill(hiBin);}
		if(numMinHFTower3 > 1){hist_hibin_HF2Th3->Fill(hiBin);}
		if(numMinHFTower3 > 2){hist_hibin_HF3Th3->Fill(hiBin);}
		if(numMinHFTower3 > 3){hist_hibin_HF4Th3->Fill(hiBin);}

		if(numMinHFTower4 > 0){hist_hibin_HF1Th4->Fill(hiBin);}
		if(numMinHFTower4 > 1){hist_hibin_HF2Th4->Fill(hiBin);}
		if(numMinHFTower4 > 2){hist_hibin_HF3Th4->Fill(hiBin);}
		if(numMinHFTower4 > 3){hist_hibin_HF4Th4->Fill(hiBin);}

		if(numMinHFTower5 > 0){hist_hibin_HF1Th5->Fill(hiBin);}
		if(numMinHFTower5 > 1){hist_hibin_HF2Th5->Fill(hiBin);}
		if(numMinHFTower5 > 2){hist_hibin_HF3Th5->Fill(hiBin);}
		if(numMinHFTower5 > 3){hist_hibin_HF4Th5->Fill(hiBin);}

		// PF Candidates
		if(hiHF_pfle > 1 && min(nCountsHFPlus_pf, nCountsHFMinus_pf) > 0){hist_hibin_PF1Th2->Fill(hiBin);}
		if(hiHF_pfle > 1 && min(nCountsHFPlus_pf, nCountsHFMinus_pf) > 1){hist_hibin_PF2Th2->Fill(hiBin);}
		if(hiHF_pfle > 1 && min(nCountsHFPlus_pf, nCountsHFMinus_pf) > 2){hist_hibin_PF3Th2->Fill(hiBin);}
		if(hiHF_pfle > 1 && min(nCountsHFPlus_pf, nCountsHFMinus_pf) > 3){hist_hibin_PF4Th2->Fill(hiBin);}

		if(hiHF_pfle > 2 && min(nCountsHFPlus_pf, nCountsHFMinus_pf) > 0){hist_hibin_PF1Th2->Fill(hiBin);}
		if(hiHF_pfle > 2 && min(nCountsHFPlus_pf, nCountsHFMinus_pf) > 1){hist_hibin_PF2Th2->Fill(hiBin);}
		if(hiHF_pfle > 2 && min(nCountsHFPlus_pf, nCountsHFMinus_pf) > 2){hist_hibin_PF3Th2->Fill(hiBin);}
		if(hiHF_pfle > 2 && min(nCountsHFPlus_pf, nCountsHFMinus_pf) > 3){hist_hibin_PF4Th2->Fill(hiBin);}

		if(hiHF_pfle > 3 && min(nCountsHFPlus_pf, nCountsHFMinus_pf) > 0){hist_hibin_PF1Th3->Fill(hiBin);}
		if(hiHF_pfle > 3 && min(nCountsHFPlus_pf, nCountsHFMinus_pf) > 1){hist_hibin_PF2Th3->Fill(hiBin);}
		if(hiHF_pfle > 3 && min(nCountsHFPlus_pf, nCountsHFMinus_pf) > 2){hist_hibin_PF3Th3->Fill(hiBin);}
		if(hiHF_pfle > 3 && min(nCountsHFPlus_pf, nCountsHFMinus_pf) > 3){hist_hibin_PF4Th3->Fill(hiBin);}

		if(hiHF_pfle > 4 && min(nCountsHFPlus_pf, nCountsHFMinus_pf) > 0){hist_hibin_PF1Th4->Fill(hiBin);}
		if(hiHF_pfle > 4 && min(nCountsHFPlus_pf, nCountsHFMinus_pf) > 1){hist_hibin_PF2Th4->Fill(hiBin);}
		if(hiHF_pfle > 4 && min(nCountsHFPlus_pf, nCountsHFMinus_pf) > 2){hist_hibin_PF3Th4->Fill(hiBin);}
		if(hiHF_pfle > 4 && min(nCountsHFPlus_pf, nCountsHFMinus_pf) > 3){hist_hibin_PF4Th4->Fill(hiBin);}

		if(hiHF_pfle > 5 && min(nCountsHFPlus_pf, nCountsHFMinus_pf) > 0){hist_hibin_PF1Th5->Fill(hiBin);}
		if(hiHF_pfle > 5 && min(nCountsHFPlus_pf, nCountsHFMinus_pf) > 1){hist_hibin_PF2Th5->Fill(hiBin);}
		if(hiHF_pfle > 5 && min(nCountsHFPlus_pf, nCountsHFMinus_pf) > 2){hist_hibin_PF3Th5->Fill(hiBin);}
		if(hiHF_pfle > 5 && min(nCountsHFPlus_pf, nCountsHFMinus_pf) > 3){hist_hibin_PF4Th5->Fill(hiBin);}

	}

	// Open, write and close the output file
	TFile *MyFile = new TFile(Form("%s", outputfile.c_str()), "RECREATE");
	if(MyFile->IsOpen()) cout << "output file: " << outputfile.c_str() << endl;
	MyFile->cd();
	MyFile->mkdir("centrality_histograms"); 
	hist_hibin_all->Write();
	hist_hibin_PV->Write();
	hist_hibin_CC->Write();
	hist_hibin_PVCC->Write();
	//TowerMaker
	hist_hibin_HF1Th2->Write();
	hist_hibin_HF2Th2->Write();
	hist_hibin_HF3Th2->Write();
	hist_hibin_HF4Th2->Write();
	hist_hibin_HF1Th3->Write();
	hist_hibin_HF2Th3->Write();
	hist_hibin_HF3Th3->Write();
	hist_hibin_HF4Th3->Write();
	hist_hibin_HF1Th4->Write();
	hist_hibin_HF2Th4->Write();
	hist_hibin_HF3Th4->Write();
	hist_hibin_HF4Th4->Write();
	hist_hibin_HF1Th5->Write();
	hist_hibin_HF2Th5->Write();
	hist_hibin_HF3Th5->Write();
	hist_hibin_HF4Th5->Write();
	//PFCandidates
	hist_hibin_PF1Th1->Write();
	hist_hibin_PF2Th1->Write();
	hist_hibin_PF3Th1->Write();
	hist_hibin_PF4Th1->Write();
	hist_hibin_PF1Th2->Write();
	hist_hibin_PF2Th2->Write();
	hist_hibin_PF3Th2->Write();
	hist_hibin_PF4Th2->Write();
	hist_hibin_PF1Th3->Write();
	hist_hibin_PF2Th3->Write();
	hist_hibin_PF3Th3->Write();
	hist_hibin_PF4Th3->Write();
	hist_hibin_PF1Th4->Write();
	hist_hibin_PF2Th4->Write();
	hist_hibin_PF3Th4->Write();
	hist_hibin_PF4Th4->Write();
	hist_hibin_PF1Th5->Write();
	hist_hibin_PF2Th5->Write();
	hist_hibin_PF3Th5->Write();
	hist_hibin_PF4Th5->Write();
	MyFile->Close();

	cout << endl;
	cout << "------------------------------------- DONE --------------------------------------" << endl;
	cout << endl;

	sec_end = clock();

	cout << "========================================" << endl;
	cout << "Total running time: " << (double)(sec_end - sec_start) / CLOCKS_PER_SEC << " [s]" << endl;
	cout << "========================================" << endl;

}

