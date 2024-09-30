#include <iostream>
#include <TH2D.h>
#include <TH1D.h>
#include <TChain.h>
#include <TMath.h>
#include <TTree.h>
#include <TFile.h>
#include <sstream> //std::ostringstsream
#include <fstream> //std::ifstream
#include <iostream> //std::cout, std::endl
#include <cmath>
#include <TGraphErrors.h>
#include <TGraph.h>
#include <TSpectrum.h>
#include <TF1.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>
#include <string>
#include <set>
#include <TVector3.h>
#include <map>
#include <vector>
#include <TDatabasePDG.h>
#include <tuple>

using namespace std;

 	vector<int> runlist = {47289, 47293, 47303, 47306, 47310, 47315, 47316, 47323, 47325, 47330};
 	vector<int> runlist1 = {47289, 47293, 47297, 47298, 47303, 47305, 47306, 47308, 47309, 47310, 47315, 47316, 47323, 47325, 47330};
 	vector<int> runlist2 = {47334, 47360, 47377, 47378, 47379, 47381, 47382, 47391};

void dijet_analysis() {

	string filename = "test_analysis.root";
	TFile *out = new TFile(filename.c_str(),"RECREATE");
	TH1F* h_vz = new TH1F("h_vz","",400, -100, 100);
	TH1F* h_njet = new TH1F("h_njet","",20,0,20);
	TH1F* h_jetspectra = new TH1F("h_jetspectra","",50,0,50);
	TH1F* h_dijetspectra = new TH1F("h_dijetspectra","",50,0,50);
  	TH1F* h_leadjet = new TH1F("h_leadjet","",50,0,50);
  	TH1F* h_subjet = new TH1F("h_subjet","",50,0,50);
  	TH1F* h_leadphi = new TH1F("h_leadphi","",50,-M_PI,M_PI);
  	TH1F* h_leadeta = new TH1F("h_leadeta","",50,-1.1,1.1);
  	TH1F* h_subphi = new TH1F("h_subphi","",50,-M_PI,M_PI);
  	TH1F* h_subeta = new TH1F("h_subeta","",50,-1.1,1.1);
  	TH1F* h_deltaphi = new TH1F("h_deltaphi","",125,-2*M_PI,2*M_PI);
  	TH1F* h_xj = new TH1F("h_xj","",20,0,1);
  	TH1F* h_zerophi_deltaeta = new TH1F("h_zerophi_deltaeta","",50,0,4);
  	TH1F* h_zerophi_deltaR = new TH1F("h_zerophi_deltaR","",50,0,4);
  	TH1F* h_pass_deltaphi = new TH1F("h_pass_deltaphi","",125,-2*M_PI,2*M_PI);
  	TH1F* h_pass_xj = new TH1F("h_pass_xj","",20,0,1);
  	TH1I* h_cut = new TH1I("h_cut","",4,0,4);

  	TH1F* h_zerophi_leadjet = new TH1F("h_zerophi_leadjet","",50,0,50);
  	TH1F* h_zerophi_subjet = new TH1F("h_zerophi_subjet","",50,0,50);
  	TH1F* h_pass_leadjet = new TH1F("h_pass_leadjet","",50,0,50);
  	TH1F* h_pass_subjet = new TH1F("h_pass_subjet","",50,0,50);

  	TH2F* h_2D_emcal = new TH2F("h_2D_emcal","",96,0,96,256,0,256);
  	TH2F* h_2D_ihcal = new TH2F("h_2D_ihcal","",24,0,24,64,0,64);
  	TH2F* h_2D_ohcal = new TH2F("h_2D_ohcal","",24,0,24,64,0,64);

  	TH2F* h_2D_emcal_time = new TH2F("h_2D_emcal_time","",96,0,96,256,0,256);
  	TH2F* h_2D_ihcal_time = new TH2F("h_2D_ihcal_time","",24,0,24,64,0,64);
  	TH2F* h_2D_ohcal_time = new TH2F("h_2D_ohcal_time","",24,0,24,64,0,64);

 	TChain chain("T");

	const char* inputDirectory = "/sphenix/tg/tg01/jets/egm2153/JetValOutput/";
	//for (int i = 0; i < runlist.size(); i++) {
	//	TString wildcardPath = TString::Format("%soutput_%d_*.root", inputDirectory, runlist[i]);
		//TString wildcardPath = TString::Format("%ssim_output_*.root", inputDirectory);
	//	chain.Add(wildcardPath);
	//}
	//string wildcardPath = "full_output.root";
	TString wildcardPath = TString::Format("%soutput_47289_509.root", inputDirectory);
	chain.Add(wildcardPath);

	bool sim = 0;
	int m_event;
	int nJet;
	float zvtx;
	vector<int> *triggerVector = nullptr;
	vector<float> *eta = nullptr;
	vector<float> *phi = nullptr;
	vector<float> *e = nullptr;
	vector<float> *pt = nullptr;

	int emcaln = 0;
	float emcale[24576] = {0.0};
	int emcalieta[24576] = {0};
	int emcaliphi[24576] = {0};
	float emcaltime[24576] = {0.0};

	int ihcaln = 0;
	float ihcale[1536] = {0.0};
	int ihcalieta[1536] = {0};
	int ihcaliphi[1536] = {0};
	float ihcaltime[1536] = {0.0};

	int ohcaln = 0;
	float ohcale[1536] = {0.0};
	int ohcalieta[1536] = {0};
	int ohcaliphi[1536] = {0};
	float ohcaltime[1536] = {0.0};

	chain.SetBranchAddress("m_event",&m_event);
	chain.SetBranchAddress("nJet",&nJet);
	chain.SetBranchAddress("zvtx",&zvtx);
	chain.SetBranchAddress("triggerVector",&triggerVector);
	chain.SetBranchAddress("eta",&eta);
	chain.SetBranchAddress("phi",&phi);
	chain.SetBranchAddress("e",&e);
	chain.SetBranchAddress("pt",&pt);

	chain.SetBranchAddress("emcaln",&emcaln);
	chain.SetBranchAddress("emcale",emcale);
	chain.SetBranchAddress("emcalieta",emcalieta);
	chain.SetBranchAddress("emcaliphi",emcaliphi);
	chain.SetBranchAddress("emcaltime",emcaltime);

	chain.SetBranchAddress("ihcaln",&ihcaln);
	chain.SetBranchAddress("ihcale",ihcale);
	chain.SetBranchAddress("ihcalieta",ihcalieta);
	chain.SetBranchAddress("ihcaliphi",ihcaliphi);
	chain.SetBranchAddress("ihcaltime",ihcaltime);

	chain.SetBranchAddress("ohcaln",&ohcaln);
	chain.SetBranchAddress("ohcale",ohcale);
	chain.SetBranchAddress("ohcalieta",ohcalieta);
	chain.SetBranchAddress("ohcaliphi",ohcaliphi);
	chain.SetBranchAddress("ohcaltime",ohcaltime);

	int eventnumber = 0;
    Long64_t nEntries = chain.GetEntries();
    std::cout << nEntries << std::endl;
    for (Long64_t entry = 0; entry < nEntries; ++entry) {
    //for (Long64_t entry = 0; entry < 50000; ++entry) {
        chain.GetEntry(entry);
    	if (eventnumber % 1000 == 0) cout << "event " << eventnumber << endl;
    	eventnumber++;

    	// additions
    	bool high_pt_jet = false;
    	for (int i = 0; i < nJet; i++) {
    		if ((*pt)[i] > 60) {
    			high_pt_jet = true;
    		}
    	}
    	if (!high_pt_jet) { continue; }

    	for (int i = 0; i < emcaln; i++) {
  			h_2D_emcal->Fill(emcalieta[i],emcaliphi[i],emcale[i]);
  			if (emcale[i] > 0.1) h_2D_emcal_time->Fill(emcalieta[i],emcaliphi[i],emcaltime[i]);
  		}

  		for (int i = 0; i < ihcaln; i++) {
  			h_2D_ihcal->Fill(ihcalieta[i],ihcaliphi[i],ihcale[i]);
  			if (ihcale[i] > 0.01) h_2D_ihcal_time->Fill(ihcalieta[i],ihcaliphi[i],ihcaltime[i]);
  		}

  		for (int i = 0; i < ohcaln; i++) {
  			h_2D_ohcal->Fill(ohcalieta[i],ohcaliphi[i],ohcale[i]);
  			if (ohcale[i] > 0.05) h_2D_ohcal_time->Fill(ohcalieta[i],ohcaliphi[i],ohcaltime[i]);
  		}
  		break;
  		//additions

  		bool jettrig = false;
  		for (auto t : *triggerVector) {
  			if (t >= 16 && t <= 23) {
  				jettrig = true;
  				break;
  			}
  		}

  		if (zvtx < -30 || zvtx > 30) { h_cut->Fill(0); }
  		else if (!jettrig) { h_cut->Fill(1); }
  		else if (nJet < 1) { h_cut->Fill(2); }

  		if (!jettrig && !sim) { continue; }
  		if (isnan(zvtx)) { continue; }
  		if (zvtx < -30 || zvtx > 30) { continue; }
  		if (nJet < 1) { continue; }

  		h_cut->Fill(3);

  		h_vz->Fill(zvtx);
  		h_njet->Fill(nJet);
  		for (int i = 0; i < nJet; i++) {
  			h_jetspectra->Fill((*pt)[i]);
  		}
  		if (nJet > 7) {
  			std::cout << "8 or 9 jets in event" << std::endl;
  			for (int i = 0; i < nJet; i++) {
  				std::cout << "Jet " << i << ": pT " << (*pt)[i] << " eta " << (*eta)[i] << " phi " << (*phi)[i] << std::endl;
  			}
  		}

  		if (nJet != 2) { continue; }
  		if (fabs((*eta)[0]) > 0.7 || fabs((*eta)[1]) > 0.7) { continue; }
  		h_dijetspectra->Fill((*pt)[0]);
  		h_dijetspectra->Fill((*pt)[1]);

  		TVector3 lead, sub;
  		if ((*pt)[0] > (*pt)[1]) {
    		lead.SetPtEtaPhi((*pt)[0], (*eta)[0], (*phi)[0]);
    		sub.SetPtEtaPhi((*pt)[1], (*eta)[1], (*phi)[1]);
  		} else {
    		lead.SetPtEtaPhi((*pt)[1], (*eta)[1], (*phi)[1]);
    		sub.SetPtEtaPhi((*pt)[0], (*eta)[0], (*phi)[0]);
  		}
  		
  		h_leadjet->Fill(lead.Pt());
  		h_subjet->Fill(sub.Pt());
  		h_leadphi->Fill(lead.Phi());
  		h_leadeta->Fill(lead.Eta());
  		h_subphi->Fill(sub.Phi());
  		h_subeta->Fill(sub.Eta());
  		h_deltaphi->Fill(lead.DeltaPhi(sub));
  		h_xj->Fill(sub.Pt()/lead.Pt());

  		// check the deltaphi = 0 
  		if (lead.DeltaPhi(sub) > -0.1 && lead.DeltaPhi(sub) < 0.1) {

  			std::cout << "dijet dphi = 0: lead.phi " << lead.Phi() << " sub.phi " << sub.Phi() << " lead.pt " << lead.Pt() << " sub.pt " << sub.Pt();
  			std::cout << " lead.eta " << lead.Eta() << " sub.eta " << sub.Eta() << std::endl;

  			h_zerophi_deltaeta->Fill(fabs(lead.Eta()-sub.Eta()));
  			h_zerophi_deltaR->Fill(lead.DeltaR(sub));
  			h_zerophi_leadjet->Fill(lead.Pt());
  			h_zerophi_subjet->Fill(sub.Pt());
  		}

  		// require deltaphi > 2.5 for atlas paper
  		// transverse region is [pi/3,2pi/3] from leading jet 
  		if (fabs(lead.DeltaPhi(sub)) > 2.5) {
  			h_pass_deltaphi->Fill(lead.DeltaPhi(sub));
  			h_pass_xj->Fill(sub.Pt()/lead.Pt());
  			h_pass_leadjet->Fill(lead.Pt());
  			h_pass_subjet->Fill(sub.Pt());
  		}

  	}

  	out->Write();
  	out->Close();

}