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
#include <TProfile.h>
#include <TProfile2D.h>

using namespace std;

void dijet_calo_analysis() {

	string filename = "dijet_calo_analysis_detroit_jet10_waveform_wAj.root";
	std::cout << filename << std::endl;
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
  	TH1F* h_pass_spectra = new TH1F("h_pass_spectra","",50,0,50);

  	TH1F* h_et_towards = new TH1F("h_et_towards","",700,-20,50);
  	TH1F* h_et_transverse = new TH1F("h_et_transverse","",700,-20,50);
  	TH1F* h_et_away = new TH1F("h_et_away","",700,-20,50);

  	TH1F* h_ue_towards = new TH1F("h_ue_towards","",700,-20,50);
  	TH1F* h_ue_transverse = new TH1F("h_ue_transverse","",700,-20,50);
  	TH1F* h_ue_away = new TH1F("h_ue_away","",700,-20,50);

  	TH2F* h_ue_2D_total = new TH2F("h_ue_2D_total","",24,-1.1,1.1,32,0,M_PI);
  	TH2F* h_ue_2D_emcal = new TH2F("h_ue_2D_emcal","",24,-1.1,1.1,32,0,M_PI);
  	TH2F* h_ue_2D_ihcal = new TH2F("h_ue_2D_ihcal","",24,-1.1,1.1,32,0,M_PI);
	TH2F* h_ue_2D_ohcal = new TH2F("h_ue_2D_ohcal","",24,-1.1,1.1,32,0,M_PI);
  	TH2F* h_ue_2D_towards = new TH2F("h_ue_2D_towards","",24,-1.1,1.1,32,0,M_PI);
  	TH2F* h_ue_2D_transverse = new TH2F("h_ue_2D_transverse","",24,-1.1,1.1,32,0,M_PI);
  	TH2F* h_ue_2D_away = new TH2F("h_ue_2D_away","",24,-1.1,1.1,32,0,M_PI);

  	TProfile* h_ue_xj_towards = new TProfile("h_ue_xj_towards","",20,0,1);
  	TProfile* h_ue_xj_transverse = new TProfile("h_ue_xj_transverse","",20,0,1);
  	TProfile* h_ue_xj_away = new TProfile("h_ue_xj_away","",20,0,1);

  	TH1F* h_ue_towards_1020 = new TH1F("h_ue_towards_1020","",700,-20,50);
  	TH1F* h_ue_transverse_1020 = new TH1F("h_ue_transverse_1020","",700,-20,50);
  	TH1F* h_ue_away_1020 = new TH1F("h_ue_away_1020","",700,-20,50);

  	TProfile* h_ue_xj_towards_1020 = new TProfile("h_ue_xj_towards_1020","",20,0,1);
  	TProfile* h_ue_xj_transverse_1020 = new TProfile("h_ue_xj_transverse_1020","",20,0,1);
  	TProfile* h_ue_xj_away_1020 = new TProfile("h_ue_xj_away_1020","",20,0,1);

  	TProfile* h_ue_pt_towards = new TProfile("h_ue_pt_towards","",100,0,50);
  	TProfile* h_ue_pt_transverse = new TProfile("h_ue_pt_transverse","",100,0,50);
  	TProfile* h_ue_pt_away = new TProfile("h_ue_pt_away","",100,0,50);

  	TH2F* h_aj_ptavg = new TH2F("h_aj_ptavg","",100,0,100,100,0,1);

 	TChain chain("T");

	const char* inputDirectory = "/sphenix/tg/tg01/jets/egm2153/JetValOutput/";
	//const char* inputDirectory = "/sphenix/user/egm2153/calib_study/jet_detroit_10GeV/";
	//for (int i = 0; i < runlist.size(); i++) {
	//TString wildcardPath = TString::Format("%soutput_%d_*.root", inputDirectory, runlist[i]);
	//TString wildcardPath = TString::Format("%ssim_output_6*.root", inputDirectory);
	//	chain.Add(wildcardPath);
	//}
	//TString wildcardPath = TString::Format("%srunlist_output_012345.root", inputDirectory);
	//TString wildcardPath = TString::Format("%sfull_output.root", inputDirectory);
	//TString wildcardPath = TString::Format("%ssim_calo_cluster_output.root", inputDirectory);
	//TString wildcardPath = TString::Format("%soutput_*.root", inputDirectory);
	//TString wildcardPath = TString::Format("%ssim_run11_nopileup_calo_cluster_output.root", inputDirectory);
	//TString wildcardPath = TString::Format("%ssim_detroit_output.root", inputDirectory);
	//TString wildcardPath = TString::Format("%ssim_detroit_calo_cluster_output.root", inputDirectory);
	//TString wildcardPath = TString::Format("%ssim_detroit_calo_nozero_output.root", inputDirectory);
	//TString wildcardPath = TString::Format("%sphpythia8_detroitUE_output.root",inputDirectory);
	TString wildcardPath = TString::Format("%ssim_detroit_jet10_topocluster_output.root", inputDirectory);

	chain.Add(wildcardPath);
	//TString wildcardPath1 = TString::Format("%soutput_48*.root", inputDirectory);
	//chain.Add(wildcardPath1);
	//TString wildcardPath2 = TString::Format("%soutput_49*.root", inputDirectory);
	//chain.Add(wildcardPath2);
	//chain.Add(wildcardPath1);
	//std::cout << wildcardPath << std::endl;

	bool sim = 1;
	int m_event;
	int nJet;
	float zvtx;
	float deltaeta = 0.0916667;
	float deltaphi = 0.0981748;
	float secteta = 2.2;
	float sectphi = (2.0*M_PI)/3.0;

	vector<int> *triggerVector = nullptr;
	vector<float> *eta = nullptr;
	vector<float> *phi = nullptr;
	vector<float> *e = nullptr;
	vector<float> *pt = nullptr;

	int emcaln = 0;
	float emcale[24576] = {0.0};
	float emcaleta[24576] = {0.0};
	float emcalphi[24576] = {0.0};

	int ihcaln = 0;
	float ihcale[24576] = {0.0};
	float ihcaleta[24576] = {0.0};
	float ihcalphi[24576] = {0.0};

	int ohcaln = 0;
	float ohcale[24576] = {0.0};
	float ohcaleta[24576] = {0.0};
	float ohcalphi[24576] = {0.0};

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
	chain.SetBranchAddress("emcaleta",emcaleta);
	chain.SetBranchAddress("emcalphi",emcalphi);

	chain.SetBranchAddress("ihcaln",&ihcaln);
	chain.SetBranchAddress("ihcale",ihcale);
	chain.SetBranchAddress("ihcaleta",ihcaleta);
	chain.SetBranchAddress("ihcalphi",ihcalphi);

	chain.SetBranchAddress("ohcaln",&ohcaln);
	chain.SetBranchAddress("ohcale",ohcale);
	chain.SetBranchAddress("ohcaleta",ohcaleta);
	chain.SetBranchAddress("ohcalphi",ohcalphi);


	int eventnumber = 0;
	int events = 0;
    Long64_t nEntries = chain.GetEntries();
    std::cout << nEntries << std::endl;
    for (Long64_t entry = 0; entry < nEntries; ++entry) {
    //for (Long64_t entry = 0; entry < 50000; ++entry) {
    
        chain.GetEntry(entry);
    	if (eventnumber % 10000 == 0) cout << "event " << eventnumber << endl;
    	eventnumber++;

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

  		bool negJet = false;
  		for (int i = 0; i < nJet; i++) {
  			if ((*e)[i] < 0) {
  				negJet = true;
  			}
  		}

  		if (!jettrig && !sim) { continue; }
  		if (isnan(zvtx)) { continue; }
  		if (zvtx < -30 || zvtx > 30) { continue; }
  		if (nJet < 1) { continue; }
  		if (negJet) { continue; }

  		h_cut->Fill(3);

  		h_vz->Fill(zvtx);
  		h_njet->Fill(nJet);
  		for (int i = 0; i < nJet; i++) {
  			h_jetspectra->Fill((*pt)[i]);
  		}

  		/*
  		if (nJet > 7 || negJet) {
  			std::cout << "8 or 9 jets in event" << std::endl;
  			for (int i = 0; i < nJet; i++) {
  				std::cout << "Jet " << i << ": pT " << (*pt)[i] << " eta " << (*eta)[i] << " phi " << (*phi)[i] << std::endl;
  			}
  		}
  		*/

  		if (nJet < 2) { continue; }

  		int ind_lead = 0;
  		int ind_sub = 0;
  		float temp_lead = 0;
  		float temp_sub = 0;
  		for (int i = 0; i < nJet; i++) {
  			if ((*pt)[i] > temp_lead) {
  				if (temp_lead != 0) {
  					temp_sub = temp_lead;
  					ind_sub = ind_lead;
  				}
  				temp_lead = (*pt)[i];
  				ind_lead = i;
  			} else if ((*pt)[i] > temp_sub) {
  				temp_sub = (*pt)[i];
  				ind_sub = i;
  			}
  		}

  		if (fabs((*eta)[ind_lead]) > 0.7 || fabs((*eta)[ind_sub]) > 0.7) { continue; }
  		if ((*pt)[ind_lead] < 12 || (*pt)[ind_sub] < 7) { continue; }

  		h_vz->Fill(zvtx);
  		h_dijetspectra->Fill((*pt)[ind_lead]);
  		h_dijetspectra->Fill((*pt)[ind_sub]);

  		TVector3 lead, sub;
    	lead.SetPtEtaPhi((*pt)[ind_lead], (*eta)[ind_lead], (*phi)[ind_lead]);
    	sub.SetPtEtaPhi((*pt)[ind_sub], (*eta)[ind_sub], (*phi)[ind_sub]);

  		h_leadjet->Fill(lead.Pt());
  		h_subjet->Fill(sub.Pt());
  		h_leadphi->Fill(lead.Phi());
  		h_leadeta->Fill(lead.Eta());
  		h_subphi->Fill(sub.Phi());
  		h_subeta->Fill(sub.Eta());
  		h_deltaphi->Fill(lead.DeltaPhi(sub));
  		h_xj->Fill(sub.Pt()/lead.Pt());

  		/*
  		// check the deltaphi = 0 
  		if (lead.DeltaPhi(sub) > -0.1 && lead.DeltaPhi(sub) < 0.1) {

  			std::cout << "dijet dphi = 0: lead.phi " << lead.Phi() << " sub.phi " << sub.Phi() << " lead.pt " << lead.Pt() << " sub.pt " << sub.Pt();
  			std::cout << " lead.eta " << lead.Eta() << " sub.eta " << sub.Eta() << std::endl;

  			h_zerophi_deltaeta->Fill(fabs(lead.Eta()-sub.Eta()));
  			h_zerophi_deltaR->Fill(lead.DeltaR(sub));
  			h_zerophi_leadjet->Fill(lead.Pt());
  			h_zerophi_subjet->Fill(sub.Pt());
  		}
  		*/

  		// require deltaphi > 2.5 for atlas paper
  		// transverse region is [pi/3,2pi/3] from leading jet 
  		if (fabs(lead.DeltaPhi(sub)) > 2.75) {
  			h_pass_deltaphi->Fill(lead.DeltaPhi(sub));
  			h_pass_xj->Fill(sub.Pt()/lead.Pt());
  			h_pass_leadjet->Fill(lead.Pt());
  			h_pass_subjet->Fill(sub.Pt());
  			h_pass_spectra->Fill(lead.Pt());
  			h_pass_spectra->Fill(sub.Pt());
  			h_aj_ptavg->Fill((lead.Pt()+sub.Pt())/2.0, (lead.Pt()-sub.Pt())/(lead.Pt()+sub.Pt()));

  			/*
  			float et_towards = 0;
  			float et_transverse = 0;
  			float et_away = 0;

  			// add in code to loop through calorimeter towers and find sum of ET and UE for towards, transverse and away
  			for (int i = 0; i < emcaln; i++) {
  				//if (emcale[i] < 0.5) { continue; }
  				TVector3 em;
  				em.SetPtEtaPhi(emcale[i]/cosh(emcaleta[i]),emcaleta[i],emcalphi[i]);
  				float dphi = lead.DeltaPhi(em);
  				h_ue_2D_total->Fill(emcaleta[i],dphi,emcale[i]/cosh(emcaleta[i]));
  				h_ue_2D_emcal->Fill(emcaleta[i],dphi,emcale[i]/cosh(emcaleta[i]));
  				if (dphi < M_PI/3.0) {
  					et_towards += emcale[i]/cosh(emcaleta[i]);
  					h_ue_2D_towards->Fill(emcaleta[i],dphi,emcale[i]/cosh(emcaleta[i]));
  				} else if (dphi > M_PI/3.0 && dphi < (2.0*M_PI)/3.0) {
					et_transverse += emcale[i]/cosh(emcaleta[i]);
  					h_ue_2D_transverse->Fill(emcaleta[i],dphi,emcale[i]/cosh(emcaleta[i]));
  				} else if (dphi > (2.0*M_PI)/3.0) {
  					et_away += emcale[i]/cosh(emcaleta[i]);
  					h_ue_2D_away->Fill(emcaleta[i],dphi,emcale[i]/cosh(emcaleta[i]));
  				}
  			}

  			for (int i = 0; i < ihcaln; i++) {
  				//if (ihcale[i] < 0.5) { continue; }
  				TVector3 ih;
  				ih.SetPtEtaPhi(ihcale[i]/cosh(ihcaleta[i]),ihcaleta[i],ihcalphi[i]);
  				float dphi = lead.DeltaPhi(ih);
  				h_ue_2D_total->Fill(ihcaleta[i],dphi,ihcale[i]/cosh(ihcaleta[i]));
  				h_ue_2D_ihcal->Fill(ihcaleta[i],dphi,ihcale[i]/cosh(ihcaleta[i]));
  				if (dphi < M_PI/3.0) {
  					et_towards += ihcale[i]/cosh(ihcaleta[i]);
  					h_ue_2D_towards->Fill(ihcaleta[i],dphi,ihcale[i]/cosh(ihcaleta[i]));
  				} else if (dphi > M_PI/3.0 && dphi < (2.0*M_PI)/3.0) {
					et_transverse += ihcale[i]/cosh(ihcaleta[i]);
  					h_ue_2D_transverse->Fill(ihcaleta[i],dphi,ihcale[i]/cosh(ihcaleta[i]));
  				} else if (dphi > (2.0*M_PI)/3.0) {
  					et_away += ihcale[i]/cosh(ihcaleta[i]);
  					h_ue_2D_away->Fill(ihcaleta[i],dphi,ihcale[i]/cosh(ihcaleta[i]));
  				}
  			}

  			for (int i = 0; i < ohcaln; i++) {
  				//if (ohcale[i] < 0.5) { continue; }
  				TVector3 oh;
  				oh.SetPtEtaPhi(ohcale[i]/cosh(ohcaleta[i]),ohcaleta[i],ohcalphi[i]);
  				float dphi = lead.DeltaPhi(oh);
  				h_ue_2D_total->Fill(ohcaleta[i],dphi,ohcale[i]/cosh(ohcaleta[i]));
  				h_ue_2D_ohcal->Fill(ohcaleta[i],dphi,ohcale[i]/cosh(ohcaleta[i]));
  				if (dphi < M_PI/3.0) {
  					et_towards += ohcale[i]/cosh(ohcaleta[i]);
  					h_ue_2D_towards->Fill(ohcaleta[i],dphi,ohcale[i]/cosh(ohcaleta[i]));
  				} else if (dphi > M_PI/3.0 && dphi < (2.0*M_PI)/3.0) {
					et_transverse += ohcale[i]/cosh(ohcaleta[i]);
  					h_ue_2D_transverse->Fill(ohcaleta[i],dphi,ohcale[i]/cosh(ohcaleta[i]));
  				} else if (dphi > (2.0*M_PI)/3.0) {
  					et_away += ohcale[i]/cosh(ohcaleta[i]);
  					h_ue_2D_away->Fill(ohcaleta[i],dphi,ohcale[i]/cosh(ohcaleta[i]));
  				}
  			}

  			h_et_towards->Fill(et_towards);
  			h_et_transverse->Fill(et_transverse);
  			h_et_away->Fill(et_away);

  			h_ue_towards->Fill(et_towards/(secteta*sectphi));
  			h_ue_transverse->Fill(et_transverse/(secteta*sectphi));
  			h_ue_away->Fill(et_away/(secteta*sectphi));

  			h_ue_xj_towards->Fill(sub.Pt()/lead.Pt(),et_towards);
  			h_ue_xj_transverse->Fill(sub.Pt()/lead.Pt(),et_transverse);
  			h_ue_xj_away->Fill(sub.Pt()/lead.Pt(),et_away);

  			h_ue_pt_towards->Fill(lead.Pt(),et_towards);
  			h_ue_pt_transverse->Fill(lead.Pt(),et_transverse);
  			h_ue_pt_away->Fill(lead.Pt(),et_away);

  			if (lead.Pt() > 10 && lead.Pt() < 20) {
  				h_ue_towards_1020->Fill(et_towards/(secteta*sectphi));
	  			h_ue_transverse_1020->Fill(et_transverse/(secteta*sectphi));
	  			h_ue_away_1020->Fill(et_away/(secteta*sectphi));

	  			h_ue_xj_towards_1020->Fill(sub.Pt()/lead.Pt(),et_towards);
	  			h_ue_xj_transverse_1020->Fill(sub.Pt()/lead.Pt(),et_transverse);
	  			h_ue_xj_away_1020->Fill(sub.Pt()/lead.Pt(),et_away);
  			}
  			*/

  			events++;
  		}

  	}

  	/*
  	h_ue_2D_total->Scale(1.0/(events*deltaeta*deltaphi));
  	h_ue_2D_towards->Scale(1.0/(events*deltaeta*deltaphi));
  	h_ue_2D_transverse->Scale(1.0/(events*deltaeta*deltaphi));
  	h_ue_2D_away->Scale(1.0/(events*deltaeta*deltaphi));

  	h_ue_2D_emcal->Scale(1.0/(events*deltaeta*deltaphi));
  	h_ue_2D_ihcal->Scale(1.0/(events*deltaeta*deltaphi));
  	h_ue_2D_ohcal->Scale(1.0/(events*deltaeta*deltaphi));

  	h_ue_xj_towards->Scale(1.0/(secteta*sectphi));
  	h_ue_xj_transverse->Scale(1.0/(secteta*sectphi));
  	h_ue_xj_away->Scale(1.0/(secteta*sectphi));

  	h_ue_xj_towards_1020->Scale(1.0/(secteta*sectphi));
  	h_ue_xj_transverse_1020->Scale(1.0/(secteta*sectphi));
  	h_ue_xj_away_1020->Scale(1.0/(secteta*sectphi));

  	h_ue_pt_towards->Scale(1.0/(secteta*sectphi));
  	h_ue_pt_transverse->Scale(1.0/(secteta*sectphi));
  	h_ue_pt_away->Scale(1.0/(secteta*sectphi));
  	*/

  	out->Write();
  	out->Close();

}