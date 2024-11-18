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
#include "TRandom.h"
#include "TH1D.h"
//#include "../roounfold/src/RooUnfoldResponse.h"
//#include "../roounfold/src/RooUnfoldBayes.h"
//#include "../roounfold/src/RooUnfoldSvd.h"
//#include "../roounfold/src/RooUnfoldTUnfold.h"

using namespace std;

R__LOAD_LIBRARY(/sphenix/user/egm2153/calib_study/JetValidation/analysis/roounfold/libRooUnfold.so)

void jet_inc_pt_unfolding() {
	
	// event and jet histograms 
  	TH1F* h_pass_deltaphi = new TH1F("h_pass_deltaphi","",125,-2*M_PI,2*M_PI);
  	TH1F* h_pass_xj = new TH1F("h_pass_xj","",20,0,1);
  	TH1F* h_pass_spectra = new TH1F("h_pass_spectra","",50,0,50);
  	TH2F* h_aj_ptavg = new TH2F("h_aj_ptavg","",100,0,100,100,0,1);

  	// should create unfolding histograms and response matrices here
  	//Needed for gaus function
  	TRandom3 obj;

  	//random number in each event
  	TRandom3 Random;

  	float ptmin = 15;
    float ptmax = 45;
    float ptbins = (ptmax - ptmin);
    bool doUnfolding = true;
   	bool sim = true;
   	bool topoclusters = true;
   	bool applyCorr = true;

  	//defining Meas and Truth Histograms
    TH1D* hMeasPT = new TH1D("hMeasPT","",ptbins,ptmin,ptmax);
    TH1D* hTruthPT = new TH1D("hTruthPT","",ptbins,ptmin,ptmax);

    // closure test histograms 
    TH1D* hMeasPTHalf = new TH1D("hMeasPTHalf","",ptbins,ptmin,ptmax);
    TH1D* hTruthPTHalf = new TH1D("hTruthPTHalf","",ptbins,ptmin,ptmax);
    TH1D* hRecoPTHalf = new TH1D("hRecoPTHalf","",ptbins,ptmin,ptmax);

    //making response matrices
    RooUnfoldResponse *resp_full = new RooUnfoldResponse(ptbins,ptmin,ptmax,ptbins,ptmin,ptmax,"resp_full","");
    RooUnfoldResponse *resp_half = new RooUnfoldResponse(ptbins,ptmin,ptmax,ptbins,ptmin,ptmax,"resp_half","");
    RooUnfoldResponse *resp_test = new RooUnfoldResponse(ptbins,ptmin,ptmax,ptbins,ptmin,ptmax,"resp_test","");

    //histograms for errors
    TH2D* hResponseTruthMeasFull = new TH2D("hResponseTruthMeasFull","",ptbins,ptmin,ptmax,ptbins,ptmin,ptmax);
    TH2D* hResponseTruthMeasHalf = new TH2D("hResponseTruthMeasHalf","",ptbins,ptmin,ptmax,ptbins,ptmin,ptmax);

 	TChain chain("T");
	const char* inputDirectory = "/sphenix/tg/tg01/jets/egm2153/JetValOutput/";
	TString wildcardPath = TString::Format("%ssim_truth_jet_output.root", inputDirectory);
	chain.Add(wildcardPath);
	std::cout << wildcardPath << std::endl;

	int m_event;
	int nJet;
	int nTruthJet;
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
	vector<float> *truthEta = nullptr;
	vector<float> *truthPhi = nullptr;
	vector<float> *truthE = nullptr;
	vector<float> *truthPt = nullptr;

	chain.SetBranchAddress("m_event",&m_event);
	chain.SetBranchAddress("nJet",&nJet);
	chain.SetBranchAddress("nTruthJet",&nTruthJet);
	chain.SetBranchAddress("zvtx",&zvtx);
	chain.SetBranchAddress("triggerVector",&triggerVector);
	chain.SetBranchAddress("eta",&eta);
	chain.SetBranchAddress("phi",&phi);
	chain.SetBranchAddress("e",&e);
	chain.SetBranchAddress("pt",&pt);
	chain.SetBranchAddress("truthEta",&truthEta);
	chain.SetBranchAddress("truthPhi",&truthPhi);
	chain.SetBranchAddress("truthE",&truthE);
	chain.SetBranchAddress("truthPt",&truthPt);

	TFile *corrFile;
  	TF1 *correction = new TF1("jet energy correction","1",0,80);
  	if(applyCorr) {
      	corrFile = new TFile("JES_IsoCorr_NumInv.root","READ");
      	corrFile -> cd();
      	correction = (TF1*)corrFile -> Get("corrFit_Iso0");
      	corrFile->Close();
    }

	int eventnumber = 0;
	int events = 0;
    Long64_t nEntries = chain.GetEntries();
    std::cout << nEntries << std::endl;
    for (Long64_t entry = 0; entry < nEntries; ++entry) {
    //for (Long64_t entry = 0; entry < 2000; ++entry) {
        chain.GetEntry(entry);
    	if (eventnumber % 10000 == 0) cout << "event " << eventnumber << endl;
    	eventnumber++;

  		bool negJet = false;
  		for (int i = 0; i < nJet; i++) {
  			if ((*e)[i] < 0) {
  				negJet = true;
  			}
  		}

  		if (isnan(zvtx)) { continue; }
  		if (zvtx < -30 || zvtx > 30) { continue; }
  		if (negJet) { continue; }
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

  		TVector3 lead, sub;
  		if (applyCorr) {
    		lead.SetPtEtaPhi((*pt)[ind_lead]*correction->Eval((*pt)[ind_lead]), (*eta)[ind_lead], (*phi)[ind_lead]);
    		sub.SetPtEtaPhi((*pt)[ind_sub]*correction->Eval((*pt)[ind_sub]), (*eta)[ind_sub], (*phi)[ind_sub]);
    	} else {
    		lead.SetPtEtaPhi((*pt)[ind_lead], (*eta)[ind_lead], (*phi)[ind_lead]);
    		sub.SetPtEtaPhi((*pt)[ind_sub], (*eta)[ind_sub], (*phi)[ind_sub]);
    	}

  		// require deltaphi > 2.75 (7pi/8)
  		// transverse region is [pi/3,2pi/3] from leading jet 
  		if (fabs(lead.DeltaPhi(sub)) > 2.75) {
  			double  choice = Random.Rndm();

  			h_pass_deltaphi->Fill(lead.DeltaPhi(sub));
  			h_pass_xj->Fill(sub.Pt()/lead.Pt());
  			h_pass_spectra->Fill(lead.Pt());
  			h_pass_spectra->Fill(sub.Pt());
  			h_aj_ptavg->Fill((lead.Pt()+sub.Pt())/2.0, (lead.Pt()-sub.Pt())/(lead.Pt()+sub.Pt()));

  			bool found_truth = false;
  			for (int i = 0; i < nTruthJet; i++) {
  				TVector3 truth;
  				truth.SetPtEtaPhi((*truthPt)[i], (*truthEta)[i], (*truthPhi)[i]);
  				if (truth.DeltaR(lead) < 0.3) {
  					hMeasPT->Fill(lead.Pt());
  					hTruthPT->Fill(truth.Pt());
  					resp_full->Fill(truth.Pt(),lead.Pt());
  					found_truth = true;
  					break;
  				}
  			}
  			if (!found_truth) {
  				resp_full->Fake(lead.Pt());
  			}

  			events++;
  		}

  	}

  	RooUnfoldBayes *full_unfold = new RooUnfoldBayes(resp_full,hMeasPT,4);
  	//full_unfold->SetNToys(100);
  	TH1D*  hRecoPT = (TH1D*)full_unfold->Hreco();
  	hRecoPT->SetName("hRecoPT");
	
  	string filename = "jet_pt_unfolding.root";
	std::cout << filename << std::endl;
	TFile *out = new TFile(filename.c_str(),"RECREATE");

  	out->cd();
  	h_pass_deltaphi->Write();
  	h_pass_xj->Write();
  	h_pass_spectra->Write();
  	h_aj_ptavg->Write();
    hMeasPT->Write();
    hTruthPT->Write();
    hRecoPT->Write();
    hMeasPTHalf->Write();
    hTruthPTHalf->Write();
    hRecoPTHalf->Write();
    hResponseTruthMeasFull->Write();
    hResponseTruthMeasHalf->Write();
    resp_full->Write();
    resp_half->Write();
    resp_test->Write();
  	out->Close();

}