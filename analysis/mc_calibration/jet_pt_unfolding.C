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

bool isInRange(float truthJetPt, float mcWeight);
void getLeadSubleadJet(std::vector<float> *pt, int &ind_lead, int &ind_sub);

void jet_pt_unfolding() {
	
	// event and jet histograms 
  	TH1F* h_pass_deltaphi = new TH1F("h_pass_deltaphi","",125,-2*M_PI,2*M_PI);
  	TH1F* h_pass_xj = new TH1F("h_pass_xj","",20,0,1);
  	TH1F* h_pass_spectra = new TH1F("h_pass_spectra","",50,0,50);
  	TH2F* h_aj_ptavg = new TH2F("h_aj_ptavg","",100,0,100,100,0,1);
  	TProfile* jes_ratio = new TProfile("jes_ratio","",50,0,50);

  	// should create unfolding histograms and response matrices here
  	//Needed for gaus function
  	TRandom3 obj;

  	//random number in each event
  	TRandom3 Random;

  	float ptmin = 15;
    float ptmax = 45;
    float ptbins = (ptmax - ptmin) / 2.0;
    bool doUnfolding = true;
   	bool sim = true;
   	bool topoclusters = true;
   	bool applyCorr = true;
   	float dRMax = 0.3;
   	float ptunfoldmin = 20;

  	//defining Meas and Truth Histograms
    TH1D* hMeasPT = new TH1D("hMeasPT","",ptbins,ptmin,ptmax);
    TH1D* hTruthPT = new TH1D("hTruthPT","",ptbins,ptmin,ptmax);

    // closure test histograms 
    TH1D* hMeasPTHalf = new TH1D("hMeasPTHalf","",ptbins,ptmin,ptmax);
    TH1D* hTruthPTHalf = new TH1D("hTruthPTHalf","",ptbins,ptmin,ptmax);

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
  			if (e->at(i) < 0) {
  				negJet = true;
  			}
  		}

  		if (isnan(zvtx)) { continue; }
  		if (zvtx < -30 || zvtx > 30) { continue; }
  		if (negJet) { continue; }
  		if (nJet < 2) { continue; }

  		int ind_lead = 0;
  		int ind_sub = 0;
  		getLeadSubleadJet(pt, ind_lead, ind_sub);

  		if (fabs(eta->at(ind_lead)) > 0.7 || fabs(eta->at(ind_sub)) > 0.7) { continue; }
  		if (pt->at(ind_lead) < 12 || pt->at(ind_sub) < 7) { continue; }

  		TVector3 lead, sub;
  		lead.SetPtEtaPhi(pt->at(ind_lead), eta->at(ind_lead), phi->at(ind_lead));
  		sub.SetPtEtaPhi(pt->at(ind_sub), eta->at(ind_sub), phi->at(ind_sub));

  		// require deltaphi > 2.75 (7pi/8)
  		// transverse region is [pi/3,2pi/3] from leading jet 
  		if (fabs(lead.DeltaPhi(sub)) > 2.75) {
  			double  choice = Random.Rndm();

  			h_pass_deltaphi->Fill(lead.DeltaPhi(sub));
  			h_pass_xj->Fill(sub.Pt()/lead.Pt());
  			h_pass_spectra->Fill(lead.Pt());
  			h_pass_spectra->Fill(sub.Pt());
  			h_aj_ptavg->Fill((lead.Pt()+sub.Pt())/2.0, (lead.Pt()-sub.Pt())/(lead.Pt()+sub.Pt()));
   
	    	int ind_truth_lead = 0;
	    	int ind_truth_sub = 0;
	    	getLeadSubleadJet(truthPt, ind_truth_lead, ind_truth_sub);

	    	TVector3 truth;
	  		truth.SetPtEtaPhi(truthPt->at(ind_truth_lead), truthEta->at(ind_truth_lead), truthPhi->at(ind_truth_lead));

	      	bool found_match = false;
			for(int i = 0; i < nJet; i++) {
				TVector3 reco;
				if (applyCorr) {
					reco.SetPtEtaPhi(pt->at(i)*correction->Eval(pt->at(i)), eta->at(i), eta->at(i));
				} else {
					reco.SetPtEtaPhi(pt->at(i), eta->at(i), eta->at(i));
				}
	  			if (truth.DeltaR(reco) < dRMax && truth.Pt() > ptunfoldmin && reco.Pt() > ptunfoldmin) {
	  				hMeasPT->Fill(reco.Pt());
	  				hTruthPT->Fill(truth.Pt());
	  				resp_full->Fill(reco.Pt(),truth.Pt());
	  				jes_ratio->Fill(truth.Pt(),reco.Pt()/truth.Pt());
	  				found_match = true;
	  				if (choice > 0.5) {
	  					hTruthPTHalf->Fill(truth.Pt());
	  					resp_half->Fill(reco.Pt(),truth.Pt());
	  				} else {
	  					hMeasPTHalf->Fill(reco.Pt());
	  				}
	  				break;
	  			} else if (truth.DeltaR(reco) < dRMax && truth.Pt() > ptunfoldmin) {
	  				hTruthPT->Fill(truth.Pt());
	  				resp_full->Miss(truth.Pt());
	  				found_match = true;
	  				if (choice > 0.5) {
	  					hTruthPTHalf->Fill(truth.Pt());
	  					resp_half->Miss(truth.Pt());
	  				}
	  				break;
	  			} else if (truth.DeltaR(reco) < dRMax && lead.Pt() > ptunfoldmin) {
	  				hMeasPT->Fill(reco.Pt());
	  				resp_full->Fake(reco.Pt());
	  				found_match = true;
	  				if (choice > 0.5) {
	  					resp_half->Fake(reco.Pt());
	  				} else {
	  					hMeasPTHalf->Fill(reco.Pt());
	  				}
	  				break;
	  			}
			} 
			if (!found_match && truth.Pt() > ptunfoldmin) {
  				hTruthPT->Fill(truth.Pt());
	  			resp_full->Miss(truth.Pt());
	  			found_match = true;
	  			if (choice > 0.5) {
	  				hTruthPTHalf->Fill(truth.Pt());
	  				resp_half->Miss(truth.Pt());
	  			}
			} else if (!found_match && lead.Pt() > ptunfoldmin) {
				resp_full->Fake(lead.Pt());
  				hMeasPT->Fill(lead.Pt());
  				if (choice > 0.5) {
  					resp_half->Fake(lead.Pt());
  				} else {
  					hMeasPTHalf->Fill(lead.Pt());
  				}
			}   

  			events++;
  		}

  	}

  	RooUnfoldBayes *full_unfold = new RooUnfoldBayes(resp_full,hMeasPT,4);
  	//full_unfold->SetNToys(100);
  	TH1D*  hRecoPT = (TH1D*)full_unfold->Hreco();
  	hRecoPT->SetName("hRecoPT");

  	RooUnfoldBayes *half_unfold = new RooUnfoldBayes(resp_half,hMeasPTHalf,4);
  	TH1D* hRecoPTHalf = (TH1D*)half_unfold->Hreco();
  	hRecoPTHalf->SetName("hRecoPTHalf");
	
  	string filename = "jet_pt_unfolding.root";
	std::cout << filename << std::endl;
	TFile *out = new TFile(filename.c_str(),"RECREATE");

  	out->cd();
  	h_pass_deltaphi->Write();
  	h_pass_xj->Write();
  	h_pass_spectra->Write();
  	h_aj_ptavg->Write();
  	jes_ratio->Write();
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

bool isInRange(float truthJetPt, float mcWeight)
{
  float ptCutOff1 = -1; 
  float ptCutOff2 = -1;
  
  
  if(abs(mcWeight/39.06e-3 - 1) < 1e-7)//Minimum bias
    {
      //std::cout << "MB event found" << std::endl;
      ptCutOff1 = 0;
      ptCutOff2 = 14;
    }
  else if(abs(mcWeight/3.210e-6 - 1) < 1e-7)
    {
      //std::cout << "10GeV event found" << std::endl;
      ptCutOff1 = 14;
      ptCutOff2 = 37;
    }
  else if(abs(mcWeight/2.178e-9 - 1) < 1e-7)
    {
      //std::cout << "30GeV event found" << std::endl;
      ptCutOff1 = 37;
      ptCutOff2 = 3000;
    }

  if(truthJetPt < ptCutOff2 && truthJetPt >= ptCutOff1) return true;
  return false;
}
 
void getLeadSubleadJet(std::vector<float> *pt, int &ind_lead, int &ind_sub)
{
  	float temp_lead = 0;
  	float temp_sub = 0;
  	for (int i = 0; i < pt->size(); i++) {
  		if (pt->at(i) > temp_lead) {
  			if (temp_lead != 0) {
  				temp_sub = temp_lead;
  				ind_sub = ind_lead;
  			}
  			temp_lead = pt->at(i);
  			ind_lead = i;
  		} else if (pt->at(i) > temp_sub) {
  			temp_sub = pt->at(i);
  			ind_sub = i;
  		}
  	}
}