R__LOAD_LIBRARY(/sphenix/user/mmeskowit/Dijet_Analysis/MDC2_Analysis/Xj_Analysis/roounfold/timroounfold/libRooUnfold.so)
void Make_Response_1D_v1000_nocalib_fakes_and_misses_just_leading_properhalf(string infile1 = "Run11_Pythia8_R04_10gev_100k.root"){



#define _USE_MATH_DEFINES

#include <math.h> 
#include <cmath>
#include <iostream>

#include<vector>
#include<array>
    
gStyle->SetOptStat(0);
gStyle->SetOptTitle(0);
        
// Make histograms needed

  //Needed for gaus function
  TRandom3 obj;

  //random number in each event
  TRandom3 Random;
   // Xj Histograms 
    TH1F* hXjTruth = new TH1F("hXjTruth","",25,0,1.5);
    TH1F* hXjReco = new TH1F("hXjReco","",25,0,1.5);

    //JES and JER histogram
    // TH2F* hJES = new TH2F("hJES","",60,0,80,20,0,1);
   
  
    const double_t pt_bins[]={0,5,10,15,20,25,30,35,40,45,50,55,60,65,70,80};
    const int pt_N = sizeof(pt_bins)/sizeof(pt_bins[0]) - 1;
    int N_unfolds = 1000;    

    

//QA Histograms No Cuts
    TH2F* TruthvsRecoEtaNoCut = new TH2F("TruthvsRecoEtaNoCut","",50,-1.2,1.2,50,-1.2,1.2);
    TH2F* TruthvsRecoPhiNoCut = new TH2F("TruthvsRecoPhiNoCut","",50,-TMath::Pi(),TMath::Pi(),50,-TMath::Pi(),TMath::Pi());
    TH2F* TruthvsRecoPtNoCut = new TH2F("TruthvsRecoPtNoCut","",70,0,70,70,0,70);
    TH1F* hTruthPhiNoCut = new TH1F("hTruthPhiNoCut","",50,-TMath::Pi(),TMath::Pi());
    TH1F* hRecoPhiNoCut = new TH1F("hRecoPhiNoCut","",50,-TMath::Pi(),TMath::Pi());
    TH1F* hTruthEtaNoCut = new TH1F("hTruthEtaNoCut","",25,-1.2,1.2);
    TH1F* hRecoEtaNoCut = new TH1F("hRecoEtaNoCut","",25,-1.2,1.2);
    TH2F* ScaleResEtaNoCut = new TH2F("ScaleResEtaNoCut","",50,-1,1,50,-.15,.15);
    TH2F* ScaleResPhiNoCut = new TH2F("ScaleResPhiNoCut","",50,-TMath::Pi(),TMath::Pi(),50,-.15,.15);
    TH2F* ScaleResPtNoCut = new TH2F("ScaleResPtNoCut","",25,10,60,50,-1.1,0.8);
    TH1D* hTruthSpectraNoCut = new TH1D("hTruthSpectraNoCut","",pt_N,pt_bins);
    TH1D* hRecoSpectraNoCut = new TH1D("hRecoSpectraNoCut","",pt_N,pt_bins);

     //defining Meas and Truth Histograms
    TH1D* hRecoPtLead = new TH1D("hRecoPtLead","",pt_N,pt_bins);
    TH1D* hTruthPtLead = new TH1D("hTruthPtLead","",pt_N,pt_bins);

    
    TH1D* hTruthPtLeadNoMisses = new TH1D("hTruthPtLeadNoMisses","",pt_N,pt_bins);

    TH1D* hRecoPtLeadHalf = new TH1D("hRecoPtLeadHalf","",pt_N,pt_bins);
    TH1D* hTruthPtLeadHalf = new TH1D("hTruthPtLeadHalf","",pt_N,pt_bins);
    TH1D* hResponseRecoPtLeadHalf = new TH1D("hResponseRecoPtLeadHalf","",pt_N,pt_bins);
   

    //making response matrices


    RooUnfoldResponse *resp_full = new RooUnfoldResponse("resp_full","");
    RooUnfoldResponse *resp_half = new RooUnfoldResponse("resp_half","");
    RooUnfoldResponse *resp_test = new RooUnfoldResponse("resp_test","");
    
    resp_full->Setup(hTruthSpectraNoCut, hRecoSpectraNoCut);
    resp_half->Setup(hTruthSpectraNoCut, hRecoSpectraNoCut);
    resp_test->Setup(hTruthSpectraNoCut, hRecoSpectraNoCut);

    //histograms for errors

    TH2D* hResponseTruthRecoFull = new TH2D("hResponseTruthRecoFull","",pt_N,pt_bins,pt_N,pt_bins);
    TH2D* hResponseTruthRecoHalf = new TH2D("hResponseTruthRecoHalf","",pt_N,pt_bins,pt_N,pt_bins);

    RooUnfoldResponse* resp_full_err[N_unfolds];  RooUnfoldResponse* resp_half_err[N_unfolds];

    for(int index = 0; index < N_unfolds;index ++){
      resp_full_err[index] = new RooUnfoldResponse("resp_full_err","");
      resp_half_err[index] = new RooUnfoldResponse("resp_half_err","");
      resp_full_err[index]->Setup(hTruthSpectraNoCut, hRecoSpectraNoCut);
      resp_half_err[index]->Setup(hTruthSpectraNoCut, hRecoSpectraNoCut);  
    }

//QA Histograms with cuts (temporarily commented out)
    /*
    TH2F* TruthvsRecoEta = new TH2F("TruthvsRecoEta","",25,-1,1,25,-1,1);
    TH2F* TruthvsRecoPhi = new TH2F("TruthvsRecoPhi","",25,-TMath::Pi(),TMath::Pi(),25,-TMath::Pi(),TMath::Pi());
    TH2F* TruthvsRecoPt = new TH2F("TruthvsRecoPt","",35,0,70,35,0,70);
    TH1F* hTruthPhi = new TH1F("hTruthPhi","",25,-TMath::Pi(),TMath::Pi());
    TH1F* hRecoPhi = new TH1F("hRecoPhi","",25,-TMath::Pi(),TMath::Pi());
    TH1F* hTruthEta = new TH1F("hTruthEta","",25,-1,1);
    TH1F* hRecoEta = new TH1F("hRecoEta","",25,-1,1);
    TH2F* ScaleResEta = new TH2F("ScaleResEta","",25,-1,1,25,-.15,.15);
    TH2F* ScaleResPhi = new TH2F("ScaleResPhi","",25,-TMath::Pi(),TMath::Pi(),25,-.15,.15);
    TH2F* ScaleResPt = new TH2F("ScaleResPt","",25,10,60,50,-1.1,0.8);
    TH1F* hTruthSpectra = new TH1F("hTruthSpectra","",20,0,80);
    TH1F* hRecoSpectra = new TH1F("hRecoSpectra","",20,0,80);
    */

   //create vectors and load in branches
   //connect variables in the code to variables in the tree
  
    int f = 2;
    

     
   TFile* _file0 = TFile::Open(infile1.c_str());
   TTree* tree = (TTree*)_file0->Get("tree");
   
    if (infile1.find("30gev") != std::string::npos) {
     f = 0;
     }
     if (infile1.find("10gev") != std::string::npos) {
     f = 1;
     }

        
	 cout << "f = " << f << endl;
	 cout << "infile is: " << infile1 << endl;

   std::vector<double> *tjet_pt ={0}; std::vector<double> *tjet_phi ={0}; std::vector<double> *tjet_eta ={0};
   std::vector<double> *rjet_pt ={0}; std::vector<double> *rjet_phi ={0}; std::vector<double> *rjet_eta ={0};
  
   tree->SetBranchAddress("tjet_pt",&tjet_pt); tree->SetBranchAddress("tjet_phi",&tjet_phi); tree->SetBranchAddress("tjet_eta",&tjet_eta);
   tree->SetBranchAddress("rjet_pt",&rjet_pt); tree->SetBranchAddress("rjet_phi",&rjet_phi); tree->SetBranchAddress("rjet_eta",&rjet_eta);
      
   
 //read in correction
   // TFile *f_calib = new TFile("/sphenix/user/vbailey/analysisclean/analysis/JS-Jet/Calibrations/MC-Calibrations/JES_IsoCorr_NumInv.root","READ");
   // TF1 *corr = (TF1*)f_calib->Get("corrFit_Iso0");


   //create variables used in cuts and calculations
   float entry = 0; float matched = 0; float count = 0; int thirty = 0; int ten = 0;
   
   Double_t tleadpt = 0; float tsubleadpt = 0;
   float tleadphi = 0; float tsubleadphi = 0; float tleadeta = 0; float tsubleadeta = 0;   
   float rleadpt = 0; float rsubleadpt = 0;
   float rleadphi = 0; float rsubleadphi = 0; float rleadeta = 0; float rsubleadeta = 0; 

   float Xjtruth = 0; float Xjreco = 0; float JES = 0; float JES2 =0; float lowpt = 0; float filled = 0;
   float nrj = 0; float ntj = 0; float ntrj = 0; double drleadmatch = 0; double drsubleadmatch = 0;  int etacount = 0; int etamiss = 0;  int highptresp = 0; int highpttest = 0; int endgame = 0; int fake = 0; float miss = 0; int misses = 0;
   float Rvalue = 1.1;
 
     if (infile1.find("R02") != std::string::npos) {
     Rvalue = 0.2;
     }

    else if (infile1.find("R03") != std::string::npos) {
     Rvalue = 0.3;
      }
 
    else if (infile1.find("R04") != std::string::npos) {
     Rvalue = 0.4;
     }
 
    else if (infile1.find("R05") != std::string::npos) {
     Rvalue = 0.5;
     }

     float etacut = 1.1 - Rvalue;
   cout << "Rvalue is: " << Rvalue << "  drmin is: " << Rvalue*0.75 << "  eta cut is: " << etacut << endl;

 //number of entries in tree
    int Nentries = tree->GetEntries(); 
    cout << "Nentries = : " << Nentries << endl;
    for (int i = 0;i<Nentries;i++){

      double  choice = Random.Rndm();

      if (i < 100){
	cout << "choice is: " << choice << endl;
      }
        //loop over entries
        tree->GetEntry(i);
	entry = entry + 1;
	 int ntjets = tjet_pt->size();  int nrjets = rjet_pt->size();
     
	 if (ntjets < 2){
	   if (nrjets >= 2){
	     nrj = nrj + 1;
             }
	   if (nrjets < 2){
	     ntrj = ntrj + 1;
	   }
	   continue;
	 }
	 
	 else if (nrjets < 2){
	   if (ntjets >= 2){
	     ntj = ntj + 1;
           }

	   continue;
	 }

	 else if (ntjets >= 2 && nrjets >= 2){ //open conditional loop
	   count = count + 1;

	  

         for (int j = 0; j < ntjets;j++){ //open ntjets loop

	    
	      
	   
              float tjetpt = tjet_pt->at(j);
	      float tjetphi = tjet_phi->at(j);
	      float tjeteta = tjet_eta->at(j);
	      
	      float drmin_lead =  Rvalue;//*0.75;
	      float drmin_sublead =  Rvalue;//*0.75;
	      Double_t reco_matched_leadpt = 0;
	      float reco_matched_leadphi = 0;
	      float reco_matched_leadeta = 0;
	      float reco_matched_subleadpt = 0;
	      float reco_matched_subleadphi = 0;
	      float reco_matched_subleadeta = 0;

	      /*  if (tjetpt < 5 &&  j != ntjets - 1){ 

			continue;
			}*/

		      if (tjetpt > tleadpt){
			tsubleadpt = tleadpt;
			tleadpt = tjetpt;

			tsubleadphi = tleadphi;
			tleadphi = tjetphi;
	        
			tsubleadeta = tleadeta;
			tleadeta = tjeteta;
		      }
		  
		     else if (tjetpt > tsubleadpt){
			tsubleadpt = tjetpt;
			tsubleadphi = tjetphi;
			tsubleadeta = tjeteta;
		      }

		  

		   
		     
		    
		     
			 
	     if (j == (ntjets - 1)){//last entry of event
		    endgame = endgame + 1;

		     for (int k=0; k < nrjets;k++){ //open nrjets loop
		       
	       		       
		      float rjetpt = rjet_pt->at(k);
		      float rjetphi = rjet_phi->at(k);
		      float rjeteta = rjet_eta->at(k);
		      
		      float detalead = fabs(rjeteta - tleadeta);
		      float dphilead = fabs(rjetphi - tleadphi);
		      float detasublead = fabs(rjeteta - tsubleadeta);
		      float dphisublead = fabs(rjetphi - tsubleadphi);

		       if (dphilead > TMath::Pi())
			{
			  dphilead = 2*TMath::Pi()-dphilead;
			}

		       if (dphisublead > TMath::Pi())
			{
			  dphisublead = 2*TMath::Pi()-dphisublead;
			}
		      float drlead = TMath::Sqrt(dphilead*dphilead + detalead*detalead);
		      float drsublead = TMath::Sqrt(dphisublead*dphisublead + detasublead*detasublead);
		      
		      if (drlead < drmin_lead){
			reco_matched_leadpt = rjetpt;
			   reco_matched_leadphi = rjetphi;
			   reco_matched_leadeta = rjeteta;
			   drmin_lead = drlead;
			   drleadmatch = drleadmatch + 1;
		      }
		      if (drsublead < drmin_sublead){
			reco_matched_subleadpt = rjetpt;
			reco_matched_subleadphi = rjetphi;
			reco_matched_subleadeta = rjeteta;
			drmin_sublead = drsublead;
			 drsubleadmatch = drsubleadmatch + 1;

		      }
		     
		     }//close nrjets loop

		     
		     if(tleadpt == 0){//fake condition
		       fake = fake + 1;

		       if(f == 0){
			 resp_full->Fake(reco_matched_leadpt,2.178*pow(10,-9));
			 resp_half->Fake(reco_matched_leadpt,2.178*pow(10,-9));
			
				       }

		      else if(f == 1){
			 resp_full->Fake(reco_matched_leadpt,3.210*pow(10,-6));
			 resp_half->Fake(reco_matched_leadpt,3.210*pow(10,-6));
		       }

		     }

		     if(abs(tleadeta) < etacut && tleadpt > 5 && reco_matched_leadpt == 0 ){//miss condition
		      
		      

		       if(f == 0  && tleadpt > 30){

			  miss = miss + 1; misses = misses + 1;

			   if (miss < 100){
			 cout << "Miss Found, Entry is: " << i << " tleadpt is: " << tleadpt << " tleadeta is: " << abs(tleadeta) << " reco pt is: " << reco_matched_leadpt <<  endl ;
		       }


			 resp_full->Miss(tleadpt,2.178*pow(10,-9));
			 hTruthPtLead->Fill(tleadpt,2.178*pow(10,-9));

			 
			 for(int k = 0; k < N_unfolds;k ++){
			   resp_full_err[k]->Miss(tleadpt,2.178*pow(10,-9));
	                    }

			 if(choice > 0.50){
			  resp_half->Miss(tleadpt,2.178*pow(10,-9));
			  hTruthPtLeadHalf->Fill(tleadpt,2.178*pow(10,-9));
			  for(int k = 0; k < N_unfolds;k ++){
			    resp_half_err[k]->Miss(tleadpt,2.178*pow(10,-9));
	                    }

			 }
			 }
				       

		      else if(f == 1 && tleadpt > 7 && tleadpt < 30){


			  miss = miss + 1; misses = misses + 1;

			   if (miss < 100){
			 cout << "Miss Found, Entry is: " << i << " tleadpt is: " << tleadpt << " tleadeta is: " << abs(tleadeta) << " reco pt is: " << reco_matched_leadpt <<  endl ;
		       } 

			 resp_full->Miss(tleadpt,3.210*pow(10,-6));
			 hTruthPtLead->Fill(tleadpt,3.210*pow(10,-6));

			  for(int k = 0; k < N_unfolds;k ++){
			    resp_full_err[k]->Miss(tleadpt,3.210*pow(10,-6));
	                    }
			  
			  if( choice > 0.50){
			 resp_half->Miss(tleadpt,3.210*pow(10,-6));
			 hTruthPtLeadHalf->Fill(tleadpt,3.210*pow(10,-6));
			 
			   for(int k = 0; k < N_unfolds;k ++){
			     resp_half_err[k]->Miss(tleadpt,3.210*pow(10,-6));
	                    }

			  }

			 
		       }

		     }
		     
		     if (abs(reco_matched_leadeta) > etacut){///etamiss condition start        
		    etamiss = etamiss + 1;
		      
		     }
		    
		  if (abs(reco_matched_leadeta) < etacut && abs(tleadeta) < etacut){///eta condition start        
		    etacount = etacount + 1;
		   
		    if (tleadpt > 5 && reco_matched_leadpt > 5)//matching
			{
			  matched = matched + 1;
			 
			  if( i < 10){
			  cout << "tleadpt is " << tleadpt << endl;
			  }


			  if(f == 0 && tleadpt > 30 ){ // 30 gev matching
			    filled = filled + 1;
			    thirty = thirty + 1;
			    Xjtruth = tsubleadpt/tleadpt;
			    Xjreco = reco_matched_subleadpt/reco_matched_leadpt;
			    // JES = reco_matched_leadpt/tleadpt;
			    // JES2 = reco_matched_subleadpt/tsubleadpt;		    
			    //Fill Histograms
			    hXjTruth->Fill(Xjtruth, 2.178*pow(10,-9));
			    hXjReco->Fill(Xjreco, 2.178*pow(10,-9));
			    // hJES->Fill(tleadpt,JES);
			    // hJES->Fill(tsubleadpt,JES2);
			    hResponseTruthRecoFull->Fill(reco_matched_leadpt,tleadpt, 2.178*pow(10,-9));
			    // hResponseTruthRecoFull->Fill(reco_matched_subleadpt,tsubleadpt, 2.178*pow(10,-9));
			    // hResponseTruthRecoFull->Fill(tleadpt,reco_matched_leadpt, 2.178*pow(10,-9));
			    // hResponseTruthRecoFull->Fill(tsubleadpt,reco_matched_subleadpt, 2.178*pow(10,-9));
			    resp_full->Fill(reco_matched_leadpt,tleadpt, 2.178*pow(10,-9));
			    hRecoPtLead->Fill(reco_matched_leadpt,2.178*pow(10,-9));
			    hTruthPtLead->Fill(tleadpt,2.178*pow(10,-9));
			    hTruthPtLeadNoMisses->Fill(tleadpt,2.178*pow(10,-9));
			    
			    if(choice > 0.50){ 
			     hResponseTruthRecoHalf->Fill(reco_matched_leadpt,tleadpt, 2.178*pow(10,-9));
			     resp_half->Fill(reco_matched_leadpt,tleadpt, 2.178*pow(10,-9));
			    // resp->Fill(tleadpt,reco_matched_leadpt, 2.178*pow(10,-9));
			    hResponseRecoPtLeadHalf->Fill(reco_matched_leadpt,2.178*pow(10,-9));
			    hTruthPtLeadHalf->Fill(tleadpt,2.178*pow(10,-9));
			    }

			    else if(choice < 0.50){
			    resp_test->Fill(reco_matched_leadpt,tleadpt, 2.178*pow(10,-9));
			    // resp_test->Fill(tleadpt,reco_matched_leadpt, 2.178*pow(10,-9));
			    hRecoPtLeadHalf->Fill(reco_matched_leadpt,2.178*pow(10,-9));
			   
			   
			    }

			    // resp->Fill(reco_matched_subleadpt,tsubleadpt, 2.178*pow(10,-9));
			    //fill QA Histograms
			    

			    TruthvsRecoEtaNoCut->Fill(tleadeta,reco_matched_leadeta,2.178*pow(10,-9));
			    TruthvsRecoPhiNoCut->Fill(tleadphi,reco_matched_leadphi, 2.178*pow(10,-9));
			    TruthvsRecoPtNoCut->Fill(tleadpt,reco_matched_leadpt, 2.178*pow(10,-9));
			    /* ScaleResEtaNoCut->Fill(tleadeta,detalead);
			    ScaleResPhiNoCut->Fill(tleadphi,dphilead);
			    ScaleResPtNoCut->Fill(tleadpt,diffpt);*/
			    hTruthPhiNoCut->Fill(tleadphi, 2.178*pow(10,-9));
			    hRecoPhiNoCut->Fill(reco_matched_leadphi, 2.178*pow(10,-9));
			    hTruthEtaNoCut->Fill(tleadeta, 2.178*pow(10,-9));
			    hTruthEtaNoCut->Fill(tsubleadeta, 2.178*pow(10,-9));
			    hRecoEtaNoCut->Fill(reco_matched_leadeta, 2.178*pow(10,-9));
			    hRecoEtaNoCut->Fill(reco_matched_subleadeta, 2.178*pow(10,-9));
			    hTruthSpectraNoCut->Fill(tleadpt,2.178*pow(10,-9));
			    hRecoSpectraNoCut->Fill(reco_matched_leadpt, 2.178*pow(10,-9));
			    hRecoSpectraNoCut->Fill(reco_matched_subleadpt, 2.178*pow(10,-9));
			    hRecoPhiNoCut->Fill(reco_matched_subleadphi, 2.178*pow(10,-9));
                            hTruthSpectraNoCut->Fill(tsubleadpt,2.178*pow(10,-9));
			    hTruthPhiNoCut->Fill(tsubleadphi, 2.178*pow(10,-9));
			  }//close 30 gev matching

			  if (f==1 && tleadpt > 7 && tleadpt < 30){//10 gev matching
			   
			    ten = ten + 1;
			    filled = filled + 1;
			    Xjtruth = tsubleadpt/tleadpt;
			    Xjreco = reco_matched_subleadpt/reco_matched_leadpt;
			    // JES = reco_matched_leadpt/tleadpt;
			    // JES2 = reco_matched_subleadpt/tsubleadpt;		    
			    //Fill Histograms
			    hXjTruth->Fill(Xjtruth, 3.210*pow(10,-6));
			    hXjReco->Fill(Xjreco, 3.210*pow(10,-6));
			    // hJES->Fill(tleadpt,JES);
			    // hJES->Fill(tsubleadpt,JES2);
			    hResponseTruthRecoFull->Fill(reco_matched_leadpt,tleadpt, 3.210*pow(10,-6));
			    // hResponseTruthRecoFull->Fill(reco_matched_subleadpt,tsubleadpt, 3.210*pow(10,-6));
			    // hResponseTruthRecoFull->Fill(tleadpt,reco_matched_leadpt, 3.210*pow(10,-6));
			    //  hResponseTruthRecoFull->Fill(tsubleadpt,reco_matched_subleadpt, 3.210*pow(10,-6));
			    
			    resp_full->Fill(reco_matched_leadpt,tleadpt, 3.210*pow(10,-6));
			    hRecoPtLead->Fill(reco_matched_leadpt,3.210*pow(10,-6));
			    hTruthPtLead->Fill(tleadpt,3.210*pow(10,-6));
			    hTruthPtLeadNoMisses->Fill(tleadpt,3.210*pow(10,-6));
			    
			    if(choice > 0.50){
			    hResponseTruthRecoHalf->Fill(reco_matched_leadpt,tleadpt, 3.210*pow(10,-6));
			    resp_half->Fill(reco_matched_leadpt,tleadpt, 3.210*pow(10,-6));
			    hResponseRecoPtLeadHalf->Fill(reco_matched_leadpt,3.210*pow(10,-6));
			    hTruthPtLeadHalf->Fill(tleadpt,3.210*pow(10,-6));
			    // resp->Fill(tleadpt,reco_matched_leadpt, 3.210*pow(10,-6));
			    }

			    else if(choice < 0.50){
			    hRecoPtLeadHalf->Fill(reco_matched_leadpt,3.210*pow(10,-6));
			    resp_test->Fill(reco_matched_leadpt,tleadpt, 3.210*pow(10,-6));
			    // resp_test->Fill(tleadpt,reco_matched_leadpt, 3.210*pow(10,-6));
			    }

			    //fill QA Histograms
			    TruthvsRecoEtaNoCut->Fill(tleadeta,reco_matched_leadeta,3.210*pow(10,-6));
			    TruthvsRecoPhiNoCut->Fill(tleadphi,reco_matched_leadphi, 3.210*pow(10,-6));
			    TruthvsRecoPtNoCut->Fill(tleadpt,reco_matched_leadpt, 3.210*pow(10,-6));
			    /* ScaleResEtaNoCut->Fill(tleadeta,detalead);
			    ScaleResPhiNoCut->Fill(tleadphi,dphilead);
			    ScaleResPtNoCut->Fill(tleadpt,diffpt);*/
			    hTruthPhiNoCut->Fill(tleadphi, 3.210*pow(10,-6));
			    hRecoPhiNoCut->Fill(reco_matched_leadphi, 3.210*pow(10,-6));
			    hTruthEtaNoCut->Fill(tleadeta, 3.210*pow(10,-6));
			    hTruthEtaNoCut->Fill(tsubleadeta, 3.210*pow(10,-6));
			    hRecoEtaNoCut->Fill(reco_matched_leadeta, 3.210*pow(10,-6));
			    hRecoEtaNoCut->Fill(reco_matched_subleadeta, 3.210*pow(10,-6));
			    hTruthSpectraNoCut->Fill(tleadpt,3.210*pow(10,-6));
			    hRecoSpectraNoCut->Fill(reco_matched_leadpt, 3.210*pow(10,-6));
			    hRecoSpectraNoCut->Fill(reco_matched_subleadpt, 3.210*pow(10,-6));
			    hRecoPhiNoCut->Fill(reco_matched_subleadphi, 3.210*pow(10,-6));
                            hTruthSpectraNoCut->Fill(tsubleadpt,3.210*pow(10,-6));
			    hTruthPhiNoCut->Fill(tsubleadphi, 3.210*pow(10,-6));


			     }// close 10 gev matching
			  			   
			     }// close matched loop
		          
		          }//close eta condition loop
		              tleadpt = 0;
		              tsubleadpt = 0;
			  }// close last entry loop   
		    
	   }//close vector loop
	 
	 }//close conditional
	 
    }//close Nentries loop
    
   
    //int maxX =  hResponseTruthRecoFull->GetNbinsX();
    //int maxY = hResponseTruthRecoFull->GetNbinsY();

   /* for(int i = 0; i < pt_N; i++){
      cout << "average p_t bin " << i << " value is " << (pt_bins[i] + pt_bins [i+1])/2 << endl;}*/ //pt average check

     //filling error for unfolding
    /*
   for (int k = 0; k < N_unfolds;k ++){
     for (int i = 1; i < maxX;i++){
       for (int j = 1; j < maxY; j++){
	   // pT values
	   double ptx = ((pt_bins[i] + pt_bins [i-1])/2);
	   double pty = ((pt_bins[j] + pt_bins [j-1])/2);
	   //full weights
	   double weightf = hResponseTruthRecoFull->GetBinContent(i,j);
	   double errf = hResponseTruthRecoFull->GetBinError(i,j);
	   double newweightf = obj.Gaus(weightf,errf);
	   //half weights
	   double weighth = hResponseTruthRecoHalf->GetBinContent(i,j);
	   double errh = hResponseTruthRecoHalf->GetBinError(i,j);
	   double newweighth = obj.Gaus(weighth,errh);
	   //fill response matrices
	   resp_full_err[k]->Fill(ptx,pty,newweightf);
	   resp_half_err[k]->Fill(ptx,pty,newweighth);

	   // if ( k >  5 and k < 10 && i == 7 && j == 7){
	   // cout << " in histogram " << k <<  " for ptx of  " << ptx << " and pty of " << pty << " the full weight is " << newweightf << " the half weight is " << newweighth << endl;} //rand working check
	   
  	
              }// close j loop

            }//close i loop
      }//close k loop
    */
    //Scale Histograms
    /*
     hXjTruth->Scale(1/filled);
     hXjReco->Scale(1/filled);
     hTruthSpectraNoCut->Scale(1/(2*filled));
     hRecoSpectraNoCut->Scale(1/(2*filled)); 
     hTruthPhiNoCut->Scale(1/(2*filled));
     hRecoPhiNoCut->Scale(1/(2*filled));
     hTruthEtaNoCut->Scale(1/(2*filled));
     hRecoEtaNoCut->Scale(1/(2*filled));
    */
     //print statistics
     cout << "Percent of entries where only nrjets passes condition: " << (nrj/entry)*100 << "%" << endl;
     cout << "Percent of entries where only ntjets passes condition: " << (ntj/entry)*100 << "%" << endl;
     cout << "Percent of entries where neither ntjets or  nrjets passes condition: " << (ntrj/entry)*100 << "%" << endl;
     cout << "total entries gone through: " << entry << endl;
     cout << "total entries where both vectors have more than 2 entries: " << count << endl;
     cout << "Percent to pass vector  conditional: " << (count/entry)*100 << "%" << endl;
     cout << " total entries to pass eta cut: " << etacount << " number of eta misses: " << etamiss << " number of dr lead matches: " << drleadmatch << " number of dr sublead matches: " << drsubleadmatch << " times entered the endgame " << endgame <<  endl;
     cout << "total jets matched: " << matched << " total misses : " << miss << " percent miss: " << miss*100/(miss+matched) << " total fakes : " << fake << endl;
     cout << "Percent of jets matched of total entries: " << (matched/entry)*100 << "%   Percent of jets matched that passed conditional: " << (matched/count)*100              << "% " << endl;
     cout <<  "times histograms filled: " << filled << endl;
     cout << "30 gev histos filled: " << thirty  << "  10 gev histos filled: " << ten << endl;
     cout << "high pt in response: " << highptresp << "high pt in test: " << highpttest;
     cout << "last file used is: " <<  _file0->GetName() << endl;

      //Then we would want to save these histograms to a file so that we can Draw them in a macro
     TString outfilename = infile1;
     outfilename.Prepend("Xj_QA_Unfold_1D_nocalib_fakes_and_misses_v1000_justlead_properhalf_");
     TFile *outfile = TFile::Open(outfilename,"RECREATE");
     
     //Write Histograms
      hXjTruth->Write();
      hXjReco->Write();
      // hJES->Write();
     
      //error hists
      hResponseTruthRecoFull->Write();
      hResponseTruthRecoHalf->Write();
      
      for (int index = 0; index < N_unfolds; index ++){
	resp_full_err[index]->Write(Form("resp_full_err%i",index));
	resp_half_err[index]->Write(Form("resp_half_err%i",index));
      }



      //QA Histograms
      
       TruthvsRecoEtaNoCut->Write();
       TruthvsRecoPhiNoCut->Write();
       TruthvsRecoPtNoCut->Write();
       ScaleResEtaNoCut->Write();
       ScaleResPhiNoCut->Write();
       ScaleResPtNoCut->Write();
       hTruthPhiNoCut->Write();
       hRecoPhiNoCut->Write();
       hTruthEtaNoCut->Write();
       hRecoEtaNoCut->Write();
       hTruthSpectraNoCut->Write();
       hRecoSpectraNoCut->Write();
       
       hRecoPtLead->Write();
       hTruthPtLead->Write();
       hRecoPtLeadHalf->Write();
       hTruthPtLeadHalf->Write();

       hResponseRecoPtLeadHalf->Write();
       
       hTruthPtLeadNoMisses->Write();
   


       resp_half->Write();
       resp_full->Write();
       resp_test->Write();
       

   
}//close macro loop
