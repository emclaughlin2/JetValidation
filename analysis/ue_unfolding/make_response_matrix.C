R__LOAD_LIBRARY(/sphenix/user/egm2153/calib_study/JetValidation/analysis/ue_unfolding/roounfold/libRooUnfold.so)
void make_response_matrix(string infile1 = "Run11_Pythia8_R04_10gev_100k.root"){

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
   
    int N_unfolds = 1000; 
    float etmin = -2;
    float etmax = 10;
    float etbins = (etmax - etmin) * 10;
    
    // histograms to make response matrix 
    TH1D* hTruthSpectraNoCut = new TH1D("hTruthSpectraNoCut","",etbins,etmin,etmax);
    TH1D* hRecoSpectraNoCut = new TH1D("hRecoSpectraNoCut","",etbins,etmin,etmax);

     //defining Meas and Truth Histograms
    TH1D* hRecoPtLead = new TH1D("hRecoPtLead","",etbins,etmin,etmax);
    TH1D* hTruthPtLead = new TH1D("hTruthPtLead","",etbins,etmin,etmax);

    TH1D* hTruthPtLeadNoMisses = new TH1D("hTruthPtLeadNoMisses","",etbins,etmin,etmax);

    // closure test histograms 
    TH1D* hRecoPtLeadHalf = new TH1D("hRecoPtLeadHalf","",etbins,etmin,etmax);
    TH1D* hTruthPtLeadHalf = new TH1D("hTruthPtLeadHalf","",etbins,etmin,etmax);
    TH1D* hResponseRecoPtLeadHalf = new TH1D("hResponseRecoPtLeadHalf","",etbins,etmin,etmax);
   

    //making response matrices
    RooUnfoldResponse *resp_full = new RooUnfoldResponse("resp_full","");
    RooUnfoldResponse *resp_half = new RooUnfoldResponse("resp_half","");
    RooUnfoldResponse *resp_test = new RooUnfoldResponse("resp_test","");
    
    resp_full->Setup(hTruthSpectraNoCut, hRecoSpectraNoCut);
    resp_half->Setup(hTruthSpectraNoCut, hRecoSpectraNoCut);
    resp_test->Setup(hTruthSpectraNoCut, hRecoSpectraNoCut);

    //histograms for errors
    TH2D* hResponseTruthRecoFull = new TH2D("hResponseTruthRecoFull","",etbins,etmin,etmax,etbins,etmin,etmax);
    TH2D* hResponseTruthRecoHalf = new TH2D("hResponseTruthRecoHalf","",etbins,etmin,etmax,etbins,etmin,etmax);

    RooUnfoldResponse* resp_full_err[N_unfolds];  
    RooUnfoldResponse* resp_half_err[N_unfolds];

    for(int index = 0; index < N_unfolds;index ++){
      resp_full_err[index] = new RooUnfoldResponse("resp_full_err","");
      resp_half_err[index] = new RooUnfoldResponse("resp_half_err","");
      resp_full_err[index]->Setup(hTruthSpectraNoCut, hRecoSpectraNoCut);
      resp_half_err[index]->Setup(hTruthSpectraNoCut, hRecoSpectraNoCut);  
    }

   //create vectors and load in branches
   //connect variables in the code to variables in the tree
   TFile* _file0 = TFile::Open(infile1.c_str());
   TTree* tree = (TTree*)_file0->Get("tree");
	 cout << "infile is: " << infile1 << endl;

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

  //create variables used in cuts and calculations
  float entry = 0; float matched = 0; float count = 0;
   
   Double_t tleadpt = 0; float tsubleadpt = 0;
   float tleadphi = 0; float tsubleadphi = 0; float tleadeta = 0; float tsubleadeta = 0;   
   float rleadpt = 0; float rsubleadpt = 0;
   float rleadphi = 0; float rsubleadphi = 0; float rleadeta = 0; float rsubleadeta = 0; 

   float Xjtruth = 0; float Xjreco = 0; float JES = 0; float JES2 =0; float lowpt = 0; float filled = 0;
   float nrj = 0; float ntj = 0; float ntrj = 0; double drleadmatch = 0; double drsubleadmatch = 0;  int etacount = 0; int etamiss = 0;  int highptresp = 0; int highpttest = 0; int endgame = 0; int fake = 0; float miss = 0; int misses = 0;
   float Rvalue = 1.1;
 
 //number of entries in tree
    int Nentries = tree->GetEntries(); 
    cout << "Nentries = : " << Nentries << endl;
    for (int i = 0;i<Nentries;i++) {
      double  choice = Random.Rndm();
      if (i < 100) { cout << "choice is: " << choice << endl; }
      //loop over entries
      tree->GetEntry(i);
			entry = entry + 1;
     	
     	

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
