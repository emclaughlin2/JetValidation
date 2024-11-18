R__LOAD_LIBRARY(/sphenix/user/mmeskowit/Dijet_Analysis/MDC2_Analysis/Xj_Analysis/roounfold/timroounfold/libRooUnfold.so)
void Unfold_1D_wError_v1000( string infile01 = "Xj_QA_Unfold_1D_nocalib_fakes_and_misses_v1000_justlead_Run11_Pythia8_R03_10gev_all.root", string infile02 = "Xj_QA_Unfold_1D_nocalib_fakes_and_misses_v1000_justlead_Run11_Pythia8_R03_30gev_all.root" ){
   
  int N_unfolds = 1000; int N_iterations = 4; int N_Toys = 1;
  
  TFile* infile = TFile::Open(infile01.c_str());
  //TFile* infile = TFile::Open("Xj_QA_Unfold_1D_nocalib_fakes_and_misses_Run11_Pythia8_R05_10gev_all.root");
  RooUnfoldResponse *resp_full= (RooUnfoldResponse*)infile->Get("resp_full");
  RooUnfoldResponse *resp_half= (RooUnfoldResponse*)infile->Get("resp_half");
  TH1D *FullClosureTruth = (TH1D*)infile->Get("hTruthPtLead");
  TH1D *FullClosureMeas = (TH1D*)infile->Get("hRecoPtLead");
  TH1D *HalfClosureTruth = (TH1D*)infile->Get("hTruthPtLeadHalf");
  TH1D *HalfClosureMeas = (TH1D*)infile->Get("hRecoPtLeadHalf");

  TH1D *FullClosureTruthNoMisses = (TH1D*)infile->Get("hTruthPtLeadNoMisses");

  TH1D* hTruthSpectra = (TH1D*)infile->Get("hTruthSpectraNoCut");
  TH1D* hRecoSpectra = (TH1D*)infile->Get("hRecoSpectraNoCut");
  
  TH2D* hResponseTruthReco = (TH2D*)infile->Get("hResponseTruthRecoFull");

  //read in error response matrices

   RooUnfoldResponse* resp_full_err_[N_unfolds];  RooUnfoldResponse* resp_half_err_[N_unfolds];
  for(int i = 0; i < N_unfolds;i ++){
    resp_full_err_[i] = (RooUnfoldResponse*)infile->Get(TString::Format("resp_full_err%i",i));
    resp_half_err_[i] = (RooUnfoldResponse*)infile->Get(TString::Format("resp_half_err%i",i));
    }

  TFile* infile2 = TFile::Open(infile02.c_str());
  // TFile* infile2 = TFile::Open("Xj_QA_Unfold_1D_nocalib_fakes_and_misses_Run11_Pythia8_R05_30gev_all.root");
  RooUnfoldResponse *resp_full2 = (RooUnfoldResponse*)infile2->Get("resp_full");
  RooUnfoldResponse *resp_half2 = (RooUnfoldResponse*)infile2->Get("resp_half");
  
  TH1D *FullClosureTruth2 = (TH1D*)infile2->Get("hTruthPtLead");
  TH1D *FullClosureMeas2 = (TH1D*)infile2->Get("hRecoPtLead");
  TH1D *HalfClosureTruth2 = (TH1D*)infile2->Get("hTruthPtLeadHalf");
  TH1D *HalfClosureMeas2 = (TH1D*)infile2->Get("hRecoPtLeadHalf");

  TH1D *FullClosureTruthNoMisses2 = (TH1D*)infile2->Get("hTruthPtLeadNoMisses");

  TH1D* hTruthSpectra2 = (TH1D*)infile2->Get("hTruthSpectraNoCut");
  TH1D* hRecoSpectra2 = (TH1D*)infile2->Get("hRecoSpectraNoCut");

  TH2D* hResponseTruthReco2 = (TH2D*)infile2->Get("hResponseTruthRecoFull");

 //read in error response matrices

   RooUnfoldResponse* resp_full_err2_[N_unfolds];  RooUnfoldResponse* resp_half_err2_[N_unfolds];
  for(int i = 0; i < N_unfolds;i ++){
    resp_full_err2_[i] = (RooUnfoldResponse*)infile2->Get(TString::Format("resp_full_err%i",i));
    resp_half_err2_[i] = (RooUnfoldResponse*)infile2->Get(TString::Format("resp_half_err%i",i));
    }

   for(int i = 0; i < N_unfolds;i ++){
     resp_full_err_[i]->Add(*resp_full_err2_[i]);
     resp_half_err_[i]->Add(*resp_half_err2_[i]); 
    }
  
  resp_full->Add(*resp_full2);
  resp_half->Add(*resp_half2);

  FullClosureTruth->Add(FullClosureTruth2);
  FullClosureMeas->Add(FullClosureMeas2);
  HalfClosureTruth->Add(HalfClosureTruth2);
  HalfClosureMeas->Add(HalfClosureMeas2);

  FullClosureTruthNoMisses->Add(FullClosureTruthNoMisses2);

  //jet spectra
  hTruthSpectra->Add(hTruthSpectra2);
  hRecoSpectra->Add(hRecoSpectra2);

  hResponseTruthReco->Add(hResponseTruthReco2);


  // TH1D* h_Meas = (TH1D*)resp_test->Hmeasured();
  // TH1D* h_true = (TH1D*)resp_test->Htruth();
 
//Needed for gaus function
  TRandom obj;
 //Error Bar Loop 1
   const double_t pt_bins[]={0,5,10,15,20,25,30,35,40,45,50,55,60,65,70,80};
   const int pt_N = sizeof(pt_bins)/sizeof(pt_bins[0]) - 1; 

   
   TH1D* hRecoUnfoldFull[N_unfolds];
   TH1D* hRecoUnfoldHalf[N_unfolds];
  for(int index = 0; index <  N_unfolds; index ++){
     hRecoUnfoldFull[index] = new TH1D("hRecoUnfoldFull","",pt_N,pt_bins);
     hRecoUnfoldHalf[index] = new TH1D("hRecoUnfoldHalf","",pt_N,pt_bins);
  }

  int maxX = FullClosureMeas->GetNbinsX();

  for (int k=0; k <  N_unfolds; k++){
    cout << "k = " << k << endl;
   for (int i=1;i <= maxX;i++){  
	double mean = FullClosureMeas->GetBinContent(i);
	double err = FullClosureMeas->GetBinError(i);
	double newmean = obj.Gaus(mean,err);
	hRecoUnfoldFull[k]->SetBinContent(i,newmean);


	//Half Closure
	double mean2 = HalfClosureMeas->GetBinContent(i);
	double err2 = HalfClosureMeas->GetBinError(i);
	double newmean2 = obj.Gaus(mean2,err2);
	hRecoUnfoldHalf[k]->SetBinContent(i,newmean2); 
   }
  }// close k loop

  /* TString outfilename = "Unfolding_Hists_Pythia_Run7_1D.root";
   TFile *f_out = new TFile(outfilename,"RECREATE");
  for(int index = 0; index <  N_unfolds; index ++){
     hRecoUnfoldFull[index]->Write(Form("hRecoUnfoldFull%i",index));
     hRecoUnfoldHalf[index]->Write(Form("hRecoUnfoldHalf%i",index));
     }*/



 
  RooUnfoldBayes* unfoldedfull[N_unfolds]; 
  RooUnfoldBayes* unfoldedhalf[N_unfolds];
  TH1D* unfolded_pt_full[N_unfolds];
  TH1D* unfolded_pt_half[N_unfolds];


  //Unfolding occurs
  for(int index = 0; index < N_unfolds; index ++){
    unfoldedfull[index] = new RooUnfoldBayes(resp_full_err_[index],hRecoUnfoldFull[index],N_iterations);
    unfoldedfull[index]->SetNToys(N_Toys);
    unfolded_pt_full[index] = (TH1D*)unfoldedfull[index]->Hreco(RooUnfold::ErrorTreatment::kCovToy);
    unfolded_pt_full[index]->SetName(Form("unfolded_pt_full%i",index));

    //half
    unfoldedhalf[index] = new RooUnfoldBayes(resp_half_err_[index],hRecoUnfoldHalf[index],N_iterations);
    unfoldedhalf[index]->SetNToys(N_Toys);
    unfolded_pt_half[index] = (TH1D*)unfoldedhalf[index]->Hreco(RooUnfold::ErrorTreatment::kCovToy);
    unfolded_pt_half[index]->SetName(Form("unfolded_pt_half%i",index));
    }

  /* old unfolding methodology
  TH1D* h_Meas = (TH1D*)resp_test->Hmeasured();
  TH1D* h_true = (TH1D*)resp_test->Htruth();
  RooUnfoldBayes *unf1 = new RooUnfoldBayes(resp,h_Meas,4);

  TH1D*  h_unfolded_pt = (TH1D*)unf1->Hreco();

  h_unfolded_pt->SetName("h_unfolded_pt");
  */
  
 //full-closure test histograms
  RooUnfoldBayes *unf1 = new RooUnfoldBayes(resp_full,FullClosureMeas,4);
  unf1->SetNToys(100);
  TH1D*  h_unfolded_pt = (TH1D*)unf1->Hreco(RooUnfold::ErrorTreatment::kCovToy);
  //TH2D*  h_unfolded_pt = (TH2D*)unf1->Hreco();
  h_unfolded_pt->SetName("h_unfolded_pt");
  //half-closure test histograms
  RooUnfoldBayes *unf2 = new RooUnfoldBayes(resp_half,HalfClosureMeas,4);
  unf2->SetNToys(100);
  TH1D*  h_unfolded_pt2 = (TH1D*)unf2->Hreco(RooUnfold::ErrorTreatment::kCovToy);
  h_unfolded_pt2->SetName("h_unfolded_pt2");
  



  TH1D* Final_Unfolded_pt_Half = new TH1D("Final_Unfolded_pt_Half","",pt_N,pt_bins);
  TH1D* Final_Unfolded_pt_Full = new TH1D("Final_Unfolded_pt_Full","",pt_N,pt_bins);
  std::vector<double> bins_half;
  std::vector<double> bins_full;

  for (int i=0;i < maxX;i++){
    
 
      for (int index = 0; index < N_unfolds; index++){

	double val_h = unfolded_pt_half[index]->GetBinContent(i);
	double val_f = unfolded_pt_full[index]->GetBinContent(i);
	bins_half.push_back(val_h);
	bins_full.push_back(val_f);
       
	if(index == N_unfolds - 1){ //Extracting Errors and Means
	  double err_h = TMath::StdDev(bins_half.begin(),bins_half.end() );
	  //double mean_h = TMath::Mean(bins_half.begin(),bins_half.end() );
	  double  mean_h =  h_unfolded_pt2->GetBinContent(i);
	  Final_Unfolded_pt_Half->SetBinContent(i,mean_h);
	  Final_Unfolded_pt_Half->SetBinError(i,err_h);
	 

	  double err_f = TMath::StdDev(bins_full.begin(),bins_full.end() );
	  // double mean_f = TMath::Mean(bins_full.begin(),bins_full.end() );
	  double mean_f = h_unfolded_pt->GetBinContent(i);
	  Final_Unfolded_pt_Full->SetBinContent(i,mean_f);
	  Final_Unfolded_pt_Full->SetBinError(i,err_f);
	  
	  bins_half.clear();
	  bins_full.clear();	  
	}// close extracting

      }// close index
      // cout << "i is: " << i << " Vector Size is: " << bins_half.size() <<" Bin error is : " << err_h << " Bin mean is : " << mean_h << endl;
  }//close i loop


  gStyle->SetOptStat(0);
  TCanvas* can = new TCanvas("can","",900,900);
  can->Divide(2,1,0,0);
  can->cd(1);
  gPad->SetLogy();
  gPad->SetLeftMargin(0.15);
  gPad->SetTopMargin(0.1);
  gPad->SetRightMargin(0.05);
  gPad->SetBottomMargin(0.1);
  FullClosureMeas->SetLineColor(2);
  FullClosureMeas->GetXaxis()->SetTitle("p_{T,1} (GeV)");
  FullClosureMeas->GetYaxis()->SetTitle("Yield [Arb. Units]");
  FullClosureMeas->SetTitle("1D Full Closure Half-Closure Test");
  FullClosureMeas->SetAxisRange(pow(10,-6),20,"Y");
  FullClosureMeas->SetAxisRange(10,59,"X");
  FullClosureMeas->SetMarkerSize(1);
  FullClosureMeas->SetMarkerStyle(21);
  FullClosureMeas->SetMarkerColor(2);
  Final_Unfolded_pt_Full->SetMarkerSize(1);
  Final_Unfolded_pt_Full->SetMarkerStyle(8);
  Final_Unfolded_pt_Full->SetMarkerColor(36);
  FullClosureTruth->SetMarkerSize(1);
  FullClosureTruth->SetMarkerStyle(4);
  FullClosureTruth->SetMarkerColor(3);
  FullClosureMeas->Draw("");
  Final_Unfolded_pt_Full->Draw("same");
  FullClosureTruth->SetLineColor(3);
  FullClosureTruth->Draw("same");
 
  
  TLegend* leg = new TLegend(0.55,0.65,0.85,0.85);
  leg->AddEntry(FullClosureMeas,"Measured","P");
  leg->AddEntry(FullClosureTruth,"Truth","P");
  leg->AddEntry(Final_Unfolded_pt_Full, "Unfolded Distribution", "P");
  leg->Draw();

   TLegend *sleg1 = new TLegend(.1,.15,.65,.45);
  sleg1->SetFillStyle(0);
  sleg1->SetBorderSize(0);
  sleg1->AddEntry("","#it{#bf{sPHENIX}} Internal","");
  sleg1->AddEntry("","PYTHIA 8 p+p #sqrt{s}=200 GeV","");
  sleg1->AddEntry("","anti-#it{k}_{#it{t}} #it{R} = 0.4, |#eta| < 0.7","");
  // sleg1->AddEntry(""," |#phi_{1} - #phi_{2}| > 7#pi/8", "");
  sleg1->Draw();

  can->cd(2);
  gPad->SetLogy();
  gPad->SetLeftMargin(0.15);
  gPad->SetTopMargin(0.1);
  gPad->SetBottomMargin(0.1);
  gPad->SetRightMargin(0.05);
  gPad->SetLogy();
  HalfClosureMeas->SetLineColor(2);
  HalfClosureMeas->GetXaxis()->SetTitle("p_{T,1} (GeV)");
  HalfClosureMeas->GetYaxis()->SetTitle("Yield [Arb. Units]");
  HalfClosureMeas->SetTitle("1D Half-Closure Test");
  HalfClosureMeas->SetAxisRange(pow(10,-6),20,"Y");
  HalfClosureMeas->SetAxisRange(10,59,"X");
  HalfClosureMeas->SetMarkerSize(1);
  HalfClosureMeas->SetMarkerStyle(21);
  HalfClosureMeas->SetMarkerColor(2);
  Final_Unfolded_pt_Half->SetMarkerSize(1);
  Final_Unfolded_pt_Half->SetMarkerStyle(8);
  Final_Unfolded_pt_Half->SetMarkerColor(36);
  HalfClosureTruth->SetMarkerSize(1);
  HalfClosureTruth->SetMarkerStyle(4);
  HalfClosureTruth->SetMarkerColor(3);
  HalfClosureMeas->Draw("");
  Final_Unfolded_pt_Half->Draw("Same");
  HalfClosureTruth->SetLineColor(3);
  HalfClosureTruth->Draw("same");
 

  
  TLegend* leg2 = new TLegend(0.55,0.65,0.85,0.85);
  leg2->AddEntry(HalfClosureMeas,"Measured","P");
  leg2->AddEntry(HalfClosureTruth,"Truth","P");
  leg2->AddEntry(Final_Unfolded_pt_Half, "Unfolded Distribution", "P");
  leg2->Draw();
 
  sleg1->Draw();

  TH1D* h_ratio_unfold_truth = (TH1D*)Final_Unfolded_pt_Full->Clone("h_ratio_unfold_truth");
   h_ratio_unfold_truth->Divide(FullClosureTruth);

  TH1D* h_ratio_unfold_truth2 = (TH1D*)Final_Unfolded_pt_Half->Clone("h_ratio_unfold_truth2");
  h_ratio_unfold_truth2->Divide(HalfClosureTruth);

 TCanvas* can2 = new TCanvas("can2","",900,900);
   can2->Divide(2,1,0,0);
  can2->cd(1);
  //gPad->SetLogz();
  gPad->SetRightMargin(0.05);
  gPad->SetLeftMargin(0.15);
  gPad->SetTopMargin(0.1);
  h_ratio_unfold_truth->SetAxisRange(10, 59,"X");
  h_ratio_unfold_truth->SetAxisRange(0.75,1.25,"Y");
  h_ratio_unfold_truth->GetXaxis()->SetTitle("p_{T,1} (GeV)");
  h_ratio_unfold_truth->GetYaxis()->SetTitle("Unfolded/Truth");
  h_ratio_unfold_truth->SetTitle("1D Full-Closure Test Ratio");
  h_ratio_unfold_truth->SetMarkerSize(1);
  h_ratio_unfold_truth->SetMarkerStyle(21);
  h_ratio_unfold_truth->SetMarkerColor(1);
  h_ratio_unfold_truth->Draw("");
  sleg1->Draw();

  can2->cd(2);
  //gPad->SetLogz();
  gPad->SetTopMargin(0.1);
  gPad->SetLeftMargin(0.15);
  gPad->SetRightMargin(0.05);
  h_ratio_unfold_truth2->SetAxisRange(10,59,"X");
  h_ratio_unfold_truth2->SetAxisRange(0.75,1.25,"Y");
  h_ratio_unfold_truth2->GetXaxis()->SetTitle("p_{T,1} (GeV)");
  h_ratio_unfold_truth2->GetYaxis()->SetTitle("Unfolded/Truth");
  h_ratio_unfold_truth2->SetTitle("1D Half-Closure Test Ratio");
  h_ratio_unfold_truth2->SetMarkerSize(1);
  h_ratio_unfold_truth2->SetMarkerStyle(21);
  h_ratio_unfold_truth2->SetMarkerColor(1);
  h_ratio_unfold_truth2->Draw("");

  TLine *l = new TLine(10,1,59,1);
  l->SetLineStyle(9);
  l->Draw("same");


  sleg1->Draw();

  // Make Error plot for Ratios
  TH1D* ratio_unfold_full_errors = new TH1D("ratio_unfold_full_errors","",pt_N,pt_bins);
  TH1D* ratio_unfold_half_errors = new TH1D("ratio_unfold_half_errors","",pt_N,pt_bins);


  for (int i = 1;i <= maxX;i++){
     double err_ratio_f =  h_ratio_unfold_truth->GetBinError(i);
     double err_ratio_h =  h_ratio_unfold_truth2->GetBinError(i);
     double pt_val = (pt_bins[i] + pt_bins [i-1])/2;

     ratio_unfold_full_errors->Fill(pt_val,err_ratio_f);
     ratio_unfold_half_errors->Fill(pt_val,err_ratio_h);
  }

  TCanvas* can_err = new TCanvas("can_err","",1600,800);
  can_err->Divide(2,1,0,0);
  can_err->cd(1);
  gPad->SetRightMargin(0.05);gPad->SetLeftMargin(0.15);gPad->SetTopMargin(0.1);
  
   ratio_unfold_full_errors->SetMarkerStyle(21);
   ratio_unfold_full_errors->Draw("Hist p");


  can_err->cd(2);
  gPad->SetRightMargin(0.05);gPad->SetLeftMargin(0.15);gPad->SetTopMargin(0.1);
   ratio_unfold_half_errors->SetMarkerStyle(21);
   ratio_unfold_half_errors->Draw("Hist P");





  TLegend* leg3 = new TLegend(0.20,0.6,0.45,0.8);
  leg3->AddEntry(HalfClosureTruth,"Truth","L");
  leg3->AddEntry(Final_Unfolded_pt_Half, "Unfolded Distribution", "L");
  leg3->SetBorderSize(0);

  
 TCanvas* can3 = new TCanvas("can3","",800,800);
  gPad->SetLeftMargin(0.13);
  gPad->SetRightMargin(0.15);
  gPad->SetTopMargin(0.1);
  gPad->SetBottomMargin(0.13);
  gPad->SetLogz();
  hResponseTruthReco->GetXaxis()->SetTitle("Measured p_{T,1}[GeV/c]");
  hResponseTruthReco->GetYaxis()->SetTitle("Truth p_{T,1} Truth [GeV/c]");
  hResponseTruthReco->SetAxisRange(0, 59,"X");
  hResponseTruthReco->SetAxisRange(0, 59,"Y");
  hResponseTruthReco->SetAxisRange(pow(10,-9), 10,"Z");
  hResponseTruthReco->SetTitle("Response Matrix Leading Spectra");
  hResponseTruthReco->Draw("ColZ");

 TCanvas* canex4 = new TCanvas("canex4","",800,800);
  gPad->SetLeftMargin(0.13);
  gPad->SetRightMargin(0.1);
  gPad->SetTopMargin(0.1);
  gPad->SetBottomMargin(0.13);
  gPad->SetLogy();
  hTruthSpectra->SetMarkerSize(1);
  hRecoSpectra->SetMarkerSize(1);
  hTruthSpectra->SetMarkerStyle(8);
  hRecoSpectra->SetMarkerStyle(21);
  hTruthSpectra->SetLineColor(1);
  hRecoSpectra->SetLineColor(2);
  hRecoSpectra->SetMarkerColor(2);
  hRecoSpectra->SetMarkerSize(1.5);
  hTruthSpectra->SetMarkerSize(1.5);
  hTruthSpectra->SetAxisRange(10,59);
  hTruthSpectra->SetAxisRange(pow(10,-9),20,"Y");
  hTruthSpectra->GetXaxis()->SetTitleOffset(1.2);
  hTruthSpectra->GetXaxis()->SetTitle(" p_{T,1} [GeV/c]");
  hTruthSpectra->GetYaxis()->SetTitle("Yield [Arbitrary Units]");
  hTruthSpectra->Draw(""); 
  hRecoSpectra->Draw("Same");
 
  TLegend *sleg2 = new TLegend(.20,.3,.50,.55);
  sleg2->SetFillStyle(0);
  sleg2->SetBorderSize(0);
  sleg2->SetTextSize(0.038);
  sleg2->AddEntry("","#it{#bf{sPHENIX}} Simulation","");
  sleg2->AddEntry("","PYTHIA 8 p+p #sqrt{s}=200 GeV","");
  sleg2->AddEntry("","anti-#it{k}_{#it{t}} #it{R} = 0.4","");
  sleg2->AddEntry("","|#eta| < 0.7","");
  sleg2->AddEntry("","p_{T,1}^{truth} > p_{T,2}^{truth}","");

  sleg2->Draw();

  TLegend *leg4 = new TLegend (0.60,0.70,0.87,0.80);
  leg4->AddEntry(hRecoSpectra   ,"Measured","P");
  leg4->AddEntry(hTruthSpectra   ,"Truth","P");

  leg4->Draw();



  TH1D* Miss_Efficiency = (TH1D*)FullClosureTruthNoMisses->Clone("Miss_Efficiency");
  Miss_Efficiency->Divide(FullClosureTruth);


 TCanvas* can5 = new TCanvas("can5","",800,800);
  //gPad->SetLogz();
  gPad->SetRightMargin(0.05);
  gPad->SetLeftMargin(0.15);
  gPad->SetTopMargin(0.1);
  Miss_Efficiency->SetAxisRange(10, 59,"X");
  Miss_Efficiency->SetAxisRange(0.75,1.25,"Y");
  Miss_Efficiency->GetXaxis()->SetTitle("p_{T,1} (GeV)");
  Miss_Efficiency->GetYaxis()->SetTitle("Matched Events/All Events");
  Miss_Efficiency->SetMarkerSize(1);
  Miss_Efficiency->SetMarkerStyle(21);
  Miss_Efficiency->SetMarkerColor(1);
  Miss_Efficiency->Draw("");
  sleg1->Draw();

  
   //Then we would want to save these histograms to a file so that we can Draw them in a macro
     TString outfilename = infile02;
     outfilename.Prepend("Unfolding_Hists_Compare_1000_");
     TFile *outfile = TFile::Open(outfilename,"RECREATE");
    
     FullClosureMeas->Write();
     FullClosureTruth->Write();
     Final_Unfolded_pt_Full->Write();

     Miss_Efficiency->Write();

     HalfClosureMeas->Write();
     HalfClosureTruth->Write();
     Final_Unfolded_pt_Half->Write();

     h_ratio_unfold_truth->Write();
     h_ratio_unfold_truth2->Write();

     hResponseTruthReco->Write();
     hTruthSpectra->Write();
     hRecoSpectra->Write();
     //write Error Hists
     ratio_unfold_full_errors->Write();
     ratio_unfold_half_errors->Write();
  
  /*
  TH1D* h_ratio_unfold_truth = (TH1D*)h_unfolded_pt->Clone("h_ratio_unfold_truth");
  h_ratio_unfold_truth->Divide(h_true);
  

  TCanvas* can2 = new TCanvas("can2","",800,800);
  h_ratio_unfold_truth->GetXaxis()->SetTitle("p_{T} (GeV)");
  h_ratio_unfold_truth->GetYaxis()->SetTitle("unfolded/truth");
  h_ratio_unfold_truth->SetTitle("1D Half-Closure Test Ratio");
  h_ratio_unfold_truth -> Draw("");

  auto* R = resp->HresponseNoOverflow();

  TCanvas* can3 = new TCanvas("can3","",800,800);
  gPad->SetLogz();
  R->SetAxisRange(pow(10,-9),1,"Z");
  R->Draw("colz");
  */

}
