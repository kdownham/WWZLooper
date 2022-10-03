void plotHist(TString histname, TString year, int nbins, float min, float max, TString tag){

     TString dir = "output/";
     TString filetype = ".root";
     std::vector<TString> WWZ_samples = {"WWZ","HZJ","ggZH","WWZJets"};

     std::vector<TH1D*> hist_vector_of;
     std::vector<TH1D*> hist_vector_sf;

     TString os_hist = "Opposite Flavor;"+histname+";Events";
     TString ss_hist = "Same Flavor;"+histname+";Events";

     THStack *hs_of = new THStack("hs_of",os_hist);
     THStack *hs_sf = new THStack("hs_sf",ss_hist);

     double yield_of;
     double yield_sf;

     TLegend* legend = new TLegend(0.60,0.60,0.90,0.90,"","NDC");
     legend->SetBorderSize(0);
     legend->SetTextFont(43);
     legend->SetTextAlign(12);
     legend->SetLineColor(1);
     legend->SetLineStyle(1);
     legend->SetLineWidth(1);
     legend->SetFillColor(0);
     legend->SetFillStyle(0);

     TLegend* legend1 = new TLegend(0.60,0.60,0.90,0.90,"","NDC");
     legend1->SetBorderSize(0);
     legend1->SetTextFont(43);
     legend1->SetTextAlign(12);
     legend1->SetLineColor(1);
     legend1->SetLineStyle(1);
     legend1->SetLineWidth(1);
     legend1->SetFillColor(0);
     legend1->SetFillStyle(0);

     TH1D* h = new TH1D("h","h", nbins, min, max);
     TH1D* h1 = new TH1D("h1","h1", nbins, min, max);

     for (auto i : WWZ_samples){
	 std::cout << "Sample = " << i << endl;
	 TFile* file = TFile::Open(dir+"output_"+i+"_"+year+filetype,"READ");
	 
	 TTree* tree_in = (TTree*)file->Get("evt");
	 TH1D* h_aux = new TH1D("h_aux","h_aux", nbins, min, max);
         TH1D* h1_aux = new TH1D("h1_aux","h1_aux", nbins, min, max);
	 
	 bool opp_flav, same_flav;
	 double m_ll, pt_l_Z1, pt_l_Z2, pt_l_W1, pt_l_W2,MT2;
	 double evt_weight;

	 tree_in->SetBranchAddress("evt_weight",&evt_weight);
	 tree_in->SetBranchAddress("opp_flav",&opp_flav);
	 tree_in->SetBranchAddress("same_flav",&same_flav);
	 tree_in->SetBranchAddress("m_ll",&m_ll);
	 tree_in->SetBranchAddress("pt_l_Z1",&pt_l_Z1);
	 tree_in->SetBranchAddress("pt_l_Z2",&pt_l_Z2);
	 tree_in->SetBranchAddress("pt_l_W1",&pt_l_W1);
	 tree_in->SetBranchAddress("pt_l_W2",&pt_l_W2);
	 tree_in->SetBranchAddress("MT2",&MT2);

	 for (Long64_t event = 0; event < tree_in->GetEntries(); ++event){

	     tree_in->GetEntry(event);
	
	     //if ( evt_weight < 0.0 ) continue;
	     //if ( !(MT2 > 0.) ) continue;
	
	     if ( same_flav ){
		  if (histname == "m_ll") h_aux->Fill(m_ll,evt_weight);
		  if (histname == "pt_l_Z1") h_aux->Fill(pt_l_Z1,evt_weight);
		  if (histname == "pt_l_Z2") h_aux->Fill(pt_l_Z2,evt_weight);
		  if (histname == "pt_l_W1") h_aux->Fill(pt_l_W1,evt_weight);
		  if (histname == "pt_l_W2") h_aux->Fill(pt_l_W2,evt_weight);
		  if (histname == "MT2") h_aux->Fill(MT2,evt_weight); 
	     }
	     if ( opp_flav ){
		  if (histname == "m_ll") h1_aux->Fill(m_ll,evt_weight);
                  if (histname == "pt_l_Z1") h1_aux->Fill(pt_l_Z1,evt_weight);
                  if (histname == "pt_l_Z2") h1_aux->Fill(pt_l_Z2,evt_weight);
                  if (histname == "pt_l_W1") h1_aux->Fill(pt_l_W1,evt_weight);
                  if (histname == "pt_l_W2") h1_aux->Fill(pt_l_W2,evt_weight);
		  if (histname == "MT2") h1_aux->Fill(MT2,evt_weight); 
	     }

	 }// Event loop  
 
	 h->Add(h_aux);
	 h1->Add(h1_aux);

     }// File loop

     h->SetFillColor(kBlue);
     h1->SetFillColor(kBlue);

     hs_of->Add(h1);
     hs_sf->Add(h);

     yield_of = h1->Integral();
     yield_sf = h->Integral(); 

     std::string of_yield = std::to_string(yield_of);
     std::string sf_yield = std::to_string(yield_sf);

     std::string of_string = "WWZ ("+of_yield+")";
     std::string sf_string = "WWZ ("+sf_yield+")";

     const char * c = sf_string.c_str();
     const char * c1 = of_string.c_str();

     TCanvas *cs = new TCanvas("cs","cs",10,10,1400,900);

     legend->AddEntry(h,c,"f");

     cs->Divide(2,1);
     cs->cd(1);
     hs_sf->Draw("hist");
     legend->Draw();

     legend1->AddEntry(h1,c1,"f");

     cs->cd(2);
     hs_of->Draw("hist");
     legend1->Draw();

     gSystem->ProcessEvents();
     TImage *img = TImage::Create();

     img->FromPad(cs);
     img->WriteImage("/home/users/kdownham/public_html/WWZ_9_28_22/"+histname+"_"+tag+".png");

}
