Int_t laser_counts_all()

{
  /*
  TFile *file1 = TFile::Open("TpcDirectLaserReconstruction_theta0_180_phi0_360_theta10_phi10_laser1.root");
  TH3D *associated_laser1 = (TH3D*)gDirectory->Get("associated_laser");
  associated_laser1->SetName("associated_laser1");
  associated_laser1->GetZaxis()->SetRangeUser(0,5);

  TFile *file2 = TFile::Open("TpcDirectLaserReconstruction_theta0_180_phi0_360_theta5_phi5_laser1.root");
  TH3D *associated_laser2 = (TH3D*)file2->Get("associated_laser");
  associated_laser2->SetName("associated_laser2");
  associated_laser2->GetZaxis()->SetRangeUser(0,5);

  TFile *file3 = TFile::Open("TpcDirectLaserReconstruction_theta0_180_phi0_360_theta2.5_phi2.5_laser1.root");
  TH3D *associated_laser3 = (TH3D*)file3->Get("associated_laser");
  associated_laser3->SetName("associated_laser3");
  associated_laser3->GetZaxis()->SetRangeUser(0,5);

  TFile *file4 = TFile::Open("TpcDirectLaserReconstruction_theta0_180_phi0_360_theta1.25_phi1.25_laser1.root");
  TH3D *associated_laser4 = (TH3D*)file4->Get("associated_laser");
  associated_laser4->SetName("associated_laser4");
  associated_laser4->GetZaxis()->SetRangeUser(0,5);

  TFile *file5 = TFile::Open("TpcDirectLaserReconstruction_theta0_180_phi0_360_theta0.625_phi0.625_laser1.root");
  TH3D *associated_laser5 = (TH3D*)file5->Get("associated_laser");
  associated_laser5->SetName("associated_laser5");
  associated_laser5->GetZaxis()->SetRangeUser(0,5);

  TFile *file6 = TFile::Open("TpcDirectLaserReconstruction_theta0_180_phi0_360_theta0.3125_phi0.3125_laser1.root");
  TH3D *associated_laser6 = (TH3D*)file6->Get("associated_laser");
  associated_laser6->SetName("associated_laser6");
  associated_laser6->GetZaxis()->SetRangeUser(0,5);
*/
  
  TFile *file1 = TFile::Open("G4sPHENIX_g4svtx_eval_theta0_90_phi0_360_theta10_phi10_laser1.root");
  TH3D *associated_laser1 = (TH3D*)gDirectory->Get("associated_laser");
  associated_laser1->SetName("associated_laser1");
  associated_laser1->GetZaxis()->SetRangeUser(0,5);

  TFile *file2 = TFile::Open("G4sPHENIX_g4svtx_eval_theta0_90_phi0_360_theta5_phi5_laser1.root");
  TH3D *associated_laser2 = (TH3D*)file2->Get("associated_laser");
  associated_laser2->SetName("associated_laser2");
  associated_laser2->GetZaxis()->SetRangeUser(0,5);

  TFile *file3 = TFile::Open("G4sPHENIX_g4svtx_eval_theta0_90_phi0_360_theta2.5_phi2.5_laser1.root");
  TH3D *associated_laser3 = (TH3D*)file3->Get("associated_laser");
  associated_laser3->SetName("associated_laser3");
  associated_laser3->GetZaxis()->SetRangeUser(0,5);

  TFile *file4 = TFile::Open("G4sPHENIX_g4svtx_eval_theta0_90_phi0_360_theta1.25_phi1.25_laser1.root");
  TH3D *associated_laser4 = (TH3D*)file4->Get("associated_laser");
  associated_laser4->SetName("associated_laser4");
  associated_laser4->GetZaxis()->SetRangeUser(0,5);

  TFile *file5 = TFile::Open("G4sPHENIX_g4svtx_eval_theta0_90_phi0_360_theta0.625_phi0.625_laser1.root");
  TH3D *associated_laser5 = (TH3D*)file5->Get("associated_laser");
  associated_laser5->SetName("associated_laser5");
  associated_laser5->GetZaxis()->SetRangeUser(0,5);

  /* TFile *file6 = TFile::Open("TpcDirectLaserReconstruction_theta0_180_phi0_360_theta0.3125_phi0.3125_laser1.root");
  TH3D *associated_laser6 = (TH3D*)file6->Get("associated_laser");
  associated_laser6->SetName("associated_laser6");
  associated_laser6->GetZaxis()->SetRangeUser(0,5); */


  if (file1 && file2 && file3 && file4 && file5->IsZombie()) { cout << "Error opening a file" << endl;
    exit(-1);} 
  else {cout << "Files are okay" << endl;}
 
  TH1F* hxy_entries1 = new TH1F("hxy_entries1","Bins multiplicity",60,0,300);
  hxy_entries1->SetFillColor(0); 
  hxy_entries1->SetLineColor(1);
  hxy_entries1->SetLineWidth(6);  
  hxy_entries1->SetXTitle("Bins multiplicity");
  hxy_entries1->SetYTitle("entries");

  TH1F* hxy_entries2 = new TH1F("hxy_entries2","Bins multiplicity",60,0,300);
  hxy_entries2->SetFillColor(0);
  hxy_entries2->SetLineColor(46);
  hxy_entries2->SetLineWidth(6);
 
 TH1F* hxy_entries3 = new TH1F("hxy_entries3","Bins multiplicity",60,0,300);
 hxy_entries3->SetFillColor(0);
 hxy_entries3->SetLineColor(41);
 hxy_entries3->SetLineWidth(6);

 TH1F* hxy_entries4 = new TH1F("hxy_entries4","Bins multiplicity",60,0,300);
 hxy_entries4->SetFillColor(0);
 hxy_entries4->SetLineColor(43);
 hxy_entries4->SetLineWidth(6);

 TH1F* hxy_entries5 = new TH1F("hxy_entries5","Bins multiplicity",60,0,300);
 hxy_entries5->SetFillColor(0);
 hxy_entries5->SetLineColor(38);
 hxy_entries5->SetLineWidth(6);
  
 /* TH1F* hxy_entries6 = new TH1F("hxy_entries6","Bins multiplicity",60,0,600);
 hxy_entries6->SetFillColor(0);
 hxy_entries6->SetLineColor(29);
 hxy_entries6->SetLineWidth(6);*/


  TH2D *hProjectionxy1 = (TH2D*) associated_laser1->Project3D("xy");
  TH2D *hProjectionxy2 = (TH2D*) associated_laser2->Project3D("xy");
  TH2D *hProjectionxy3 = (TH2D*) associated_laser3->Project3D("xy");
  TH2D *hProjectionxy4 = (TH2D*) associated_laser4->Project3D("xy");
  TH2D *hProjectionxy5 = (TH2D*) associated_laser5->Project3D("xy");
  // TH2D *hProjectionxy6 = (TH2D*) associated_laser6->Project3D("xy");


  Int_t nentries1 = hProjectionxy1->GetEntries();
  Int_t nentries2 = hProjectionxy2->GetEntries();
  Int_t nentries3 = hProjectionxy3->GetEntries();
  Int_t nentries4 = hProjectionxy4->GetEntries();
  Int_t nentries5 = hProjectionxy5->GetEntries();
  // Int_t nentries6 = hProjectionxy6->GetEntries();


  Int_t nxbins1 = hProjectionxy1->GetNbinsX();
  Int_t nxbins2 = hProjectionxy2->GetNbinsX();
  Int_t nxbins3 = hProjectionxy3->GetNbinsX();
  Int_t nxbins4 = hProjectionxy4->GetNbinsX();
  Int_t nxbins5 = hProjectionxy5->GetNbinsX();
  // Int_t nxbins6 = hProjectionxy6->GetNbinsX();

  Int_t nybins1 = hProjectionxy1->GetNbinsY();
  Int_t nybins2 = hProjectionxy2->GetNbinsY();
  Int_t nybins3 = hProjectionxy3->GetNbinsY();
  Int_t nybins4 = hProjectionxy4->GetNbinsY();
  Int_t nybins5 = hProjectionxy5->GetNbinsY();
  // Int_t nybins6 = hProjectionxy6->GetNbinsY();

  Int_t nbins1 =0;
  Int_t nbins2 =0;
  Int_t nbins3 =0;
  Int_t nbins4 =0;
  Int_t nbins5 =0;
  //Int_t nbins6 =0;

  //Loop over entries                                                                                    
  for (Int_t i=0; i<= nxbins1; i++)
    {
      for (Int_t j=0; j<=nybins1; j++)
        {
	  hxy_entries1->Fill(hProjectionxy1->GetBinContent(i,j));
	  if ((hProjectionxy1->GetBinContent(i,j)) > 0){
	    nbins1++;	                                                                                                              
	  }                                                                                                                   
        }
    }

  for (Int_t i=0; i<= nxbins2; i++)
    {
      for (Int_t j=0; j<=nybins2; j++)
	{
          hxy_entries2->Fill(hProjectionxy2->GetBinContent(i,j));
          if ((hProjectionxy2->GetBinContent(i,j)) > 0){
            nbins2++;
	  }
	}
    }

  for (Int_t i=0; i<= nxbins3; i++)
    {
      for (Int_t j=0; j<=nybins3; j++)
	{
          hxy_entries3->Fill(hProjectionxy3->GetBinContent(i,j));
          if ((hProjectionxy3->GetBinContent(i,j)) > 0){
            nbins3++;
	  }
	}
    }

  for (Int_t i=0; i<= nxbins4; i++)
    {
      for (Int_t j=0; j<=nybins4; j++)
	{
          hxy_entries4->Fill(hProjectionxy4->GetBinContent(i,j));
          if ((hProjectionxy4->GetBinContent(i,j)) > 0){
            nbins4++;
	  }
	}
    }

  for (Int_t i=0; i<= nxbins5; i++)
    {
      for (Int_t j=0; j<=nybins5; j++)
	{
          hxy_entries5->Fill(hProjectionxy5->GetBinContent(i,j));
          if ((hProjectionxy5->GetBinContent(i,j)) > 0){
            nbins5++;
	  }
	}
    }

  /*for (Int_t i=0; i<= nxbins6; i++)
    {
      for (Int_t j=0; j<=nybins6; j++)
	{
          hxy_entries6->Fill(hProjectionxy6->GetBinContent(i,j));
          if ((hProjectionxy6->GetBinContent(i,j)) > 0){
            nbins6++;
          }
        }
    }*/


 
  TCanvas *c1 = new TCanvas("c1","Bins multiplicity");
  gStyle->SetOptStat(false);
  c1->SetRightMargin(0.1);
  c1->SetTopMargin(0.1);
  c1->SetFillColor(0);
  c1->SetLogy();

  hxy_entries1->SetTitle("#delta#theta = #delta#phi=10");
  hxy_entries2->SetTitle("#delta#theta = #delta#phi=5");
  hxy_entries3->SetTitle("#delta#theta = #delta#phi=2.5");
  hxy_entries4->SetTitle("#delta#theta = #delta#phi=1.25");
  hxy_entries5->SetTitle("#delta#theta = #delta#phi=0.625");
  // hxy_entries6->SetTitle("#delta#theta = #delta#phi=0.3125");
 

  hxy_entries1->Draw();
  hxy_entries2->Draw("same");
  hxy_entries3->Draw("same");
  hxy_entries4->Draw("same");
  hxy_entries5->Draw("same");
  //hxy_entries6->Draw("same");

  c1->BuildLegend();
  c1->SaveAs("anglesteps.pdf");

  cout<<"End of the program."<<endl;
  return 0;
}
