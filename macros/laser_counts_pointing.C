#include <iostream>
#include "TString.h"
#include "TFile.h"
#include "TSystem.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TNtuple.h"
#include "TTree.h"


void laser_counts_pointing(int run = 190930){ //this input doesn't matter - but maybe we can use it as a tag to append to output files later?

  std::cout << "laser counts - run: " << run << std::endl;

  TNtuple* h_angles = new TNtuple("angles","Laser pointing angles, central membrane","#theta:#phi");
  
  TFile* infile = new TFile("G4sPHENIX_g4svtx_eval_theta0_90_phi0_360_theta0.625_phi0.625_laser1.root");
  TFile* outfile = new TFile("theta_phi_laser.root", "RECREATE");
 

 TH1F* hxy_entries = new TH1F("hxy_entries","Bins multiplicity",150,0,150);
  hxy_entries->SetXTitle("Bins multiplicity");
  hxy_entries->SetYTitle("entries");

  TH3D *associated_laser = (TH3D*)infile->Get("associated_laser");
  associated_laser->GetZaxis()->SetRangeUser(0,5);

  TH2D *hProjectionxy = (TH2D*) associated_laser->Project3D("xy");
  hProjectionxy->SetXTitle("r (mm)");
  hProjectionxy->SetYTitle("#phi");
 
  TH2F* h_phi_theta_laser = new TH2F("h_phi_theta_laser","Corresponding laser angle",91, -0.5, 90.5, 361,-0.5,360.5);
  h_phi_theta_laser->SetXTitle("#theta (deg.)");
  h_phi_theta_laser->SetYTitle("#phi (deg.)");

  Int_t nentries = hProjectionxy->GetEntries();
  Int_t nxbins = hProjectionxy->GetNbinsX();
  Int_t nybins = hProjectionxy->GetNbinsY();
  Int_t nbins =0;

  //  cout << "nxbins:  " << nxbins << " nybins:  " << nybins << endl;                                      //  cout << "nentries:  " << nentries << endl;                                                                                           

    //Loop over entries
  for (Int_t i=0; i<= nxbins; i++)
    {          
      for (Int_t j=0; j<=nybins; j++)
	{ 
	   hxy_entries->Fill(hProjectionxy->GetBinContent(i,j));
          
	   //   if ((hProjectionxy->GetBinContent(i,j)) > 0){ //Turn on this for only considering the angles that hit the bins of  
           //the r vs phi distribution. Turn off to found out the required angles to fill all of the r vs phi space, even considering the blind spot.

	      nbins++;   
	      //    cout << "nbins:  " << nbins << endl;                                                                                     
	      float rh = hProjectionxy->GetXaxis()->GetBinCenter(i)/10.;
              float phih = hProjectionxy->GetYaxis()->GetBinCenter(j);

              //cout<<"rh = "<<rh*10<<" , phih = "<<phih<<endl;        

	      float xh=rh*TMath::Cos(phih), yh=rh*TMath::Sin(phih), zh=0.25; //cm (central membrane hit coordinate)                                                                                                 
	      float x0=57.9556, y0=15.5291, z0=105.5; //cm (laser coordinates)                                                                                                                       

	      TVector3 v(xh-x0, yh-y0, zh-z0);

	      const float rad_deg = 180./M_PI;

	      float phi = v.Phi();
	      float theta = v.Theta();

	      static constexpr float m_phimin = 0;
	      static constexpr float m_phimax = 2.*M_PI;

	      while( phi < m_phimin ) phi += 2.*M_PI;
	      while( phi >= m_phimax ) phi -= 2.*M_PI;

	      theta = M_PI - theta;
              h_phi_theta_laser->Fill(theta*rad_deg,(phi*rad_deg - 15),hProjectionxy->GetBinContent(i,j));
              h_angles->Fill(theta*rad_deg,(phi*rad_deg - 15));
              

	      //  } //Turn on for only considering the angles that produce the r vs phi distribution
	  // cout <<"bin i: "<< i << " bins multiplicity :  " << hProjectionxy->GetBinContent(i,j) << endl;
                                                                         
        }

    } 

  
     h_angles->Write();                                                    
     // h_angles->Draw("#phi:#theta");                                                    
  
  // h_angles->Draw("#phi:#theta","","colz");
   //hxy_entries->Draw();                                                        
   //h_phi_theta_laser->Draw("colz");
 
 //since plot_function is an int, need to return something
 return 0;
 }
