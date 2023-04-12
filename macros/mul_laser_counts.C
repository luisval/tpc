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



int plot_function(TString, TH2F*,Float_t, Float_t);

void mul_laser_counts(int run = 190930){ //this input doesn't matter - but maybe we can use it as a tag to append to output files later?

  std::cout << "laser counts - run: " << run << std::endl;

  TH2F* hxy_entries = new TH2F("hxy_entries","Saturation",40,0,20,40,0,20);
  hxy_entries->SetXTitle("#delta #phi");
  hxy_entries->SetYTitle("#delta #theta ");
  //////////////////arrray of step sizes in theta and phi///////////////////


  //  Float_t step_sizes[6] = {0.3125, 0.625, 1.25, 2.5 , 5 , 10};                                  
  Float_t step_sizes[5] = {0.625, 1.25, 2.5 , 5 , 10};
                                                                                                               
  ////////////////Loop over theta and phi steps/////////////////////////////
  for (Int_t i = 0; i < 5; i++ )
    {          
      for (Int_t j = 0; j < 5; j++ )
	{ 
	  // cout << "i :  " << i << " j: "<< j<< endl;

	  //  const TString filename( Form( "TpcDirectLaserReconstruction_theta0_180_phi0_360_theta%g_phi%g_laser1.root", step_sizes[i],step_sizes[j] ) );
	  const TString filename( Form( "G4sPHENIX_g4svtx_eval_theta0_90_phi0_360_theta%g_phi%g_laser1.root", step_sizes[i],step_sizes[j] ) );
           std::cout << "current file: " << filename << std::endl;
	  TFile *infile = new TFile(filename);
	  
	   if (infile->IsZombie()) continue;
    
	   plot_function(filename, hxy_entries, step_sizes[i],step_sizes[j]);	                                                                                                                                                                              
                                                                          
        }
    }                                                                                                           
   
      hxy_entries->Draw();                                
   
     TH1D *projh2Xphi = hxy_entries->ProjectionX();
     TH1D *projh2Ytheta = hxy_entries->ProjectionY();

     //   projh2Xphi->Draw();  
     // projh2Ytheta->Draw();
     
      TCanvas *c1 = new TCanvas("c1","#delta #phi and #delta #theta");
      c1->Divide(2);
      c1->SetLogx();
      c1->cd(1);  
      projh2Xphi->Draw();
   
      c1->cd(2);
      projh2Ytheta->Draw();
                
  }
   ///////////////Plot function////////////////////////////////////////////
   int plot_function( TString filename = "TpcDirectLaserReconstruction_theta0_180_phi0_360_theta5_phi5_laser1.root", TH2F* hxy_entries = new TH2F("hxy_entries","Bins > 0", 5,0,10,5,0,10), Float_t i = 0, Float_t j = 0 ) {

     
     // get the filename                                                                         
     TFile *infile = new TFile(filename);
     //if (infile->IsZombie()) continue;    
                                                                               
     TH3D *associated_laser = (TH3D*)infile->Get("associated_laser");
     associated_laser->GetZaxis()->SetRangeUser(0,5);

     TH2D *hProjectionxy = (TH2D*) associated_laser->Project3D("xy");
     hProjectionxy->SetXTitle("r (mm)");
     hProjectionxy->SetYTitle("#phi");
     
     Int_t nentries = hProjectionxy->GetEntries();
     Int_t nxbins = hProjectionxy->GetNbinsX();
     Int_t nybins = hProjectionxy->GetNbinsY();
     Int_t nbins =0;

     std::cout<<"Num x bins in hProjectionxy = "<<nxbins<<" , Num y bins in hProjectionxy = "<<nybins<<std::endl;
 
     //Loop over entries                                                                                                                 
     for (Int_t k=0; k<= nxbins; k++)
       {
	 for (Int_t l=0; l<=nybins; l++)
	   {
	     //   std::cout<<"in the loop, k: "<<k<<" l: "<<l<<std::endl;
	     if ((hProjectionxy->GetBinContent(k,l)) > 0){
	               nbins++;
	     		                                                                                          
	     }

	     // cout <<"bin k: "<< k << " bins multiplicity :  " << hProjectionxy->GetBinContent(k,l) << endl;                             

	   }

       }
     std::cout<<"i = "<<i<<" j = "<<j<<std::endl;
     hxy_entries->Fill(i,j,nbins);                                                                                 


 
  //since plot_function is an int, need to return something
 return 0;
 }
