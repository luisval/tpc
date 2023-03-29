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


void laser_counts(int run = 190930){ //this input doesn't matter - but maybe we can use it as a tag to append to output files later?

  std::cout << "laser counts - run: " << run << std::endl;

  
  TFile* infile = new TFile("TpcDirectLaserReconstruction_theta0_180_phi0_360_theta10_phi10_laser1.root");
 
 TH1F* hxy_entries = new TH1F("hxy_entries","Bins multiplicity",150,0,1500);
  hxy_entries->SetXTitle("Bins multiplicity");
  hxy_entries->SetYTitle("entries");

  TH3D *associated_laser = (TH3D*)infile->Get("associated_laser");
  associated_laser->GetZaxis()->SetRangeUser(0,5);

  TH2D *hProjectionxy = (TH2D*) associated_laser->Project3D("xy");
  hProjectionxy->SetXTitle("r (mm)");
  hProjectionxy->SetYTitle("#phi");
 
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
          
	   if ((hProjectionxy->GetBinContent(i,j)) > 0){
	      nbins++;   
	      //    cout << "nbins:  " << nbins << endl;                                                                                
	    }
	  // cout <<"bin i: "<< i << " bins multiplicity :  " << hProjectionxy->GetBinContent(i,j) << endl;
                                                                         
        }

    } 

                                                       
                                                      
   
   hxy_entries->Draw();                                                        

  //since plot_function is an int, need to return something
 return 0;
 }
