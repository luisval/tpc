#include <iostream>
#include "TString.h"
#include "TFile.h"
#include "TSystem.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TNtuple.h"
#include "TTree.h"

//need to declare this here - otherwise code doesn't know about function
int plot_function( TString, TH2F*, int, int );

void files_gem(int run = 190930){ //this input doesn't matter - but maybe we can use it as a tag to append to output files later?

  std::cout << "files_gem - run: " << run << std::endl;

  //create TH2 histogram which we fill with ->GetEntries() from GEMS_hit TH1 histo in each file
  TH2F *max_angles_hit = new TH2F("max_angles_hit","Total Number of GEMS hit for event w/ 8 lasers aimed at #theta, #phi",8,0,80,9,0,350);
  max_angles_hit->SetXTitle("Laser #theta (deg.)");
  max_angles_hit->SetYTitle("Laser #phi (deg.)");

  //set up the loop, here we loop over i (theta angles) and j (phi angles) - important: since declared as int, must use whole numbers
  for( int i:{10,20,30,40,50,60,70,80} ) 
    for( int j:{10,45,90,135,180,225,270,315,350} )
    {
      {
	const TString filename( Form( "TpcDirectLaserReconstruction_theta%i_phi%i.root", i,j ) );
        std::cout << "files_gem current file: " << filename << std::endl;

        plot_function(filename, max_angles_hit, i, j); //pass in the file name, and histogram you want to fill and angles looping over
      }
    }

  max_angles_hit->Draw("colz"); //Draw the histogrram !!!!!!

}

//////////////////////////////Plot function//////////////////////////

//in C++ 11, you need the definition to have default cases or it does not compile
int plot_function( TString filename = "TpcDirectLaserReconstruction_theta10_ph45.root", TH2F* max_angles_hit = new TH2F("max_angles_hit","Total Number of GEMS hit for event w/ 8 lasers aimed at #theta, #phi",3,-0.5,50.5,3,-0.5,50.5), int theta = 10, int phi = 45 ) {

  // get the filename
  TFile *infile = new TFile(filename);
  // get the GEMS_hit histo
  TH1F *GEM_histo = (TH1F*)infile->Get("GEMS_hit");
  
  //fill the max_angles_hit with the GetEntries() from GEMS_hit at the current theta and phi you are looping over !
  max_angles_hit->Fill(theta,phi,GEM_histo->GetEntries());

  //since plot_function is an int, need to return something
 return 0;
 }
