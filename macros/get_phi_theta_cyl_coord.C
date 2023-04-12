#include <iostream>
#include "TString.h"
#include "TFile.h"
#include "TSystem.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TNtuple.h"
#include "TTree.h"


void get_phi_theta_cyl_coord(){

  float rh=70. ,  phih = 6.;

  float xh=rh*TMath::Cos(phih), yh=rh*TMath::Sin(phih), zh=0.25; //cm
  float x0=57.9556, y0=15.5291, z0=105.5; //cm
 
  TVector3 v(xh-x0, yh-y0, zh-z0); 

  const float rad_deg = 180./M_PI;

  float phi = v.Phi();
  float theta = v.Theta();

 
  static constexpr float m_phimin = 0;
  static constexpr float m_phimax = 2.*M_PI;

  while( phi < m_phimin ) phi += 2.*M_PI;
  while( phi >= m_phimax ) phi -= 2.*M_PI;

  theta = M_PI - theta;
  

  std::cout<<"Phi Pointing Angle: "<< (phi*rad_deg) - 15 <<", Theta Pointing Angle: "<< theta*rad_deg << std::endl;


}
