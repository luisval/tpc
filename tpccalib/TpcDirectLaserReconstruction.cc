/**
 * \file TpcDirectLaserReconstruction.cc
 * \brief performs the reconstruction of TPC direct laser tracks
 * \author Hugo Pereira Da Costa <hugo.pereira-da-costa@cea.fr>
 */

#include "TpcDirectLaserReconstruction.h"

#include "TpcSpaceChargeMatrixContainerv1.h"

#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/getClass.h>

#include <g4detectors/PHG4TpcCylinderGeom.h>
#include <g4detectors/PHG4TpcCylinderGeomContainer.h>

#include <trackbase/ActsTrackingGeometry.h>
#include <trackbase_historic/SvtxTrack.h>
#include <trackbase_historic/SvtxTrackMap.h>
#include <trackbase_historic/SvtxTrackState_v1.h>
#include <trackbase/TrkrHitSet.h>
#include <trackbase/TrkrHitSetContainer.h>
#include <trackbase/TrkrHitv2.h>  // for TrkrHit
#include <trackbase/TpcDefs.h>

#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TVector3.h>
#include <TNtuple.h>

#include <cassert>

namespace
{

  //! range adaptor to be able to use range-based for loop
  template<class T> class range_adaptor
  {
    public:
    range_adaptor( const T& range ):m_range(range){}
    inline const typename T::first_type& begin() {return m_range.first;}
    inline const typename T::second_type& end() {return m_range.second;}
    private:
    T m_range;
  };

  //! convenience square method
  template<class T>
    inline constexpr T square( const T& x ) { return x*x; }

  //! get radius from x and y
  template<class T>
    inline constexpr T get_r( const T& x, const T& y ) { return std::sqrt( square(x) + square(y) ); }

  // calculate intersection between line and circle
  std::pair<double, double> line_circle_intersection( const TVector3& p, const TVector3& d, double radius )
  {
    std::pair<double,double> ret = std::make_pair(-1, -1);

    const double A = square(d.x()) + square(d.y());
    const double B = 2*p.x()*d.x() + 2*p.y()*d.y();
    const double C = square(p.x()) + square(p.y()) - square(radius);
    const double delta = square(B)-4*A*C;
    if( delta < 0 ) return ret;

    // we want the first intersection
    const double tup = (-B + std::sqrt(delta))/(2*A);
    const double tdn = (-B-sqrt(delta))/(2*A);
    // std::cout << " tup " << tup << " tdn " << tdn << std::endl;
    ret = std::make_pair(tup, tdn);

    return ret;

    /*
    if( tup > 0 && tup < tdn) 
      return tup;
    if(tdn > 0 && tdn < tup)
       return tdn;

    // no valid extrapolation
    return -1;
    */
  }
    
  /// TVector3 stream
  inline std::ostream& operator << (std::ostream& out, const TVector3& vector )
  {
    out << "( " << vector.x() << ", " << vector.y() << ", " << vector.z() << ")";
    return out;
  }

  /// calculate delta_phi between -pi and pi
  template< class T>
    inline constexpr T delta_phi( const T& phi )
  {
    if( phi >= M_PI ) return phi - 2*M_PI;
    else if( phi < -M_PI ) return phi + 2*M_PI;
    else return phi;
  }

  // phi range
  static constexpr float m_phimin = 0;
  static constexpr float m_phimax = 2.*M_PI;

  // TODO: could try to get the r and z range from TPC geometry
  // r range
  static constexpr float m_rmin = 20;
  static constexpr float m_rmax = 78;

  // z range
  static constexpr float m_zmin = -105.5;
  static constexpr float m_zmax = 105.5;

}

                                                                     
                      
//_____________________________________________________________________
TpcDirectLaserReconstruction::TpcDirectLaserReconstruction( const std::string& name ):
  SubsysReco( name)
  , PHParameterInterface(name)
  , m_matrix_container( new TpcSpaceChargeMatrixContainerv1 )
{
  InitializeParameters();
}

//_____________________________________________________________________
int TpcDirectLaserReconstruction::Init(PHCompositeNode*)
{
  m_total_hits = 0;
  m_matched_hits = 0;
  m_accepted_clusters = 0;
  
  if( m_savehistograms ) create_histograms();
  return Fun4AllReturnCodes::EVENT_OK;
}

//_____________________________________________________________________
int TpcDirectLaserReconstruction::InitRun(PHCompositeNode* )
{
  UpdateParametersWithMacro();
  m_max_dca = get_double_param( "directlaser_max_dca" );
  m_max_drphi = get_double_param( "directlaser_max_drphi" );
  m_max_dz = get_double_param( "directlaser_max_dz" );

  // print
  if( Verbosity() )
  {
    std::cout
      << "TpcDirectLaserReconstruction::InitRun\n"
      << " m_outputfile: " << m_outputfile << "\n"
      << " m_max_dca: " << m_max_dca << "\n"
      << " m_max_drphi: " << m_max_drphi << "\n"
      << " m_max_dz: " << m_max_dz << "\n"
      << std::endl;

    // also identify the matrix container
    m_matrix_container->identify();
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

//_____________________________________________________________________
int TpcDirectLaserReconstruction::process_event(PHCompositeNode* topNode)
{
  // load nodes
  const auto res = load_nodes(topNode);
  if( res != Fun4AllReturnCodes::EVENT_OK ) return res;

  process_tracks();
  return Fun4AllReturnCodes::EVENT_OK;
}

//_____________________________________________________________________
int TpcDirectLaserReconstruction::End(PHCompositeNode* )
{
  // save matrix container in output file
  if( m_matrix_container )
  {
    std::unique_ptr<TFile> outputfile( TFile::Open( m_outputfile.c_str(), "RECREATE" ) );
    outputfile->cd();
    m_matrix_container->Write( "TpcSpaceChargeMatrixContainer" );
  }

  // write evaluation histograms to output
  if( m_savehistograms && m_histogramfile )
  {
    m_histogramfile->cd();
    for(const auto& o:std::initializer_list<TObject*>({ h_dca_layer, h_deltarphi_layer_south,h_deltarphi_layer_north, h_deltaz_layer, h_deltar_r,h_deltheta_delphi,h_deltheta_delphi_1,h_deltheta_delphi_2,h_deltheta_delphi_3,h_deltheta_delphi_4,h_deltheta_delphi_5,h_deltheta_delphi_6,h_deltheta_delphi_7,h_deltheta_delphi_8, h_entries,h_hits,h_clusters,h_origins,h_assoc_hits,h_associated_laser,h_GEMs_hit,h_layers_hit,h_xy,h_xz,h_xy_pca,h_xz_pca,h_dca_path,h_zr,h_zr_pca, h_dz_z }))
    { if( o ) o->Write(); }
    m_histogramfile->Close();
  }

  // print counters
  std::cout << "TpcDirectLaserReconstruction::End - m_total_hits: " << m_total_hits << std::endl;
  std::cout << "TpcDirectLaserReconstruction::End - m_matched_hits: " << m_matched_hits << std::endl;
  std::cout << "TpcDirectLaserReconstruction::End - m_accepted_clusters: " << m_accepted_clusters << std::endl;
  std::cout << "TpcDirectLaserReconstruction::End - fraction cluster/hits: " << m_accepted_clusters/m_total_hits << std::endl;

  return Fun4AllReturnCodes::EVENT_OK;
}

//___________________________________________________________________________
void TpcDirectLaserReconstruction::SetDefaultParameters()
{

  // DCA cut, to decide whether a cluster should be associated to a given laser track or not
  set_default_double_param( "directlaser_max_dca", 20.0 );

  
//   // residual cuts, used to decide if a given cluster is used to fill SC reconstruction matrices
//   set_default_double_param( "directlaser_max_drphi", 0.5 );
//   set_default_double_param( "directlaser_max_dz", 0.5 );

  set_default_double_param( "directlaser_max_drphi", 2. );
  set_default_double_param( "directlaser_max_dz", 2. );
}

//_____________________________________________________________________
void TpcDirectLaserReconstruction::set_grid_dimensions( int phibins, int rbins, int zbins )
{ m_matrix_container->set_grid_dimensions( phibins, rbins, zbins ); }

//_____________________________________________________________________
int TpcDirectLaserReconstruction::load_nodes( PHCompositeNode* topNode )
{
  m_geom_container =  findNode::getClass<PHG4TpcCylinderGeomContainer>(topNode, "CYLINDERCELLGEOM_SVTX");
  assert( m_geom_container );

  // acts geometry
  m_tGeometry = findNode::getClass<ActsGeometry>(topNode, "ActsGeometry");
  assert( m_tGeometry );

  // tracks
  m_track_map = findNode::getClass<SvtxTrackMap>(topNode, "SvtxTrackMap");
  assert(m_track_map);

 // get node containing the digitized hits
  m_hit_map = findNode::getClass<TrkrHitSetContainer>(topNode, "TRKR_HITSET");
  assert(m_hit_map);

  return Fun4AllReturnCodes::EVENT_OK;
}

//_____________________________________________________________________
void TpcDirectLaserReconstruction::create_histograms()
{
  std::cout << "TpcDirectLaserReconstruction::makeHistograms - writing evaluation histograms to: " << m_histogramfilename << std::endl;
  m_histogramfile.reset( new TFile(m_histogramfilename.c_str(), "RECREATE") );
  m_histogramfile->cd();

  // residuals vs layers
  h_dca_layer = new TH2F( "dca_layer", ";radius; DCA (cm)", 78, 0, 78, 500, 0, 20 );
  h_deltarphi_layer_north = new TH2F( "deltarphi_layer_north", ";radius; r.#Delta#phi_{track-cluster} (cm)", 78, 0, 78, 2000, -5, 5 );
  h_deltarphi_layer_south = new TH2F( "deltarphi_layer_south", ";radius; r.#Delta#phi_{track-cluster} (cm)", 78, 0, 78, 2000, -5, 5 );
  h_deltaz_layer = new TH2F( "deltaz_layer", ";radius; #Deltaz_{track-cluster} (cm)", 78, 0, 78, 2000, -20, 20 );
  h_deltar_r = new TH2F("deltar_r",";radius;#Deltar_{track-cluster} (cm)", 78,0,78,2000,-3,3);

  h_xy = new TH2F("h_xy"," x vs y", 320,-80,80,320,-80,80);
  h_xz = new TH2F("h_xz"," x vs z", 320,-80,80,440,-110,110);
  h_xy_pca = new TH2F("h_xy_pca"," x vs y pca", 320,-80,80,320,-80,80);
  h_xz_pca = new TH2F("h_xz_pca"," x vs z pca", 320,-80,80,440,-110,110);
  h_dca_path = new TH2F("h_dca_path"," dca vs pathlength", 440,0,110,100,0,20);
  h_zr = new TH2F("h_zr"," z vs r", 440,-110,110,1000,28,80);
  h_zr->GetXaxis()->SetTitle("z");
  h_zr->GetYaxis()->SetTitle("rad");
  h_zr_pca = new TH2F("h_zr_pca"," z vs r pca", 440,-110,110,1000,28,80);
  h_dz_z = new TH2F("h_dz_z"," dz vs z", 440,-110,110, 1000, -20, 20);
  h_hits = new TNtuple("hits","raw hits","x:y:z");
  h_assoc_hits = new TNtuple("assoc_hits","hits associated with tracks (dca cut)","x:y:z");
  h_clusters = new TNtuple("clusters","associated clusters","x:y:z");
  h_origins = new TNtuple("origins","track origins","x:y:z");
  

  h_deltheta_delphi = new TH2F("deltheta_delphi","#Delta#theta, #Delta#phi for separation b/w TPC volume hits and ALL laser start points", 181, -10.5,180.5,361, -180.5, 180.5);
  h_deltheta_delphi->SetXTitle("#Delta#theta");
  h_deltheta_delphi->SetYTitle("#Delta#phi");

  h_deltheta_delphi_1 = new TH2F("deltheta_delphi_1","#Delta#theta, #Delta#phi for separation b/w TPC volume hits and LASER 0 only", 181, -10.5,180.5,361, -180.5, 180.5);
  h_deltheta_delphi_1->SetXTitle("#Delta#theta");
  h_deltheta_delphi_1->SetYTitle("#Delta#phi");

  h_deltheta_delphi_2 = new TH2F("deltheta_delphi_2","#Delta#theta, #Delta#phi for separation b/w TPC volume hits and LASER 1 only", 181, -10.5,180.5,361, -180.5, 180.5);
  h_deltheta_delphi_2->SetXTitle("#Delta#theta");
  h_deltheta_delphi_2->SetYTitle("#Delta#phi");

  h_deltheta_delphi_3 = new TH2F("deltheta_delphi_3","#Delta#theta, #Delta#phi for separation b/w TPC volume hits and LASER 2 only", 181, -10.5,180.5,361, -180.5, 180.5);
  h_deltheta_delphi_3->SetXTitle("#Delta#theta");
  h_deltheta_delphi_3->SetYTitle("#Delta#phi");

  h_deltheta_delphi_4 = new TH2F("deltheta_delphi_4","#Delta#theta, #Delta#phi for separation b/w TPC volume hits and LASER 3 only", 181, -10.5,180.5,361, -180.5, 180.5);
  h_deltheta_delphi_4->SetXTitle("#Delta#theta");
  h_deltheta_delphi_4->SetYTitle("#Delta#phi");

  h_deltheta_delphi_5 = new TH2F("deltheta_delphi_5","#Delta#theta, #Delta#phi for separation b/w TPC volume hits and LASER 4 only", 181, -10.5,180.5,361, -180.5, 180.5);
  h_deltheta_delphi_5->SetXTitle("#Delta#theta");
  h_deltheta_delphi_5->SetYTitle("#Delta#phi");

  h_deltheta_delphi_6 = new TH2F("deltheta_delphi_6","#Delta#theta, #Delta#phi for separation b/w TPC volume hits and LASER 5 only", 181, -10.5,180.5,361, -180.5, 180.5);
  h_deltheta_delphi_6->SetXTitle("#Delta#theta");
  h_deltheta_delphi_6->SetYTitle("#Delta#phi");

  h_deltheta_delphi_7 = new TH2F("deltheta_delphi_7","#Delta#theta, #Delta#phi for separation b/w TPC volume hits and LASER 6 only", 181, -10.5,180.5,361, -180.5, 180.5);
  h_deltheta_delphi_7->SetXTitle("#Delta#theta");
  h_deltheta_delphi_7->SetYTitle("#Delta#phi");

  h_deltheta_delphi_8 = new TH2F("deltheta_delphi_8","#Delta#theta, #Delta#phi for separation b/w TPC volume hits and LASER 7 only", 181, -10.5,180.5,361, -180.5, 180.5);
  h_deltheta_delphi_8->SetXTitle("#Delta#theta");
  h_deltheta_delphi_8->SetYTitle("#Delta#phi");

  h_GEMs_hit = new TH1F("GEMS_hit","Number of Unique GEM Modules hit for each laser",8,0,8);
  h_layers_hit = new TH1F("layers_hit","Number of Unique Layers hit for each laser",8,8,8);

  char GEM_bin_label[128];
  for(int GEMhistiter = 0; GEMhistiter < 8; GEMhistiter++)
  {  // (pos z) laser 1 {0,60}, laser 2 {60,0}, laser 3 {0,-60}, laser 4 {-60,0}, (neg z) laser 5 {0,60}, laser 2 {60,0}, laser 3 {0,-60}, laser 4 {-60,0}
    sprintf(GEM_bin_label,"laser %i",GEMhistiter+1);
    h_GEMs_hit->GetXaxis()->SetBinLabel(GEMhistiter+1,GEM_bin_label);
    h_layers_hit->GetXaxis()->SetBinLabel(GEMhistiter+1,GEM_bin_label);
  }
  h_GEMs_hit->SetYTitle("Number of Unique GEM Modules Hit");
  h_layers_hit->SetYTitle("Number of Unique Layers Hit");
  
  // <move code here>
  //put code here                                                                                                                                                                                                                    
  //___________________________________________________________________                                                                                                                                                                                                                                                                                                                                                                                                
  //Associated hits to the laser (binning stuff)                                                                                                                                                                                                                                                                                                                                                                                                                        
  float cm = 10.0;  //um=1e-3,//changed to make 'cm' 1.0, for convenience.                                                                                                                                                           
  const int nz = 72;
  double z_rdo = 108 * cm;

  const int r_bins_N = 66;  //51;                                                                                                                                                                                                   
                                                                                                                                                                                                                                     
  double r_bins[r_bins_N + 1] = {217.83,
                                 223.83, 229.83, 235.83, 241.83, 247.83, 253.83, 259.83, 265.83, 271.83, 277.83, 283.83, 289.83, 295.83, 301.83, 306.83,
                                 311.05, 317.92, 323.31, 329.27, 334.63, 340.59, 345.95, 351.91, 357.27, 363.23, 368.59, 374.55, 379.91, 385.87, 391.23, 397.19, 402.49,
                                 411.53, 421.70, 431.90, 442.11, 452.32, 462.52, 472.73, 482.94, 493.14, 503.35, 513.56, 523.76, 533.97, 544.18, 554.39, 564.59, 574.76,
                                 583.67, 594.59, 605.57, 616.54, 627.51, 638.48, 649.45, 660.42, 671.39, 682.36, 693.33, 704.30, 715.27, 726.24, 737.21, 748.18, 759.11};


  const int nphi = 205;

  double phi_bins[nphi + 1] = { 0., 6.3083-2 * M_PI, 6.3401-2 * M_PI, 6.372-2 * M_PI, 6.4039-2 * M_PI, 6.4358-2 * M_PI, 6.4676-2 * M_PI, 6.4995-2 * M_PI, 6.5314-2 * M_PI,
                                0.2618, 0.2937, 0.3256, 0.3574, 0.3893, 0.4212, 0.453, 0.4849, 0.5168, 0.5487, 0.5805, 0.6124, 0.6443, 0.6762, 0.7081,
                                0.7399, 0.7718, 0.7854, 0.8173, 0.8491, 0.881, 0.9129, 0.9448, 0.9767, 1.0085, 1.0404, 1.0723, 1.1041, 1.136, 1.1679,
                                1.1998, 1.2317, 1.2635, 1.2954, 1.309, 1.3409, 1.3727, 1.4046, 1.4365, 1.4684, 1.5002, 1.5321, 1.564, 1.5959, 1.6277,
                                1.6596, 1.6915, 1.7234, 1.7552, 1.7871, 1.819, 1.8326, 1.8645, 1.8963, 1.9282, 1.9601, 1.992, 2.0238, 2.0557, 2.0876,
                                2.1195, 2.1513, 2.1832, 2.2151, 2.247, 2.2788, 2.3107, 2.3426, 2.3562, 2.3881, 2.42, 2.4518, 2.4837, 2.5156, 2.5474,
                                2.5793, 2.6112, 2.6431, 2.6749, 2.7068, 2.7387, 2.7706, 2.8024, 2.8343, 2.8662, 2.8798, 2.9117, 2.9436, 2.9754, 3.0073,
                                3.0392, 3.0711, 3.1029, 3.1348, 3.1667, 3.1986, 3.2304, 3.2623, 3.2942, 3.326, 3.3579, 3.3898, 3.4034, 3.4353, 3.4671,
                                3.499, 3.5309, 3.5628, 3.5946, 3.6265, 3.6584, 3.6903, 3.7221, 3.754, 3.7859, 3.8178, 3.8496, 3.8815, 3.9134, 3.927,
                                3.9589, 3.9907, 4.0226, 4.0545, 4.0864, 4.1182, 4.1501, 4.182, 4.2139, 4.2457, 4.2776, 4.3095, 4.3414, 4.3732, 4.4051,
                                4.437, 4.4506, 4.4825, 4.5143, 4.5462, 4.5781, 4.61, 4.6418, 4.6737, 4.7056, 4.7375, 4.7693, 4.8012, 4.8331, 4.865,
                                4.8968, 4.9287, 4.9606, 4.9742, 5.0061, 5.0379, 5.0698, 5.1017, 5.1336, 5.1654, 5.1973, 5.2292, 5.2611, 5.2929, 5.3248,
                                5.3567, 5.3886, 5.4204, 5.4523, 5.4842, 5.4978, 5.5297, 5.5615, 5.5934, 5.6253, 5.6572, 5.689, 5.7209, 5.7528, 5.7847,
                                5.8165, 5.8484, 5.8803, 5.9122, 5.944, 5.9759, 6.0078, 6.0214, 6.0533, 6.0851, 6.117, 6.1489, 6.1808, 6.2127, 6.2445,
                                6.2764, 2 * M_PI};

  double z_bins[2 * nz + 1];
  for (int z = 0; z <= (2 * nz); z++)
    {
      z_bins[z] = -z_rdo + z_rdo / nz * z;
    }

  h_associated_laser = new TH3F("associated_laser", "Number of associated hits", nphi, phi_bins, r_bins_N, r_bins, 2 * nz, z_bins);

  //_____________________________________________________________________     
 
  // entries vs cell grid
  /* histogram dimension and axis limits must match that of TpcSpaceChargeMatrixContainer */
  if( m_matrix_container )
  {
    int phibins = 0;
    int rbins = 0;
    int zbins = 0;
    m_matrix_container->get_grid_dimensions( phibins, rbins, zbins );
    h_entries = new TH3F( "entries", ";#phi;r (cm);z (cm)", phibins, m_phimin, m_phimax, rbins, m_rmin, m_rmax, zbins, m_zmin, m_zmax );
  }
}                                                                                
                                                                                                                           
  //__________________________________________________________
void TpcDirectLaserReconstruction::process_tracks()
{
  if( !( m_track_map && m_hit_map ) ) return;

  // loop over tracks and process
  for( auto iter = m_track_map->begin(); iter != m_track_map->end(); ++iter )
  { process_track( iter->second ); }
}

//_____________________________________________________________________
void TpcDirectLaserReconstruction::process_track( SvtxTrack* track )
{
  std::multimap<unsigned int, std::pair<float, TVector3>> cluspos_map;
  std::set<unsigned int> layer_bin_set;

  // get track parameters
  const TVector3 origin( track->get_x(), track->get_y(), track->get_z() );
  const TVector3 direction( track->get_px(), track->get_py(), track->get_pz() );

  if(h_origins)
  {
    h_origins->Fill(track->get_x(), track->get_y(), track->get_z());
  }

  const unsigned int trkid = track->get_id();
  if( Verbosity() )
  { 
    std::cout << "----------  processing track " << track->get_id() << std::endl;
    std::cout << "TpcDirectLaserReconstruction::process_track - position: " << origin << " direction: " << direction << std::endl; 
  }
  
  int GEM_Mod_Arr[72] = {0};

  // loop over hits
  TrkrHitSetContainer::ConstRange hitsetrange = m_hit_map->getHitSets(TrkrDefs::TrkrId::tpcId);

  for(auto hitsetitr = hitsetrange.first; hitsetitr != hitsetrange.second; ++hitsetitr)
  {
    const TrkrDefs::hitsetkey& hitsetkey = hitsetitr->first;
    const int side = TpcDefs::getSide(hitsetkey);
    
    auto hitset = hitsetitr->second;
    const unsigned int layer = TrkrDefs::getLayer(hitsetkey);
    const auto layergeom = m_geom_container->GetLayerCellGeom(layer);
    const auto layer_center_radius = layergeom->get_radius();
    
    // maximum drift time.
    /* it is needed to calculate a given hit position from its drift time */
    static constexpr double AdcClockPeriod = 53.0;   // ns 
    const unsigned short NTBins = (unsigned short)layergeom->get_zbins();
    const float tdriftmax =  AdcClockPeriod * NTBins / 2.0;

    // get corresponding hits
    TrkrHitSet::ConstRange hitrangei = hitset->getHits();
    
    for (auto hitr = hitrangei.first; hitr != hitrangei.second; ++hitr)
      {
	++m_total_hits;

	const unsigned short phibin = TpcDefs::getPad(hitr->first);
  const unsigned short zbin = TpcDefs::getTBin(hitr->first);

	const double phi = layergeom->get_phicenter(phibin);
  const double x = layer_center_radius * cos(phi);
  const double y = layer_center_radius * sin(phi);
  
  const double zdriftlength = layergeom->get_zcenter(zbin)*m_tGeometry->get_drift_velocity();
  double z  =  tdriftmax*m_tGeometry->get_drift_velocity() - zdriftlength;
  if(side == 0)  z *= -1;
      
	const TVector3 global(x,y,z);

        if(h_hits)
          {
            h_hits->Fill(x,y,z);
          }
//_________________________________________

//_________________________________________

	float adc = (hitr->second->getAdc()) - m_pedestal; 
	
	// calculate dca
	// origin is track origin, direction is track direction
	const TVector3 oc( global.x()-origin.x(), global.y()-origin.y(), global.z()-origin.z()  );  // vector from track origin to cluster
	auto t = direction.Dot( oc )/square( direction.Mag() );
	auto om = direction*t;     // vector from track origin to PCA
	const auto dca = (oc-om).Mag();

        // relative angle histogram - only fill for hits in the same quadrant as origin
        h_deltheta_delphi->Fill(oc.Theta()*(180./M_PI), oc.Phi()*(180./M_PI)) ;

        if( trkid == 0) //only fill for the first laser !!!
        {
          h_deltheta_delphi_1->Fill(oc.Theta()*(180./M_PI), oc.Phi()*(180./M_PI)) ;
        }
        if( trkid == 1) //only fill for the first laser !!!
        {
          h_deltheta_delphi_2->Fill(oc.Theta()*(180./M_PI), oc.Phi()*(180./M_PI)) ;
        }
        if( trkid == 2) //only fill for the first laser !!!
        {
          h_deltheta_delphi_3->Fill(oc.Theta()*(180./M_PI), oc.Phi()*(180./M_PI)) ;
        }
        if( trkid == 3) //only fill for the first laser !!!
        {
          h_deltheta_delphi_4->Fill(oc.Theta()*(180./M_PI), oc.Phi()*(180./M_PI)) ;
        }
        if( trkid == 4) //only fill for the first laser !!!
        {
          h_deltheta_delphi_5->Fill(oc.Theta()*(180./M_PI), oc.Phi()*(180./M_PI)) ;
        }
        if( trkid == 5) //only fill for the first laser !!!
        {
          h_deltheta_delphi_6->Fill(oc.Theta()*(180./M_PI), oc.Phi()*(180./M_PI)) ;
        }
        if( trkid == 6) //only fill for the first laser !!!
        {
          h_deltheta_delphi_7->Fill(oc.Theta()*(180./M_PI), oc.Phi()*(180./M_PI)) ;
        }
        if( trkid == 7) //only fill for the first laser !!!
        {
          h_deltheta_delphi_8->Fill(oc.Theta()*(180./M_PI), oc.Phi()*(180./M_PI)) ;
        }

	// do not associate if dca is too large
	if( dca > m_max_dca ) continue;

  ++m_matched_hits;

        if(h_assoc_hits){
          h_assoc_hits->Fill( x,y,z);
        }

        //for locating the associated hits   
        float r2 = std::sqrt((x*x) + (y*y));
        float phi2 = phi;
        float z2 = z;
        while( phi2 < m_phimin ) phi2 += 2.*M_PI; 
        while( phi2 >= m_phimax ) phi2 -= 2.*M_PI; 

        const int locateid = Locate(r2 , phi2 , z2); //find where the cluster is 
      
        if( z2 > m_zmin || z2 < m_zmax ) GEM_Mod_Arr[locateid-1]++; // the array ath the cluster location - counts the number of clusters in each array, (IFF its in the volume !!!)

	// bin hits by layer
	const auto cluspos_pair = std::make_pair(adc, global); 	
	cluspos_map.insert(std::make_pair(layer, cluspos_pair));	
	layer_bin_set.insert(layer);
	//_________________________________________                                                                
	//Associated laser stuff 
        if(h_associated_laser){ h_associated_laser->Fill(phi2, r2*10, z2*10);}//If the object exists, fill it

	//std::cout<<"phi:  "<< phi2 << "r:  "<< r2*10 << "z:  "<< z2*10<<std::endl;   

	//std::cout<<"The histogram exists, Luis 03.14!!!!!!!!!!!!"<<std::endl;
	//cout << "r: " << r2 <<::endl;	
       //_________________________________________                                                       
           
      }
  }

  for(int GEMS_iter = 0; GEMS_iter < 72; GEMS_iter++)
  {
    if(GEM_Mod_Arr[GEMS_iter] > 0){  //laser 0
      h_GEMs_hit->Fill(trkid + 0.5);
    }
  }

  
  // We will bring this back here when we fix the clustering
  //int GEM_Mod_Arr[72] = {0};

  // we have all of the residuals for this track added to the map, binned by layer
  // now we calculate the centroid of the clusters in each bin

  for(auto layer : layer_bin_set)
    {
      PHG4TpcCylinderGeom *layergeom = m_geom_container->GetLayerCellGeom(layer);
      const auto layer_center_radius = layergeom->get_radius();
      const auto layer_inner_radius = layer_center_radius - layergeom->get_thickness() / 2.0;
      const auto layer_outer_radius = layer_center_radius + layergeom->get_thickness() / 2.0;

      h_layers_hit->Fill(trkid + 0.5);
          
      // does track pass completely through this layer? If not do not use the hits
      auto tupdn = line_circle_intersection(origin, direction, layer_outer_radius);
      if( tupdn.first <= 0 && tupdn.second <= 0)
	{
	  std::cout << " punt:  layer " << layer << " layer outer radius " << layer_outer_radius << " tup " << tupdn.first << " tdn " << tupdn.second << std::endl; 
	  continue;
	}
      double layer_entry = tupdn.first;
      if(tupdn.second >= 0 && tupdn.second < tupdn.first) layer_entry = tupdn.second; 

      tupdn = line_circle_intersection(origin, direction, layer_inner_radius);
      if( tupdn.first <= 0 && tupdn.second <= 0)
	{
	  std::cout << " punt:  layer " << layer << " layer inner radius " << layer_inner_radius << " tup " << tupdn.first << " tdn " << tupdn.second << std::endl; 
	  continue;
	}
      double layer_exit = tupdn.first;
      if(tupdn.second > 0 && tupdn.second < tupdn.first) layer_exit = tupdn.second;
      if(Verbosity() > 2) std::cout << " layer " << layer << " layer entry " << layer_entry << " layer exit " << layer_exit << std::endl;     

      // calculate track intersection with layer center
      tupdn = line_circle_intersection(origin, direction, layer_center_radius);
      double tup = tupdn.first;
      double tdn = tupdn.second;
      // take the first one if there are two intersections
      double t = tup;
      if(tdn > 0 && tdn < tup) t = tdn;
      if( t < 0 )
	{
	  std::cout << " punt:  layer " << layer << " layer center radius " << layer_center_radius << " t " << t << " tup " << tupdn.first << " tdn " << tupdn.second << std::endl; 
	  continue;
	}

      // get position of layer center on track
      auto om = direction*t;      

      // get vector pointing to the track intersection with layer center
      const auto projection = origin + om;
      
      double zmax = -999.;
      double zmin = 999.;
      
      auto bin_residual = cluspos_map.equal_range(layer);

      TVector3 clus_centroid(0,0,0);
      float wt = 0;
      for( auto mit = bin_residual.first; mit != bin_residual.second; ++mit)
	{
	  float adc =  (float) mit->second.first;
	  TVector3 cluspos =mit->second.second;

	  if(fabs(cluspos.z() - projection.z()) > m_max_zrange) continue;  // to reject clusters from possible second traverse of layer
 
	  if(Verbosity() > 2)
	    {
	      std::cout << "  layer " << layer << " adc " << adc << std::endl;
	      std::cout << "            cluspos " << cluspos.x() << "  " << cluspos.y() << "  " << cluspos.z() 
			<< " clus radius " << sqrt(square(cluspos.x()) + square(cluspos.y())) << std::endl;
	    }

	  clus_centroid += cluspos * adc;
	  wt += adc;

	  if(cluspos.z() < zmin) zmin = cluspos.z();
	  if(cluspos.z() > zmax) zmax = cluspos.z();
	}

      clus_centroid.SetX(clus_centroid.x()/ wt);
      clus_centroid.SetY(clus_centroid.y()/ wt);
      clus_centroid.SetZ(clus_centroid.z()/ wt);

      double zrange =  zmax-zmin;
      if(zrange > m_max_zrange)
	{
	  std::cout << "    exeeded  max zrange:  zrange " << zrange << " max zrange " << m_max_zrange << std::endl; 
	  continue;
	}

      // get the distance of the cluster centroid to the track-layer intersection point

      const TVector3 oc( clus_centroid.x()-origin.x(), clus_centroid.y()-origin.y(), clus_centroid.z()-origin.z()  );  // vector from track origin to cluster

      const auto dca = (oc-om).Mag();

      // path length
      const auto pathlength = om.Mag();

      // Correct cluster z for the track transit time using the pathlength 
      double ns_per_cm = 1e9 / 3e10;
      double dt = pathlength * ns_per_cm;
      double transit_dz = dt * m_tGeometry->get_drift_velocity();
      if(origin.z() > 0)
	clus_centroid.SetZ(clus_centroid.z() + transit_dz);
      else
	clus_centroid.SetZ(clus_centroid.z() - transit_dz);

      if(Verbosity() > 0)
	{
	  std::cout << "  layer " << layer << " radius " << layer_center_radius << " wt " << wt << " t " << t << " dca " << dca 
		    << " pathlength " << pathlength << " transit_dz " << transit_dz << std::endl;
	  std::cout << "      clus_centroid " << clus_centroid.x() << "  " << clus_centroid.y() << "  " << clus_centroid.z() << " zrange " << zrange << std::endl;
	  std::cout << "      projection " << projection.x() << "  " << projection.y() << "  " << projection.z() << " dz " << clus_centroid.z() - projection.z() << std::endl;
	}

      // create relevant state vector and assign to track
      SvtxTrackState_v1 state( pathlength );
      state.set_x( projection.x() );
      state.set_y( projection.y() );
      state.set_z( projection.z() );
      
      state.set_px( direction.x());
      state.set_py( direction.y());
      state.set_pz( direction.z());
      track->insert_state( &state );
      
      // also associate cluster to track
      //track->insert_cluster_key( key );
      
      // cluster r, phi and z
      const auto cluster_r = get_r(clus_centroid.x(), clus_centroid.y());
      const auto cluster_phi = std::atan2(clus_centroid.y(),clus_centroid.x());
      const auto cluster_z = clus_centroid.z();

      // We will bring this back when we fix clustering !!
      /*
      float r2 = cluster_r;
      float phi2 = cluster_phi;
      float z2 = cluster_z;
      while( phi2 < m_phimin ) phi2 += 2.*M_PI; 
      while( phi2 >= m_phimax ) phi2 -= 2.*M_PI; 

      const int locateid = Locate(r2 , phi2 , z2); //find where the cluster is 
      
      if( z2 > m_zmin || z2 < m_zmax ) GEM_Mod_Arr[locateid-1]++; // the array ath the cluster location - counts the number of clusters in each array, (IFF its in the volume !!!)
      */

      // cluster errors
      const auto cluster_rphi_error = 0.015;
      const auto cluster_z_error = 0.075;
      
      // track position
      const auto track_phi = std::atan2( projection.y(), projection.x() );
      const auto track_z = projection.z();
      
      // track angles
      const auto cosphi( std::cos( track_phi ) );
      const auto sinphi( std::sin( track_phi ) );
      const auto track_pphi = -state.get_px()*sinphi + state.get_py()*cosphi;
      const auto track_pr = state.get_px()*cosphi + state.get_py()*sinphi;
      const auto track_pz = state.get_pz();
      const auto talpha = -track_pphi/track_pr;
      const auto tbeta = -track_pz/track_pr;
      
      // sanity check
      if( std::isnan(talpha) )
	{
	  std::cout << "TpcDirectLaserReconstruction::process_track - talpha is nan" << std::endl;
	  continue;
	}
      
      if( std::isnan(tbeta) )
	{
	  std::cout << "TpcDirectLaserReconstruction::process_track - tbeta is nan" << std::endl;
	  continue;
	}
      
      // residuals
      const auto drp = cluster_r*delta_phi( cluster_phi - track_phi );
      const auto dz = cluster_z - track_z;
      
      // sanity checks
      if( std::isnan(drp) )
	{
	  std::cout << "TpcDirectLaserReconstruction::process_track - drp is nan" << std::endl;
	  continue;
	}
      
      if( std::isnan(dz) )
	{
	  std::cout << "TpcDirectLaserReconstruction::process_track - dz is nan" << std::endl;
	  continue;
	}
      
      if(m_savehistograms)
	{
	  const float r = get_r( projection.x(), projection.y() );
	  const float dr = cluster_r - r;
	  if(h_dca_layer) h_dca_layer->Fill(r, dca);
	  if(h_deltarphi_layer_south && clus_centroid.z() < 0) h_deltarphi_layer_south->Fill(r, drp);
	  if(h_deltarphi_layer_north && clus_centroid.z() > 0) h_deltarphi_layer_north->Fill(r, drp);
	  if(h_deltaz_layer) h_deltaz_layer->Fill(r, dz);
	  if(h_deltar_r) h_deltar_r->Fill(r, dr);
	  if(h_entries)
	    {
	      auto phi = cluster_phi;
	      while( phi < m_phimin ) phi += 2.*M_PI;
	      while( phi >= m_phimax ) phi -= 2.*M_PI;
	      h_entries->Fill( phi, cluster_r, cluster_z );
	    }
	  if(h_xy) h_xy->Fill(clus_centroid.x(), clus_centroid.y());
	  if(h_xz) h_xz->Fill(clus_centroid.x(), clus_centroid.z());
	  if(h_xy_pca) h_xy_pca->Fill(projection.x(), projection.y());
	  if(h_xz_pca) h_xz_pca->Fill(projection.x(), projection.z());
	  if(h_dca_path) h_dca_path->Fill(pathlength,dca);
	  if(h_zr) h_zr->Fill(clus_centroid.z(), cluster_r);
	  if(h_zr_pca) h_zr_pca->Fill(projection.z(), r);
	  if(h_dz_z)h_dz_z->Fill(projection.z(), clus_centroid.z()- projection.z());
          if(h_clusters)h_clusters->Fill(clus_centroid.x(),clus_centroid.y(),clus_centroid.z());
	}
      
      //       // check against limits
      //       if( std::abs( drp ) > m_max_drphi ) continue;
      //       if( std::abs( dz ) > m_max_dz ) continue;
      
      // residual errors squared
      const auto erp = square(cluster_rphi_error);
      const auto ez = square(cluster_z_error);
      
      // sanity check
      if( std::isnan( erp ) )
	{
	  std::cout << "TpcDirectLaserReconstruction::process_track - erp is nan" << std::endl;
	  continue;
	}
      
      if( std::isnan( ez ) )
	{
	  std::cout << "TpcDirectLaserReconstruction::process_track - ez is nan" << std::endl;
	  continue;
	}
      
      // get cell
      const auto i = get_cell_index( clus_centroid );
      if( i < 0 )
	{
	  if( Verbosity() )
	    {
	      std::cout << "TpcDirectLaserReconstruction::process_track - invalid cell index"
			<< " r: " << cluster_r
			<< " phi: " << cluster_phi
			<< " z: " << cluster_z
			<< std::endl;
	    }
	  continue;
	}
      
      // update matrices
      // see https://indico.bnl.gov/event/7440/contributions/43328/attachments/31334/49446/talk.pdf for details
      m_matrix_container->add_to_lhs(i, 0, 0, 1./erp );
      m_matrix_container->add_to_lhs(i, 0, 1, 0 );
      m_matrix_container->add_to_lhs(i, 0, 2, talpha/erp );
      
      m_matrix_container->add_to_lhs(i, 1, 0, 0 );
      m_matrix_container->add_to_lhs(i, 1, 1, 1./ez );
      m_matrix_container->add_to_lhs(i, 1, 2, tbeta/ez );
      
      m_matrix_container->add_to_lhs(i, 2, 0, talpha/erp );
      m_matrix_container->add_to_lhs(i, 2, 1, tbeta/ez );
      m_matrix_container->add_to_lhs(i, 2, 2, square(talpha)/erp + square(tbeta)/ez );
      
      m_matrix_container->add_to_rhs(i, 0, drp/erp );
      m_matrix_container->add_to_rhs(i, 1, dz/ez );
      m_matrix_container->add_to_rhs(i, 2, talpha*drp/erp + tbeta*dz/ez );
      
      // update entries in cell
      m_matrix_container->add_to_entries(i);
      
      // increment number of accepted clusters
      ++m_accepted_clusters;      
    }

// We will eventually bring this back when we fix clustering
/*
  for(int GEMS_iter = 0; GEMS_iter < 72; GEMS_iter++)
  {
    if(GEM_Mod_Arr[GEMS_iter] > 0){  //laser 0
      h_GEMs_hit->Fill(trkid + 0.5);
    }
  }
*/

}


//_____________________________________________________________________
int TpcDirectLaserReconstruction::get_cell_index( const TVector3& global ) const
{
  // get grid dimensions from matrix container
  int phibins = 0;
  int rbins = 0;
  int zbins = 0;
  m_matrix_container->get_grid_dimensions( phibins, rbins, zbins );

  // phi
  // bound check
  float phi = std::atan2( global.y(), global.x() );
  while( phi < m_phimin ) phi += 2.*M_PI;
  while( phi >= m_phimax ) phi -= 2.*M_PI;
  int iphi = phibins*(phi-m_phimin)/(m_phimax-m_phimin);

  // radius
  const float r = get_r( global.x(), global.y() );
  if( r < m_rmin || r >= m_rmax ) return -1;
  int ir = rbins*(r-m_rmin)/(m_rmax-m_rmin);

  // z
  const float z = global.z();
  if( z < m_zmin || z >= m_zmax ) return -1;
  int iz = zbins*(z-m_zmin)/(m_zmax-m_zmin);

  return m_matrix_container->get_cell_index( iphi, ir, iz );
}
//_____________________________________________________________________
int TpcDirectLaserReconstruction::Locate(float r , float phi , float z)
{

  ////////////////////////////////////////////////////////////////////////
  //          North Label Conventions                                  //
  ///////////////////////////////////////////////////////////////////////
  //   1 - 00_R1   16 - 05_R1   31 - 10_R1    
  //   2 - 00_R2   17 - 05_R2   32 - 10_R2   
  //   3 - 00_R3   18 - 05_R3   33 - 10_R3 
  //   4 - 01_R1   19 - 06_R1   34 - 11_R1
  //   5 - 01_R2   20 - 06_R2   35 - 11_R2
  //   6 - 01_R3   21 - 06_R3   36 - 11_R3
  //   7 - 02_R1   22 - 07_R1
  //   8 - 02_R2   23 - 07_R2
  //   9 - 02_R3   24 - 07_R3
  //  10 - 03_R1   25 - 08_R1
  //  11 - 03_R2   26 - 08_R2
  //  12 - 03_R3   27 - 08_R3
  //  13 - 04_R1   28 - 09_R1
  //  14 - 04_R2   29 - 09_R2
  //  15 - 04_R3   30 - 09_R3

  //////////////////////////////////////////////////////////////////////// 
  //          South Label Conventions                                   //  
  ///////////////////////////////////////////////////////////////////////  
  //  37 - 12_R1   52 - 17_R1   67 - 22_R1    
  //  38 - 12_R2   53 - 17_R2   68 - 22_R2   
  //  39 - 12_R3   54 - 17_R3   69 - 22_R3 
  //  40 - 13_R1   55 - 18_R1   70 - 23_R1
  //  41 - 13_R2   56 - 18_R2   71 - 23_R2
  //  42 - 13_R3   57 - 18_R3   72 - 23_R3
  //  43 - 14_R1   58 - 19_R1
  //  44 - 14_R2   59 - 19_R2
  //  45 - 14_R3   60 - 19_R3
  //  46 - 15_R1   61 - 20_R1
  //  47 - 15_R2   62 - 20_R2
  //  48 - 15_R3   63 - 20_R3
  //  49 - 16_R1   64 - 21_R1
  //  50 - 16_R2   65 - 21_R2
  //  51 - 16_R3   66 - 21_R3


  int GEM_Mod_ID; //integer from 1 to 72

  const float Angle_Bins[13] = {23*M_PI/12,1*M_PI/12,3*M_PI/12,5*M_PI/12,7*M_PI/12,9*M_PI/12,11*M_PI/12,13*M_PI/12,15*M_PI/12,17*M_PI/12,19*M_PI/12,21*M_PI/12,23*M_PI/12}; //12 angle bins on each side

  const float r_bins[4] = {30,46,62,78}; //4 r bins 

  int angle_id = 0; // default to 0
  int r_id = 0; //default to 0
  int side_id = 0; //default to 0 (north side)

  for(int r_iter = 0; r_iter < 3; r_iter++ ){
    if( r > r_bins[r_iter] && r < r_bins[r_iter+1] )
    {
      r_id = r_iter+1;
      break; //break out of the for loop if you found where the hit is in r
    }
  }

  for(int ang_iter = 0; ang_iter < 12; ang_iter++ ){
    if( phi > Angle_Bins[ang_iter] && phi < Angle_Bins[ang_iter+1] ) 
    {
      angle_id = ang_iter;
      break; //break out the for loop if you found where the hit is in phi
    }
  }
  if(z < 0)side_id=1; //(south side)

  return GEM_Mod_ID = (36*side_id) + (3*angle_id) + r_id;
}
