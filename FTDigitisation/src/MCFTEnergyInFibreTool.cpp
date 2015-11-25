// Include files 

// from Gaudi
#include "GaudiKernel/ToolFactory.h" 
// local
#include "MCFTEnergyInFibreTool.h"

// from Event
#include "Event/MCHit.h"


//-----------------------------------------------------------------------------
// Implementation file for class : MCFTEnergyInFibreTool
//
// 2015-10-27 : Tobias Tekampe, based on Julian Wishahis MCSciFiEnergyInFibreCoreCalc
//-----------------------------------------------------------------------------

// Declaration of the Tool Factory
DECLARE_TOOL_FACTORY( MCFTEnergyInFibreTool )


//=============================================================================
// Standard constructor, initializes variables
//=============================================================================
MCFTEnergyInFibreTool::MCFTEnergyInFibreTool( const std::string& type,
                                              const std::string& name,
                                              const IInterface* parent )
  : GaudiHistoTool ( type, name , parent )
{
  declareInterface<IMCFTEnergyInFibreTool>(this);


  declareProperty( "FibreDiameter", m_fibre_diameter = 0.250 * Gaudi::Units::mm, 
                   "Diameter of the fibre in the fibre mats");

  declareProperty( "FibreCoreDiameter", m_fibre_core_diameter = 0.220 * Gaudi::Units::mm, 
                   "Diameter of the fibre core");

  declareProperty( "FibrePitch", m_fibre_pitch = 0.270 * Gaudi::Units::mm, 
                   "Distance between two neighboring fibre centers in x direction");  

  declareProperty( "DeltaZ", m_delta_z = 0.205 * Gaudi::Units::mm, 
                   "Pitch in z-direction");  

  declareProperty( "Nlayers", m_nlayers = 6, 
                   "Number of layers");

  m_inner_mat_thickness = FibreMatThickness(m_nlayers);
  


}
//=============================================================================
// Destructor
//=============================================================================
MCFTEnergyInFibreTool::~MCFTEnergyInFibreTool() {

}
//=============================================================================
// Calculate the EnergyInFibre
//=============================================================================

StatusCode MCFTEnergyInFibreTool::initializeTool(){
  
  return StatusCode::SUCCESS;
}


std::vector<std::pair<Gaudi::XYZPoint,double>> MCFTEnergyInFibreTool::EffectivePathFracInCores(
        const double x_in,
        const double y_in,
        const double z_in,
        const double dx_over_dz,
        const double dy_over_dz,
        const double zpath,
        const double outermat_thickness) const {
  
  // set in and out vector in 2d plane 
  Vec2 vec_in(x_in, z_in);
  Vec2 vec_out(x_in+dx_over_dz*zpath, z_in+zpath); 

  // flight and slopw direction
  int zpath_sign = (zpath > 0) ? 1: ((zpath < 0) ? -1 : 0); 
  int dx_over_dz_sign = (dx_over_dz > 0) ? 1 : ((dx_over_dz < 0) ? -1 : 0);

  // total path length (not entirely correct). 
  // TODO: What about paths that are not fully contained in a mat?
  double path_x = dx_over_dz*zpath;
  double path_y = dy_over_dz*zpath;
  double path_z = zpath;
  double total_path = sqrt(path_x*path_x+path_y*path_y+path_z*path_z);

  double init_glue_thickness = (outermat_thickness - m_inner_mat_thickness)/2.;
  
  // shift in coordinate systems to have 0 at mid of mat and not at beginning of mat
  double mat_mid_shift_z = -1.*outermat_thickness/2.;

  // x position of hit path at inner mat boundary
  double x_in_at_inner = x_in + zpath_sign * dx_over_dz * init_glue_thickness;

  // initial fibre (right) of hit position for positive (negative) x_in_at_inner
  double fibre_x = static_cast<int>(x_in_at_inner/m_fibre_pitch) * m_fibre_pitch;
  fibre_x -= zpath_sign*dx_over_dz_sign*m_fibre_pitch; // go one right/ left
  double fibre_z = zpath_sign*(mat_mid_shift_z + (init_glue_thickness + m_fibre_diameter/2.));

  // intialise vector of the fibre center
  Vec2 vec_fibre_center(fibre_x,fibre_z);
  // std::cout << " First fibre at (" << vec_fibre_center(0) << ", " 
  //    << vec_fibre_center(1) << ")." << std::endl;


  // max number of fibres to consider
  unsigned int nfibres = static_cast<int>(std::abs(dx_over_dz)*m_fibre_diameter/m_fibre_pitch)+4;


  // helper variables (get rid of unnecessary initialisations in loop)
  double fibre_core_radius = m_fibre_core_diameter/2.;
  double perp_dist = 0.; 
  double dist_in_core = 0.;
  
  std::vector<std::pair<Gaudi::XYZPoint,double>> fibrepos_and_pathfracs;
  // loop over each fibre in each layer
  for (unsigned int layer_i = 0; layer_i < m_nlayers; ++layer_i) {
    for (unsigned int fibre_i = 0; fibre_i < nfibres; ++fibre_i) {

      // calculate distance to fibre
      perp_dist = MinDistanceFromLineSegmentToPoint(vec_in, vec_out, vec_fibre_center);
      
      dist_in_core = 0.;
      // distance in fibres only if distance is smaller than the core radius
      if (perp_dist < fibre_core_radius) {
        dist_in_core = 2.*sqrt(fibre_core_radius*fibre_core_radius - perp_dist*perp_dist);
        fibrepos_and_pathfracs.push_back({
          Gaudi::XYZPoint(vec_fibre_center(0),y_in+(vec_fibre_center(1)-z_in)*dy_over_dz, vec_fibre_center(1)),
          dist_in_core/total_path
        });
      }
      // std::cout << dist_in_core << "." << std::endl;
      // set fibre x position for next iteration
      vec_fibre_center(0) += (m_fibre_pitch*dx_over_dz_sign*zpath_sign);
    }
    // get x position of hit at next fibre mat layer
    x_in_at_inner += (dx_over_dz*zpath_sign*m_delta_z);

    vec_fibre_center(0) = (static_cast<int>(x_in_at_inner/m_fibre_pitch)*m_fibre_pitch);
    // move by half a pitch when switching to an uneven layer in the next step
    if (layer_i%2 == 0)  vec_fibre_center(0) -= zpath_sign*dx_over_dz_sign*m_fibre_pitch/2.;
    vec_fibre_center(1) += zpath_sign*m_delta_z;
  }

  return fibrepos_and_pathfracs;

}
  


double DistanceToSquared( const Vec2& v, const Vec2& p )
{
  const double dX = p(0) - v(0); //p(0) ^= p_x
  const double dY = p(1) - v(1);

  return dX * dX + dY * dY;
}

double DistanceTo( const Vec2& v, const Vec2& p )
{
  return sqrt( DistanceToSquared(v, p) );
}

double MinDistanceFromLineSegmentToPoint( const Vec2& v, const Vec2& w, const Vec2& p )
{
  const double distSq = DistanceToSquared(v, w);                                // i.e. |w-v|^2 ... avoid a sqrt
  if ( distSq == 0.0 ) return DistanceTo(v, p);                                 // v == w case


  // consider the line extending the segment, parameterized as v + t (w - v)
  // we find projection of point p onto the line
  // it falls where t = [(p-v) . (w-v)] / |w-v|^2

  const double t = ROOT::Math::Dot( p-v, w-v ) / distSq;
  if ( t < 0.0 ) return DistanceTo(v, p);                                       // beyond the v end of the segment
  else if ( t > 1.0 ) return DistanceTo(w, p);                                  // beyond the w end of the segment

  // projection falls on the segment
  const Vec2 projection = v + ( ( w - v ) * t );

  return DistanceTo(p, projection);
}
