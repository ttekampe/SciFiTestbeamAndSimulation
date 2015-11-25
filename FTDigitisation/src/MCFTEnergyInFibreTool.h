#ifndef MCFTENERGYINFIBRETOOL_H 
#define MCFTENERGYINFIBRETOOL_H 1

// Include files
// from Gaudi
#include "GaudiAlg/GaudiHistoTool.h"

// from STL
#include <iostream>
#include <algorithm>
#include <string>
#include <vector>
#include <cmath>

// from boost
#include <boost/lambda/lambda.hpp>

// from ROOT
#include "Math/SVector.h"

#include "GaudiKernel/Point3DTypes.h"


#include "IMCFTEnergyInFibreTool.h"            // Interface

/** @class MCFTEnergyInFibreTool MCFTEnergyInFibreTool.h
 *  
 *  This tool calculates the fraction of the total path length
 *  of a track in the core of each fibre that was hit along with
 *  the position of the fibre.
 *
 *  @author Tobias Tekampe, based on Julian Wishahis MCSciFiEnergyInFibreCoreCalc
 *  @date   2015-10-27
 */


typedef ROOT::Math::SVector<double, 2> Vec2;
typedef ROOT::Math::SVector<double, 3> Vec3;

class MCFTEnergyInFibreTool : public GaudiHistoTool, virtual public IMCFTEnergyInFibreTool {

public: 

  MCFTEnergyInFibreTool(const std::string& type, 
                        const std::string& name,
                        const IInterface* parent); 

  ~MCFTEnergyInFibreTool();

  StatusCode initializeTool();

  std::vector<std::pair<Gaudi::XYZPoint,double>> EffectivePathFracInCores(
            const double x_in,
            const double y_in,
            const double z_in,
            const double dx_over_dz,
            const double dy_over_dz,
            const double zpath,
            const double outermat_thickness) const;

protected:

private:

  double FibreMatThickness(unsigned int nlayers) {
    return (nlayers-1)*m_delta_z+m_fibre_diameter;
  }

  double m_fibre_diameter;
  double m_fibre_core_diameter;
  double m_fibre_pitch;
  double m_delta_x;
  double m_delta_z;
  unsigned int m_nlayers;
  double m_inner_mat_thickness;
  
};

double DistanceToSquared( const Vec2& v, const Vec2& p );
double DistanceTo( const Vec2& v, const Vec2& p );
double MinDistanceFromLineSegmentToPoint( const Vec2& v, const Vec2& w, const Vec2& p );

#endif // MCFTENERGYINFIBRETOOL_H
