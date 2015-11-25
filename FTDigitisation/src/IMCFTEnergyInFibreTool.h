#ifndef IMCFTENERGYINFIBRETOOL_H 
#define IMCFTENERGYINFIBRETOOL_H 1

// Include files
// from STL
#include <string>

// from Gaudi
#include "GaudiKernel/IAlgTool.h"

static const InterfaceID IID_IMCFTEnergyInFibreTool ( "IMCFTEnergyInFibreTool", 1, 0 );

/** @class IMCFTEnergyInFibreTool IMCFTEnergyInFibreTool.h
 *  
 *  Interface for the tool that calculates the coordinates and fractional path length
 *  in the core of each fibre
 *
 *  @author Tobais Tekampe
 *  @date   2015-10-27
 */

namespace LHCb{
  class MCHit;
}


class IMCFTEnergyInFibreTool : virtual public IAlgTool {
public: 

  // Return the interface ID
  static const InterfaceID& interfaceID() { return IID_IMCFTEnergyInFibreTool; }

  virtual std::vector<std::pair<Gaudi::XYZPoint,double>> EffectivePathFracInCores(
  			const double x_in, 
  			const double y_in, 
  			const double z_in, 
  			const double dx_over_dz, 
  			const double dy_over_dz,
        const double zpath,
  			const double outermat_thickness) const = 0;

  virtual StatusCode initializeTool() = 0;
  
protected:

private:

};
#endif // IMCFTENERGYINFIBRETOOL_H
