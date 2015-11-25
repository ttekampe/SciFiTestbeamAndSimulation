#ifndef MCFTDEPOSITPOSITIONTOOL_H
#define MCFTDEPOSITPOSITIONTOOL_H 1

// Include files
// from Gaudi
#include "GaudiAlg/GaudiHistoTool.h"
// FTDet
#include "FTDet/DeFTDetector.h"

#include "IMCFTDepositPositionTool.h"            // Interface

#include "IMCFTEnergyInFibreTool.h"


/** @class MCFTDepositPositionTool MCFTDepositPositionTool.h
 *
 *  Tool that distributes the deposited energy among the SiPM channels
 *  and handles the light sharing
 *
 *  @author Tobias Tekampe
 *  @date   2015-10-23
 */


class MCFTDepositPositionTool : public GaudiHistoTool, virtual public IMCFTDepositPositionTool {

public:
  /// Standard constructor
  MCFTDepositPositionTool( const std::string& type,
                           const std::string& name,
                           const IInterface* parent);

  virtual ~MCFTDepositPositionTool( ); ///< Destructor

  /// Initialize the transmission map
  StatusCode initializeTool();

  /// Calculate fired SiPM channels and energy taking into account light sharing
  void firedChannels(const LHCb::MCHit* ftHit, const DeFTFibreMat& fibMat, FTDoublePairs& channels);

protected:

private:
  // Tool that calculates the distance the particle travelled within each fibre core that was hit
  IMCFTEnergyInFibreTool* m_energyInFibreTool;
  std::string m_energyInFibreToolName;


};
#endif // MCFTDEPOSITPOSITIONTOOL_H
