#ifndef IMCFTDEPOSITPOSITIONTOOL_H
#define IMCFTDEPOSITPOSITIONTOOL_H 1

// Include files
// from STL
#include <string>

// from Gaudi
#include "GaudiKernel/IAlgTool.h"

static const InterfaceID IID_IMCFTDepositPositionTool ( "IMCFTDepositPositionTool", 1, 0 );

/** @class IMCFTDepositPositionTool IMCFTDepositPositionTool.h
 *
 *  Interface for the tool that distributes the deposited energy among the SiPM channels
 *  and handles the light sharing
 *  
 *  @author Tobias Tekampe
 *  @date   2015-10-23
 */

namespace LHCb{
  class MCHit;
}




class IMCFTDepositPositionTool : virtual public IAlgTool {

public:

  typedef std::pair<LHCb::FTChannelID, double> FTDoublePair;
  typedef std::vector< FTDoublePair > FTDoublePairs;

  // Return the interface ID
  static const InterfaceID& interfaceID() { return IID_IMCFTDepositPositionTool; }

  virtual void firedChannels(const LHCb::MCHit* ftHit, const DeFTFibreMat& fibMat, FTDoublePairs& channels) = 0;
  virtual StatusCode initializeTool() = 0;

protected:

private:

};
#endif // IMCFTDEPOSITPOSITIONTOOL_H
