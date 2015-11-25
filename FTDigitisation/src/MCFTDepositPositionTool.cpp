// Include files

// from Gaudi
#include "GaudiKernel/ToolFactory.h"
// local
#include "MCFTDepositPositionTool.h"

// from Event
#include "Event/MCHit.h"

//from STL
#include <map>
#include <vector>
#include <utility>


//-----------------------------------------------------------------------------
// Implementation file for class : MCFTDepositPositionTool
//
// 2015-10-23 : Tobias Tekampe
//-----------------------------------------------------------------------------

// Declaration of the Tool Factory
DECLARE_TOOL_FACTORY( MCFTDepositPositionTool )


//=============================================================================
// Standard constructor, initializes variables
//=============================================================================
MCFTDepositPositionTool::MCFTDepositPositionTool( const std::string& type,
                                                  const std::string& name,
                                                  const IInterface* parent )
  : GaudiHistoTool ( type, name , parent )
{
  declareInterface<IMCFTDepositPositionTool>(this);

  declareProperty( "EnergyInFibreToolName", m_energyInFibreToolName = "MCFTEnergyInFibreTool");

  

}
//=============================================================================
// Destructor
//=============================================================================

struct lightShareInfo{
  double fractionLeft;
  double fractionSelf;
  double fractionRight;
};

MCFTDepositPositionTool::~MCFTDepositPositionTool() {

}
//=========================================================================
// Initialise needed tools
//=========================================================================
StatusCode MCFTDepositPositionTool::initializeTool(){


  m_energyInFibreTool = tool<IMCFTEnergyInFibreTool>( m_energyInFibreToolName, this );
  if( m_energyInFibreTool == nullptr ){
    error() << "Could not find: " << m_energyInFibreToolName << endmsg;
  }


  return StatusCode::SUCCESS;


}

//=============================================================================
// Calculate fired channels
//=============================================================================
void MCFTDepositPositionTool::firedChannels(const LHCb::MCHit* ftHit,
                                            const DeFTFibreMat& fibMat, 
                                            FTDoublePairs& channels){

    //transform coordinates of the MCHit to the local FibreMat coordinates
  Gaudi::XYZPoint enPLocal = fibMat.geometry()->toLocal( ftHit->entry() );
  Gaudi::XYZPoint exPLocal = fibMat.geometry()->toLocal( ftHit->exit()  );

  //Calculate the positions of the fibres that were hit along with the fractional path lengths
  //within these


  double dz = exPLocal.Z() - enPLocal.Z(); 
  double dx_over_dz = ( exPLocal.X() - enPLocal.X() ) / dz;
  double dy_over_dz = ( exPLocal.Y() - enPLocal.Y() ) / dz;
  double fibMatThickness = fibMat.layerMaxZ() - fibMat.layerMinZ();

  //std::cout << "fibMatThickness is " << fibMatThickness << "\n";

  

  std::vector<std::pair<Gaudi::XYZPoint,double>> fibreCentersAndFracPathLengths 
  = m_energyInFibreTool->EffectivePathFracInCores(enPLocal.X(), 
                                                  enPLocal.Y(), 
                                                  enPLocal.Z(), 
                                                  dx_over_dz, 
                                                  dy_over_dz,
                                                  dz,
                                                  fibMatThickness);

  //std::cout << "Frac path length:\n";
  double totFracPath = 0.;
  for(const auto& fibCoordAndFracPath : fibreCentersAndFracPathLengths){
    totFracPath += fibCoordAndFracPath.second;
  //  std::cout << fibCoordAndFracPath.second << " in (" 
  //            << fibCoordAndFracPath.first.X() << ", "
  //            << fibCoordAndFracPath.first.Y() << ", "
  //            << fibCoordAndFracPath.first.Z() << ")\n";
  }

  plot(totFracPath, "fracPathLength", "Fraction of pathlength in active material ; Fraction;Number of hits" , -0.01 , 1.01 );

  //if(totFracPath == 0){
  //  std::cout << "Frac path is 0\n"
  //            << "In Local coordinates:\n"
  //            << "\tEntry point: (" << enPLocal.X() << ", " << enPLocal.Y() << ", " << enPLocal.Z() << ")\n"
  //            << "\tExit point: (" << exPLocal.X() << ", " << exPLocal.Y() << ", " << exPLocal.Z() << ")\n"
  //            << "In MCHit coordinates:\n"
  //            << "\tEntry point: (" << ftHit->entry().X() << ", " << ftHit->entry().Y() << ", " << ftHit->entry().Z() << ")\n"
  //            << "\tExit point: (" << ftHit->exit().X() << ", " << ftHit->exit().Y() << ", " << ftHit->exit().Z() << ")\n";
  //}

  unsigned int sipmID;
  unsigned int cellID; 
  double fracDistCellCenter;
  double sipmREdgeX;

  std::map<unsigned int, double> chanIdsandEnergieFracs;

  std::map<unsigned int, lightShareInfo> chanIDsAndSharedLight;

  //loop over all fibres that were hit

  LHCb::FTChannelID currentChannel;
  for(const auto& fibCenterAndFracPathLength : fibreCentersAndFracPathLengths){

    fibMat.cellIDCoordinates(fibCenterAndFracPathLength.first.X(), 
                             fibMat.quarter(),
                             sipmID,
                             cellID, 
                             fracDistCellCenter,
                             sipmREdgeX
                             );

    //distribute deposit weights among the channels

    currentChannel = fibMat.createChannel(fibMat.layer(), fibMat.module(), fibMat.FibreMatID()%10, sipmID, cellID);

    chanIdsandEnergieFracs[currentChannel.channelID()] += fibCenterAndFracPathLength.second;

    //std::cout << "adding channel with id " << currentChannel.channelID() << " and weight " << fibCenterAndFracPathLength.second << "\n";

    //strore information on lightsharing

    if (fracDistCellCenter > 0) {
      chanIDsAndSharedLight[currentChannel.channelID()] = {
        //currentChannel, //donnator,
        0.68 * fracDistCellCenter + 0.16, //fractionLeft
        -0.36 * fracDistCellCenter + 0.68, //fractionSelf,
        -0.32 * fracDistCellCenter + 0.16 //fractionRight
      };
    }
    else {
      chanIDsAndSharedLight[currentChannel.channelID()] = {
        //currentChannel, //donnator,
        0.32 * fracDistCellCenter + 0.16, //fractionLeft,
        0.36 * fracDistCellCenter + 0.68, //fractionSelf,
        -0.68 * fracDistCellCenter + 0.16 //fraction right
      };

    }
  }
  
  std::map<unsigned int, double> chanIdsandEnergieFracsWithLightShare;

  unsigned int currentChanID, chan2RightID, chan2LeftID;

  // apply light sharing
  for(const auto& chanIdandSharedLight : chanIDsAndSharedLight){
    currentChanID = chanIdandSharedLight.first;
    chan2LeftID = fibMat.nextChannelLeft(LHCb::FTChannelID(currentChanID)).channelID();
    chan2RightID = fibMat.nextChannelRight(LHCb::FTChannelID(currentChanID)).channelID();

    //std::cout << "Adding energy to channels " << chan2LeftID << ", " << currentChanID << " and " << chan2RightID << "\n";

    chanIdsandEnergieFracsWithLightShare[chan2LeftID] += chanIdsandEnergieFracs[currentChanID] * chanIdandSharedLight.second.fractionLeft;
    chanIdsandEnergieFracsWithLightShare[currentChanID] += chanIdsandEnergieFracs[currentChanID] * chanIdandSharedLight.second.fractionSelf;
    chanIdsandEnergieFracsWithLightShare[chan2RightID] += chanIdsandEnergieFracs[currentChanID] * chanIdandSharedLight.second.fractionRight;
  }


  //transform chanIdsandEnergieFracs into vector of FTDoublePairs and multiply weight by total energy deposit
  channels = FTDoublePairs(chanIdsandEnergieFracsWithLightShare.size());
  unsigned int i=0;
  //double totalEfrac = 0.;
  for(const auto& chanIdAndEnergyFrac : chanIdsandEnergieFracsWithLightShare){
    //totalEfrac += chanIdAndEnergyFrac.second;
    channels[i] = std::make_pair(
      LHCb::FTChannelID(chanIdAndEnergyFrac.first),
      chanIdAndEnergyFrac.second * ftHit->energy()
    );
    if(channels[i].second < 0.){
      std::cerr << "Warning: negative energy in SiPM channel!\n";
    }
    ++i;
  }
  // std::cout << "Deposited energy that is in fibre core: " << totalEfrac << "\n"
  //           << "Total fraction of parh in active material: " <<totFracPath << "\n";
   
}

