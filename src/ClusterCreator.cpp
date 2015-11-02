//from here
#include "ClusterCreator.h"
#include "Cluster.h"


//from std
#include <iostream>
#include <utility>


ClusterCreator::~ClusterCreator(){
  for(auto& clPtr : clusters){
    delete clPtr;
  }
}

std::vector<Cluster*> ClusterCreator::FindClustersInEventMax(
                            const Event& event
                            ,const double neighbour_threshold
                            ,const double seed_threshold
                            ,const double sum_threshold
                            )
{
  std::vector<Cluster*> foundClusters;
  Cluster* currentCluster = nullptr;

  for(Event::const_iterator chan = event.begin(); chan != event.end(); ++chan){
    if( chan->AdcValue > neighbour_threshold){ //look for channels exceeding cluster threshold
      if(!currentCluster){
        currentCluster = new Cluster();
        currentCluster->AddChannel(chan->Uplink, chan->ChannelNumber, chan->AdcValue);
      }
      else{
        currentCluster->AddChannel(chan->Uplink, chan->ChannelNumber, chan->AdcValue);
      }
    }
    if( chan->AdcValue <= neighbour_threshold && currentCluster){
      if(currentCluster->GetMaximumAdcValue() > seed_threshold && currentCluster->GetSumOfAdcValues() > sum_threshold){
           clusters.push_back(currentCluster);
           foundClusters.push_back(currentCluster);
       }
      currentCluster = nullptr;
    }
  }
  return foundClusters;
}


std::vector<Cluster*> ClusterCreator::FindClustersInEventBoole(const Event& event,
                              const double neighbourThreshold,
                              const double seedThreshold,
                              const double sumThreshold,
                              const int maxClusterSize,
                              bool debug)
{
  std::vector<Cluster*> foundClusters;
  std::vector<Channel>::const_iterator lastStopDigitIter = event.begin(); // end digit of last cluster, to prevent overlap

  // Since Digit Container is sorted wrt channelID, clusters are defined searching for bumps of ADC Count
  std::vector<Channel>::const_iterator seedDigitIter = event.begin();

  while(seedDigitIter != event.end()){
      // loop over digits

//    std::cout << "Checking channel number " << seedDigitIter->ChannelNumber << std::endl;

    Channel seedDigit = *seedDigitIter; // the seed of the cluster

    if(seedDigit.AdcValue >= seedThreshold){
      // ADC above seed : start clustering

    //  debug() << " ---> START NEW CLUSTER WITH SEED @ " << seedDigit->ChannelNumber << endmsg;

      std::vector<Channel>::const_iterator startDigitIter = seedDigitIter; // begin channel of cluster
      std::vector<Channel>::const_iterator stopDigitIter  = seedDigitIter; // end   channel of cluster

      // vector of MCHits that contributed to the cluster
    //  std::vector<const LHCb::MCHit*> clusterHitDistribution;

      // total energy in the cluster from MC
//      double totalEnergyFromMC = 0;

      // map of contributing MCParticles (MCHits) with their relative energy deposit
    //  std::map< const LHCb::MCParticle*, double> mcContributionMap;
    //  std::map< const LHCb::MCHit*, double> mcHitContributionMap;


      //
      // test neighbours to define starting and ending channels of Cluster
      bool ContinueStartLoop = true;
      bool ContinueStopLoop  = true;

      while(((stopDigitIter-startDigitIter)< maxClusterSize )
            && (ContinueStartLoop || ContinueStopLoop) ) {
        // IF cluster size =< maxClusterSize SiPM Channels


        // EXTEND TO THE RIGHT
        if(ContinueStopLoop && ((stopDigitIter+1) != event.end()) ) {
          // IF the next digit exists: try to extend cluster to the 'right'

          const Channel stopDigit = *(stopDigitIter+1); // the next digit

          // next digit should be in the same SiPM, the channel next door, and above adcThreshold
          if((stopDigit.Uplink == seedDigit.Uplink)
             &&(stopDigit.ChannelNumber==((*stopDigitIter).ChannelNumber+1))
             && (stopDigit.AdcValue >=neighbourThreshold)) {

            // if ADC of next digit > the current seed digit, redefine the seed
            if(stopDigit.AdcValue > seedDigit.AdcValue) {

              seedDigitIter = stopDigitIter+1; // increment loop iterator
              seedDigit = *seedDigitIter;      // increment seed channel
    //          debug() << "  ==> Redefined the seed @ " << seedDigit.ChannelNumber << endmsg;

              // set min and max channel of cluster to the new seed (i.e. reset the cluster finding)
              startDigitIter = seedDigitIter;
              stopDigitIter  = seedDigitIter;
              ContinueStartLoop = true;

            } else {
              // IF next digit ADC < current seed, but passes clustering requirements
              stopDigitIter++; // extend cluster to the 'right'
            }


          } else {
            // IF next digit does not satisfy clustering requirements
            ContinueStopLoop = false;
          }

        } else {
          // IF the next digit does not exist in the container (i.e. done with all ClusterCreator)
          ContinueStopLoop = false;
        }


        // So far, we have extended the cluster to the 'right' as far as we could,
        //  redefining the seed as we go (pulling our left-side 'tail' with us).
        // We now need to extend to the 'left' side.


        // EXTEND TO THE LEFT
        if(ContinueStartLoop && ((startDigitIter) != event.begin())) {
          // IF the previous digit exists: try to extend cluster to the 'left'

          const Channel startDigit = *(startDigitIter-1); // the 'previous' digit

          // previous digit should be in the same SiPM, the channel next door, above adcThreshold,
          //  and also after the ending channel of the previous cluster
          if((startDigit.Uplink == seedDigit.Uplink)
             &&(startDigit.ChannelNumber==((*startDigitIter).ChannelNumber-1))
             && (startDigit.AdcValue >=neighbourThreshold)
             &&((startDigitIter-1) > lastStopDigitIter)) {

            // extend cluster to the 'left'
            startDigitIter--;

          } else {
            // previous channel does not satisfy clustering requirements
            ContinueStartLoop = false;
          }

        } else {
          // there is no previous digit in the container
          ContinueStartLoop = false;
        }

      }
      // MaxClusterSize has been reached, or iterator stop was set due to criteria
      // Cluster spans from startDigitIter to stopDigitIter


      // Check if cluster size < maxWidth. If not: choose highest ADC sum for cluster and shrink
      // This is possible because we check the size, then extend right, then extend left, and repeat.
      while((stopDigitIter - startDigitIter) >= maxClusterSize) {
//        debug() << "Cluster size exceeded threshold, Shrinking cluster." << endmsg;
        if((*stopDigitIter).AdcValue > (*startDigitIter).AdcValue)   startDigitIter++;
        else stopDigitIter--;
      }



      lastStopDigitIter = stopDigitIter; // update the 'previous cluster last channel'


        Cluster* currentCluster = new Cluster();
        for(std::vector<Channel>::const_iterator iterDigit = startDigitIter; iterDigit <(stopDigitIter+1); ++iterDigit) {
          currentCluster->AddChannel(iterDigit->Uplink, iterDigit->ChannelNumber, iterDigit->AdcValue);
        }

        // calculate the cluster charge / mean position
  //      debug() << " ---> Done with cluster finding, now calculating charge / frac.pos" << endmsg;

        //
        // Before Cluster is recorded, check that total ADC > threshold and size > minSize

        if( ( currentCluster->GetSumOfAdcValues() >= sumThreshold) &&
          ((stopDigitIter-startDigitIter+1) >= 1) ){


          // Define new cluster
          //  FTCluster::FTCluster( LHCb::FTChannelID &id, double fraction, int size, int charge )



          clusters.push_back(currentCluster);
          foundClusters.push_back(currentCluster);


        } // end of Cluster satisfies charge / size requirements
        else delete currentCluster;
        // Set the loop iterator to the end digit, to continue looking for next cluster without overlap.
        // Will get +1 at end of loop.
        seedDigitIter = stopDigitIter;

      } // END of clustering ( if(seedDigit->adcCount() > m_adcThreshold) )

      // Prepare for next cluster
      ++seedDigitIter;

    } // END of loop over Digits

    return foundClusters;
}
