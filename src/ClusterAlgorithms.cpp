//from here
#include "Cluster.h"

//from std
#include <iostream>
#include <vector>




void FindClustersInEvent(std::vector<Cluster*>& clusterVector,
                         const std::vector<Channel>& event,
                         const double neighbourThreshold,
                         const double seedThreshold,
                         const double sumThreshold,
                         const unsigned int maxClusterSize,
                         bool debug){

  if(debug) std::cout << "event has " << event.size() << " channels" << std::endl;
  bool possibleMultipleClusters = false;
  Cluster* currentCluster = nullptr;
  for(std::vector<Channel>::const_iterator chan = event.begin(); chan != event.end(); ++chan){
    if(debug) std::cout << "search at channel " << chan->ChannelNumber << " with adc value " << chan->AdcValue << std::endl;
    std::vector<Channel>::const_iterator storeCurrentChannel;
    if (chan->AdcValue > seedThreshold && currentCluster == nullptr){ //create a new cluster if there is none atm  (&& currentCluster == nullptr should have no effect here)
      currentCluster = new Cluster();
      possibleMultipleClusters = false;
      currentCluster->AddChannel(chan->Uplink, chan->ChannelNumber, chan->AdcValue);
      storeCurrentChannel = chan;
      if(chan != event.begin()){ //add all neighbours in backward direction
        do{
          --chan;
          if(debug) std::cout << "going down to channel " << chan->ChannelNumber << " with adc value " << chan->AdcValue << std::endl;
          if(chan->AdcValue > neighbourThreshold) currentCluster->AddChannel(chan->Uplink, chan->ChannelNumber, chan->AdcValue);
        }while(chan->AdcValue > neighbourThreshold && chan != event.begin());
      }
      chan = storeCurrentChannel;
      if(chan != event.end()-1){ //add all neighbours in backward direction
        do{
          ++chan;
          if(debug) std::cout << "going up to channel " << chan->ChannelNumber << " with adc value " << chan->AdcValue << std::endl;
          if(chan->AdcValue > neighbourThreshold) currentCluster->AddChannel(chan->Uplink, chan->ChannelNumber, chan->AdcValue);
          if(chan->AdcValue > seedThreshold && (chan-1)->AdcValue < seedThreshold) possibleMultipleClusters = true;
        }while(chan->AdcValue > neighbourThreshold && chan != event.end()-1);
      }
    }
    if(currentCluster == nullptr) continue;
    if(currentCluster->GetSumOfAdcValues() > sumThreshold){
      if(debug) std::cout << "storing cluster" << std::endl;
      //if cluster size is greater than max cluster size
      if(currentCluster->GetClusterSize() >= 2*maxClusterSize && possibleMultipleClusters){// try to split the clusters
        std::cout << "cluster splitting needed!" << std::endl;
        }
      if(currentCluster->GetClusterSize() > maxClusterSize){//just reduce the cluster
        currentCluster->Resize(maxClusterSize);
       }
      clusterVector.push_back(currentCluster); //store the cluster
    }
    currentCluster = nullptr; //reset the pointer
  }
}




void FindClustersInEventFTStyle(std::vector<Cluster*>& clusterVector,
                                const std::vector<Channel>& event,
                                const double neighbourThreshold,
                                const double seedThreshold,
                                const double sumThreshold,
                                const unsigned int maxClusterSize,
                                bool debug)
  {
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
          // IF the next digit does not exist in the container (i.e. done with all clusterisation)
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



          clusterVector.push_back(currentCluster);


        } // end of Cluster satisfies charge / size requirements
        else delete currentCluster;
        // Set the loop iterator to the end digit, to continue looking for next cluster without overlap.
        // Will get +1 at end of loop.
        seedDigitIter = stopDigitIter;

      } // END of clustering ( if(seedDigit->adcCount() > m_adcThreshold) )

      // Prepare for next cluster
      ++seedDigitIter;

    } // END of loop over Digits
  }




