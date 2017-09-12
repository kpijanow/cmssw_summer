#include <iostream>
#include <algorithm>
#include <strstream>

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/FileInPath.h"
#include "FWCore/Utilities/interface/Exception.h"

#include "CondFormats/L1TObjects/interface/L1TMuonOverlapParams.h"

#include "L1Trigger/L1TMuonOverlap/interface/OMTFProcessor.h"
#include "L1Trigger/L1TMuonOverlap/interface/GoldenPattern.h"
#include "L1Trigger/L1TMuonOverlap/interface/OMTFinput.h"
#include "L1Trigger/L1TMuonOverlap/interface/GoldenPatternResult.h"
#include "L1Trigger/L1TMuonOverlap/interface/OMTFSorter.h"
#include "L1Trigger/L1TMuonOverlap/interface/GhostBuster.h"

#include "SimDataFormats/Track/interface/SimTrack.h"
///////////////////////////////////////////////
///////////////////////////////////////////////
OMTFProcessor::OMTFProcessor(): ProcessorBase<GoldenPattern>()  {
  setSorter(new OMTFSorter()); //initialize with the default sorter
  setGhostBuster(new GhostBuster()); //initialize with the default sorter
};

OMTFProcessor::~OMTFProcessor(){

}

///////////////////////////////////////////////
///////////////////////////////////////////////
//const std::vector<OMTFProcessor::resultsMap> &
const void OMTFProcessor::processInput(unsigned int iProcessor,
    const OMTFinput & aInput){
/*  for(auto & itRegion: myResults)
  	for(auto & itKey: itRegion) 
  		itKey.second.reset();*/

  for(auto& itGP: theGPs) {
    for(auto& result : itGP->getResults()) {
      result.reset();
    }
  }

  //////////////////////////////////////
  //////////////////////////////////////  
  std::bitset<128> refHitsBits = aInput.getRefHits(iProcessor);
  if(refHitsBits.none())
    return; // myResults;

  for(unsigned int iLayer=0;iLayer<myOmtfConfig->nLayers();++iLayer) {
    const OMTFinput::vector1D & layerHits = aInput.getLayerData(iLayer);

    if(!layerHits.size()) continue;
    ///Number of reference hits to be checked. 
    unsigned int nTestedRefHits = myOmtfConfig->nTestRefHits();
    for(unsigned int iRefHit = 0; iRefHit < myOmtfConfig->nRefHits(); ++iRefHit) { //loop over all possible refHits, i.e. 128
      if(!refHitsBits[iRefHit]) continue;
      if(nTestedRefHits-- == 0) break;

      const RefHitDef& aRefHitDef = myOmtfConfig->getRefHitsDefs()[iProcessor][iRefHit];

      int phiRef = aInput.getLayerData(myOmtfConfig->getRefToLogicNumber()[aRefHitDef.iRefLayer])[aRefHitDef.iInput];
      int etaRef = aInput.getLayerData(myOmtfConfig->getRefToLogicNumber()[aRefHitDef.iRefLayer],true)[aRefHitDef.iInput];
      unsigned int iRegion = aRefHitDef.iRegion;

      if(myOmtfConfig->getBendingLayers().count(iLayer)) //this iLayer is a banding layer
        phiRef = 0;  //then in the delta_phi in process1Layer1RefLayer one obtains simply the iLayer phi

      const OMTFinput::vector1D restrictedLayerHits = restrictInput(iProcessor, iRegion, iLayer,layerHits);
      //std::cout<<"iLayer "<<iLayer<<" refHitNum "<<myOmtfConfig->nTestRefHits()-nTestedRefHits-1<<" iRefHit "<<iRefHit;
      //std::cout<<" nTestedRefHits "<<nTestedRefHits<<" aRefHitDef "<<aRefHitDef<<std::endl;

      int refLayerLogicNumber = myOmtfConfig->getRefToLogicNumber()[aRefHitDef.iRefLayer];
      int refLayerPhiB = 0;
      if(refLayerLogicNumber < 6) {//is DT layer TODO - check
        refLayerPhiB = aInput.getLayerData(refLayerLogicNumber+1)[aRefHitDef.iInput]; //corresponding bending layer has number +1 versus the phi DT layer - this is configured in the XML
        if(myOmtfConfig->getBendingLayers().count(refLayerLogicNumber+1) == 0) {
          throw cms::Exception("not good: the layer is not bending layer");
        }
      }

      for(auto& itGP: theGPs) {
        if(itGP->key().thePt == 0) //empty pattern
          continue;
        GoldenPattern::layerResult aLayerResult = itGP->process1Layer1RefLayer(aRefHitDef.iRefLayer, iLayer,
            phiRef,
            restrictedLayerHits,
            refLayerPhiB);
        int phiRefSt2 = itGP->propagateRefPhi(phiRef, etaRef, aRefHitDef.iRefLayer);
/*        myResults[myOmtfConfig->nTestRefHits()-nTestedRefHits-1][itGP.second->key()].setRefPhiRHits(aRefHitDef.iRefLayer, phiRef);
        myResults[myOmtfConfig->nTestRefHits()-nTestedRefHits-1][itGP.second->key()].addResult(aRefHitDef.iRefLayer, iLayer,
            aLayerResult.first,
            phiRefSt2, etaRef);*/

        //myResults.at(myOmtfConfig->nTestRefHits()-nTestedRefHits-1).at(itGP->key()).set(aRefHitDef.iRefLayer, phiRefSt2, etaRef, phiRef, iLayer, aLayerResult.first);
        itGP->getResults().at(myOmtfConfig->nTestRefHits()-nTestedRefHits-1).set(aRefHitDef.iRefLayer, phiRefSt2, etaRef, phiRef, iLayer, aLayerResult.first);
      }
    }
  }
  //////////////////////////////////////
  ////////////////////////////////////// 
  {
    /*unsigned int iRefHitNum  = 0;
    for(auto & itRefHit: myResults) {
      for(auto & itKey: itRefHit) {
        itKey.second.finalise();
        if(itKey.second.isValid())
          std::cout<<"iRefHitNum "<<iRefHitNum<<" "<<itKey.first<<"\t"<<itKey.second<<std::endl;
      }
      iRefHitNum++;
    }*/

    for(auto& itGP: theGPs) {
      itGP->finalise();
    }
  }

  std::ostringstream myStr;
  myStr<<"iProcessor: "<<iProcessor<<std::endl;
  myStr<<"Input: ------------"<<std::endl;
  myStr<<aInput<<std::endl; 
  edm::LogInfo("OMTF processor")<<myStr.str();


  return; //myResults;
}   

///////////////////////////////////////////////////////
///////////////////////////////////////////////////////
bool OMTFProcessor::checkHitPatternValidity(unsigned int hits){

  ///FIXME: read the list from configuration so this can be controlled at runtime.
  std::vector<unsigned int> badPatterns = {99840, 34304, 3075, 36928, 12300, 98816, 98944, 33408, 66688, 66176, 7171, 20528, 33856, 35840, 4156, 34880};

  for(auto aHitPattern: badPatterns){
    if(hits==aHitPattern) return false;
  }

  return true;
}
///////////////////////////////////////////////////////
///////////////////////////////////////////////////////

std::vector<AlgoMuon> OMTFProcessor::sortResults(int charge) {
  std::vector<AlgoMuon> algoCandidates = sorter->sortResults(getPatterns(), charge);
  return algoCandidates;
}
///////////////////////////////////////////////////////
///////////////////////////////////////////////////////
std::vector<l1t::RegionalMuonCand> OMTFProcessor::getFinalcandidates(unsigned int iProcessor, l1t::tftype mtfType, const std::vector<AlgoMuon> & algoCands)
{

  std::vector<l1t::RegionalMuonCand> result;

  for(auto myCand: algoCands){
    l1t::RegionalMuonCand candidate;
    candidate.setHwPt(myCand.getPt());
    candidate.setHwEta(myCand.getEta());

    int phiValue = myCand.getPhi();
    if(phiValue>= int(myOmtfConfig->nPhiBins()) )
      phiValue -= myOmtfConfig->nPhiBins();
    ///conversion factor from OMTF to uGMT scale: 5400/576
//    phiValue/=9.375;
    phiValue *= (437./pow(2,12));    // ie. use as in hw: 9.3729977
    candidate.setHwPhi(phiValue);

    candidate.setHwSign(myCand.getCharge()<0 ? 1:0  );
    candidate.setHwSignValid(1);

    unsigned int quality = checkHitPatternValidity(myCand.getFiredLayerBits()) ? 0 | (1 << 2) | (1 << 3)
                                                                     : 0 | (1 << 2);
    if (    abs(myCand.getEta()) == 115
        && (    static_cast<unsigned int>(myCand.getFiredLayerBits()) == std::bitset<18>("100000001110000000").to_ulong()
             || static_cast<unsigned int>(myCand.getFiredLayerBits()) == std::bitset<18>("000000001110000000").to_ulong()
             || static_cast<unsigned int>(myCand.getFiredLayerBits()) == std::bitset<18>("100000000110000000").to_ulong()
             || static_cast<unsigned int>(myCand.getFiredLayerBits()) == std::bitset<18>("100000001100000000").to_ulong()
             || static_cast<unsigned int>(myCand.getFiredLayerBits()) == std::bitset<18>("100000001010000000").to_ulong()
           )
       ) quality =4;

//  if (abs(myCand.getEta()) == 121) quality = 4;
    if (abs(myCand.getEta()) == 121) quality = 0; // changed on request from HI

    candidate.setHwQual (quality);

    std::map<int, int> trackAddr;
    trackAddr[0] = myCand.getFiredLayerBits();
    trackAddr[1] = myCand.getRefLayer();
    trackAddr[2] = myCand.getDisc();
    candidate.setTrackAddress(trackAddr);
    candidate.setTFIdentifiers(iProcessor,mtfType);
    if (candidate.hwPt() >= 0)  result.push_back(candidate);
  }
  return result;
}
///////////////////////////////////////////////////////
///////////////////////////////////////////////////////

