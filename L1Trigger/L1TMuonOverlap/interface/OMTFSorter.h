#ifndef OMTF_OMTFSorter_H
#define OMTF_OMTFSorter_H

#include <tuple>
#include <vector>

#include "DataFormats/L1TMuon/interface/RegionalMuonCand.h"
#include "DataFormats/L1TMuon/interface/RegionalMuonCandFwd.h"

#include <L1Trigger/L1TMuonOverlap/interface/GoldenPatternResult.h>
#include "L1Trigger/L1TMuonOverlap/interface/OMTFProcessor.h"
#include "L1Trigger/L1TMuonOverlap/interface/AlgoMuon.h"
#include "L1Trigger/L1TMuonOverlap/interface/GhostBuster.h"

class OMTFSorter{

 public:

  void setNphiBins(unsigned int phiBins) { nPhiBins = phiBins;}

  void sortRefHitResults(const std::vector<IGoldenPattern*>& gPatterns,
        std::vector<AlgoMuon> & refHitCleanCands,
        int charge=0);


  //convert algo muon to outgoing Candidates
  //FIXME (MK): nothing to do with sorter, move it from here!
   std::vector<l1t::RegionalMuonCand> candidates(
                 unsigned int iProcessor, l1t::tftype mtfType,
                 const std::vector<AlgoMuon> & algoCands);

  ///Sort results from a single reference hit.
  ///Select candidate with highest number of hit layers
  ///Then select a candidate with largest likelihood value and given charge
  ///as we allow two candidates with opposite charge from single 10deg region
  AlgoMuon sortRefHitResults(unsigned int iRefHit, const std::vector<IGoldenPattern*>& gPatterns,
				int charge=0);

 private:

  ///Check if the hit pattern of given OMTF candite is not on the list
  ///of invalid hit patterns. Invalid hit patterns provode very little
  ///to efficiency, but gives high contribution to rate.
  ///Candidate with invalid hit patterns is assigned quality=0.
  ///Currently the list of invalid patterns is hardcoded.
  ///This has to be read from configuration.
  bool checkHitPatternValidity(unsigned int hits);

  ///Find a candidate with best parameters for given GoldenPattern
  ///Sorting is made amongs candidates with different reference layers
  ///The output tuple contains (nHitsMax, pdfValMax, refPhi, refLayer, hitsWord, refEta)
  ///hitsWord codes number of layers hit: hitsWord= sum 2**iLogicLayer,
  ///where sum runs over layers which were hit
  //AlgoMuon findBestSingleGpResult(const GoldenPatternResult & aResult); //TODO should be in GoldenPatters


  unsigned int nPhiBins;
};

#endif
