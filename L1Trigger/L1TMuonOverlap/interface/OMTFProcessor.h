#ifndef OMTF_OMTFProcessor_H
#define OMTF_OMTFProcessor_H

#include <map>

#include "L1Trigger/L1TMuonOverlap/interface/IGoldenPattern.h"
#include "L1Trigger/L1TMuonOverlap/interface/GoldenPattern.h"
#include <L1Trigger/L1TMuonOverlap/interface/GoldenPatternResult.h>
#include "L1Trigger/L1TMuonOverlap/interface/OMTFConfiguration.h"


class L1TMuonOverlapParams;
class OMTFinput;

class SimTrack;

namespace edm{
class ParameterSet;
}

class OMTFProcessor {

 public:

  //typedef std::map<Key,GoldenPatternResult> resultsMap;

  OMTFProcessor() {
    myOmtfConfig = 0;
  };

  ~OMTFProcessor();
  
  ///Just sets the myOmtfConfig
  void configure(const OMTFConfiguration* omtfParams) {
    myOmtfConfig = omtfParams;
  }

  ///Fill GP map with patterns from CondFormats object
  bool configure(const OMTFConfiguration * omtfParams, const L1TMuonOverlapParams * omtfPatterns);

  ///Process input data from a single event
  ///Input data is represented by hits in logic layers expressed in local coordinates
  ///Vector index: number of the ref hit (from 0 to nTestRefHits i.e. 4)
  ///Map key: GoldenPattern key
  //const std::vector<OMTFProcessor::resultsMap> &
  const void processInput(unsigned int iProcessor,
							      const OMTFinput & aInput);

  ///Return vector of GoldenPatterns
  const std::vector<IGoldenPattern*> & getPatterns() const;

  ///Fill counts for a GoldenPattern of this
  ///processor unit. Pattern key is selcted according 
  ///to the SimTrack parameters.
  void fillCounts(unsigned int iProcessor,
		  const OMTFinput & aInput,
		  const SimTrack* aSimMuon);

  ///Average patterns. Use same meanDistPhi for two
  ///patterns neighboring in pt code.
  ///Averaging is made saparately fo each charge
  void averagePatterns(int charge);
  
  ///Add GoldenPattern to pattern map.
  ///If GP key already exists in map, a new entry is ignored
  bool addGP(IGoldenPattern *aGP);

 private:
  ///Reset all configuration parameters
  void resetConfiguration();



  ///Shift pdf indexes by differecne between averaged and
  ///original meanDistPhi
  void shiftGP(GoldenPattern *aGP,
	       const GoldenPattern::vector2D & meanDistPhiNew,
	       const GoldenPattern::vector2D & meanDistPhiOld);

  ///Fill map of used inputs.
  ///FIXME: using hack from OMTFConfiguration
  void fillInputRange(unsigned int iProcessor,
		      unsigned int iCone,
		      const OMTFinput & aInput);

  void fillInputRange(unsigned int iProcessor,
		      unsigned int iCone,
		      unsigned int iRefLayer,
		      unsigned int iHit);
    
  ///Remove hits whis are outside input range
  ///for given processor and cone
  OMTFinput::vector1D restrictInput(unsigned int iProcessor,
				    unsigned int iCone,
				    unsigned int iLayer,
				    const OMTFinput::vector1D & layerHits);

  ///vector holding Golden Patterns
  std::vector<IGoldenPattern*> theGPs;

  ///Map holding results on current event data
  ///for each GP. 
  ///Reference hit number is isued as a vector index.
  //std::vector<OMTFProcessor::resultsMap> myResults;

  ///Configuration of the algorithm. This object
  ///does not contain the patterns data.
  const OMTFConfiguration  * myOmtfConfig;
};


#endif
