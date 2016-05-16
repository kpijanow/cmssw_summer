
#ifndef EMTFTokens_h
#define EMTFTokens_h

#include "DataFormats/L1TMuon/interface/RegionalMuonCand.h"
#include "DataFormats/L1TMuon/interface/EMTFOutput.h"
#include "DataFormats/CSCDigi/interface/CSCCorrelatedLCTDigiCollection.h"

#include "EventFilter/L1TRawToDigi/interface/PackerTokens.h"

namespace l1t {
  namespace stage2 {
    class EMTFTokens : public PackerTokens {
    public:
      EMTFTokens(const edm::ParameterSet&, edm::ConsumesCollector&);
      
      inline const edm::EDGetTokenT<RegionalMuonCandBxCollection>& getRegionalMuonCandToken() const { return regionalMuonCandToken_; };
      inline const edm::EDGetTokenT<EMTFOutputCollection>& getEMTFOutputToken() const { return EMTFOutputToken_; };
      inline const edm::EDGetTokenT<CSCCorrelatedLCTDigiCollection>& getEMTFLCTToken() const { return EMTFLCTToken_; }
      
    private:
      
      edm::EDGetTokenT<RegionalMuonCandBxCollection> regionalMuonCandToken_;
      edm::EDGetTokenT<EMTFOutputCollection> EMTFOutputToken_;
      edm::EDGetTokenT<CSCCorrelatedLCTDigiCollection> EMTFLCTToken_;
      
    };
  }
}

#endif
