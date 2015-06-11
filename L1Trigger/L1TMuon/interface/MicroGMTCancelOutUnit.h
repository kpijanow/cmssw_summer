#ifndef __l1microgmtcanceloutunit_h
#define __l1microgmtcanceloutunit_h

#include "MicroGMTConfiguration.h"
#include "MicroGMTMatchQualLUT.h"

namespace l1t {
  class MicroGMTCancelOutUnit {
    public: 
      explicit MicroGMTCancelOutUnit (const edm::ParameterSet&);
      virtual ~MicroGMTCancelOutUnit ();
      /// Cancel out between sectors/wedges in one track finder
      void setCancelOutBits(L1TGMTInternalWedges&, tftype trackFinder);
      /// Cancel-out between overlap and barrel track finders
      void setCancelOutBitsOverlapBarrel(L1TGMTInternalWedges&, L1TGMTInternalWedges&);
      /// Cancel-out between overlap and endcap track finders
      void setCancelOutBitsOverlapEndcap(L1TGMTInternalWedges&, L1TGMTInternalWedges&);
    private:
      /// Compares all muons from coll1 with all muons from coll2 and sets the cancel-bits
      void getCancelOutBits(std::vector<L1TGMTInternalMuon*>& coll1, std::vector<L1TGMTInternalMuon*>& coll2);

      MicroGMTMatchQualLUT m_boPosMatchQualLUT;
      MicroGMTMatchQualLUT m_boNegMatchQualLUT;
      MicroGMTMatchQualLUT m_foPosMatchQualLUT;
      MicroGMTMatchQualLUT m_foNegMatchQualLUT;
      MicroGMTMatchQualLUT m_brlSingleMatchQualLUT;
      MicroGMTMatchQualLUT m_ovlPosSingleMatchQualLUT;
      MicroGMTMatchQualLUT m_ovlNegSingleMatchQualLUT;
      MicroGMTMatchQualLUT m_fwdPosSingleMatchQualLUT;
      MicroGMTMatchQualLUT m_fwdNegSingleMatchQualLUT;
      std::map<int, MicroGMTMatchQualLUT*> m_lutDict;
  };
}
#endif /* defined(__l1microgmtcanceloutunit_h) */