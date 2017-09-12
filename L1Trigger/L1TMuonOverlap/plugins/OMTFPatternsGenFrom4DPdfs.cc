#include <iostream>
#include <strstream>
#include <iomanip>

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "CondFormats/DataRecord/interface/L1TMuonOverlapParamsRcd.h"
#include "CondFormats/L1TObjects/interface/L1TMuonOverlapParams.h"

#include "L1Trigger/L1TMuonOverlap/plugins/OMTFPatternsGenFrom4DPdfs.h"
#include "L1Trigger/L1TMuonOverlap/interface/OMTFinputMaker.h"
#include "L1Trigger/L1TMuonOverlap/interface/OMTFinput.h"
#include "L1Trigger/L1TMuonOverlap/interface/OMTFConfiguration.h"
#include "L1Trigger/L1TMuonOverlap/interface/OMTFConfigMaker.h"
#include "L1Trigger/L1TMuonOverlap/interface/XMLConfigWriter.h"
#include "L1Trigger/L1TMuonOverlap/interface/GoldenPatternPdf4D.h"

#include "L1Trigger/RPCTrigger/interface/RPCConst.h"

#include "SimDataFormats/Track/interface/SimTrack.h"

#include "Math/VectorUtil.h"

#include "L1Trigger/RPCTrigger/interface/RPCConst.h"

#include <TH2F.h>
#include "TFile.h"


OMTFPatternsGenFrom4DPdfs::OMTFPatternsGenFrom4DPdfs(const edm::ParameterSet& cfg):
theConfig(cfg),
g4SimTrackSrc(cfg.getParameter<edm::InputTag>("g4SimTrackSrc")) {
  inputTokenDTPh = consumes<L1MuDTChambPhContainer>(theConfig.getParameter<edm::InputTag>("srcDTPh"));
  inputTokenDTTh = consumes<L1MuDTChambThContainer>(theConfig.getParameter<edm::InputTag>("srcDTTh"));
  inputTokenCSC = consumes<CSCCorrelatedLCTDigiCollection>(theConfig.getParameter<edm::InputTag>("srcCSC"));
  inputTokenRPC = consumes<RPCDigiCollection>(theConfig.getParameter<edm::InputTag>("srcRPC"));
  inputTokenSimHit = consumes<edm::SimTrackContainer>(theConfig.getParameter<edm::InputTag>("g4SimTrackSrc"));
  myInputMaker = new OMTFinputMaker();

  makeGoldenPatterns = theConfig.getParameter<bool>("makeGoldenPatterns");
  makeConnectionsMaps = theConfig.getParameter<bool>("makeConnectionsMaps");
  mergeXMLFiles = theConfig.getParameter<bool>("mergeXMLFiles");
  myOMTFConfig = 0;
}
/////////////////////////////////////////////////////
/////////////////////////////////////////////////////
OMTFPatternsGenFrom4DPdfs::~OMTFPatternsGenFrom4DPdfs(){

  delete myOMTFConfig;
  //delete myOMTFConfigMaker;
  delete processor;

}

/////////////////////////////////////////////////////
/////////////////////////////////////////////////////
void OMTFPatternsGenFrom4DPdfs::beginRun(edm::Run const& run, edm::EventSetup const& iSetup) {

  const L1TMuonOverlapParamsRcd& omtfParamsRcd = iSetup.get<L1TMuonOverlapParamsRcd>();

  edm::ESHandle<L1TMuonOverlapParams> omtfParamsHandle;
  omtfParamsRcd.get(omtfParamsHandle);

  const L1TMuonOverlapParams* omtfParams = omtfParamsHandle.product();
  if (!omtfParams) {
    edm::LogError("L1TMuonOverlapTrackProducer") << "Could not retrieve parameters from Event Setup" << std::endl;
  }
  std::cout<<__FUNCTION__<<"\t"<<__LINE__<<std::endl;
  ///Initialise XML writer with default pdf.
  myWriter = new XMLConfigWriter(myOMTFConfig);

/*  ///For making the patterns use extended pdf width in phi, as pdf are later shifted by the mean value
  ///For low pt muons non shifted pdfs would go out of the default pdf range.
  L1TMuonOverlapParams omtfParamsMutable = *omtfParams;
  std::vector<int> generalParams = *omtfParamsMutable.generalParams();
  nPdfAddrBits = omtfParams->nPdfAddrBits();

  //if(!mergeXMLFiles) //the LUTs are then too big in the GoldenPatternPdf4D
  generalParams[L1TMuonOverlapParams::GENERAL_ADDRBITS] = nPdfAddrBits + 2; //2*nPdfAddrBits;
  omtfParamsMutable.setGeneralParams(generalParams);
  */
  myOMTFConfig->configure(omtfParams);
  processor->configure(myOMTFConfig, omtfParams);

/*  int ptCode = theConfig.getParameter<int>("ptCode"); //assuming that here the legay PAC pt code is given
  int charge = theConfig.getParameter<int>("charge");

  //converting to the uGMT ptCode
  double pt = RPCConst::ptFromIpt(ptCode);
  unsigned int patNum = myOMTFConfig->getPatternNum(pt, charge);
  ptCode = omtfParams->ptLUT()->data(patNum );
  configureProcesor(myOMTFConfig, omtfParams, ptCode, charge, patNum);*/
  //myOMTFConfigMaker = new OMTFConfigMaker(myOMTFConfig);

  std::cout<<"OMTFPatternsGenFrom4DPdfs::beginRun: myOMTFConfig "<<*myOMTFConfig;
}
/////////////////////////////////////////////////////
/////////////////////////////////////////////////////
void OMTFPatternsGenFrom4DPdfs::beginJob(){

  myOMTFConfig = new OMTFConfiguration();
  processor = new PatternsGeneratorProcessor();
}
/////////////////////////////////////////////////////
/////////////////////////////////////////////////////  
void OMTFPatternsGenFrom4DPdfs::endJob(){
  processor->generatePatterns();
  writeGPs();
}
/////////////////////////////////////////////////////
/////////////////////////////////////////////////////
void OMTFPatternsGenFrom4DPdfs::writeGPs() {
  myWriter->initialiseXMLDocument("OMTF");
  const std::vector<GoldenPattern*> & myGPs = processor->getPatterns();

  GoldenPattern *dummy = new GoldenPattern(Key(0,0,0), myOMTFConfig);
  dummy->reset();

  OMTFConfiguration::vector2D mergedPartters = myOMTFConfig->getMergedPartters();
  for(unsigned int iGroup = 0; iGroup < mergedPartters.size(); iGroup++) {
    std::vector<GoldenPattern*> gps(4, dummy);
    for(unsigned int i = 0; i < mergedPartters[iGroup].size(); i++) {
      GoldenPattern* gp = dynamic_cast<GoldenPattern*>(myGPs.at(mergedPartters[iGroup][i]));
      gps[i] =  gp;
    }
    myWriter->writeGPData(*gps[0],*gps[1], *gps[2], *gps[3]);
  }
  std::string fName = "GPs.xml";
  myWriter->finaliseXMLDocument(fName);
}
/////////////////////////////////////////////////////
/////////////////////////////////////////////////////
void OMTFPatternsGenFrom4DPdfs::analyze(const edm::Event& iEvent, const edm::EventSetup& evSetup){
  return;
}
/////////////////////////////////////////////////////
/////////////////////////////////////////////////////  
/*const SimTrack * OMTFPatternsGenFrom4DPdfs::findSimMuon(const edm::Event &ev, const edm::EventSetup &es, const SimTrack * previous){

  const SimTrack * result = 0;
  edm::Handle<edm::SimTrackContainer> simTks;
  ev.getByToken(inputTokenSimHit,simTks);

  for (std::vector<SimTrack>::const_iterator it=simTks->begin(); it< simTks->end(); it++) {
    const SimTrack & aTrack = *it;
    if ( !(aTrack.type() == 13 || aTrack.type() == -13) )continue;
    if(previous && ROOT::Math::VectorUtil::DeltaR(aTrack.momentum(),previous->momentum())<0.07) continue;
    if ( !result || aTrack.momentum().pt() > result->momentum().pt()) result = &aTrack;
  }
  return result;
}*/
/////////////////////////////////////////////////////
/////////////////////////////////////////////////////  
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(OMTFPatternsGenFrom4DPdfs);
