/*
 * PatternsGeneratorProcessor.h
 *
 *  Created on: Aug 3, 2017
 *      Author: kbunkow
 */

#ifndef OMTF_PATTERNSGENERATORPROCESSOR_H_
#define OMTF_PATTERNSGENERATORPROCESSOR_H_

#include "L1Trigger/L1TMuonOverlap/interface/ProcessorBase.h"
#include "L1Trigger/L1TMuonOverlap/interface/GoldenPattern.h"
#include "TROOT.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TString.h"

class PatternsGeneratorProcessor: public ProcessorBase<GoldenPattern> {
public:
  PatternsGeneratorProcessor();
  virtual ~PatternsGeneratorProcessor();

  void generatePatterns();

private:
  void generateDTPattern(int iGroup, int iRefLayer, int iLayer);
  void copyDTPattern(int iGroup, int iRefLayer, int iLayer);
  void generateNonDTPattern(int iGroup, int iRefLayer, int iLayer);
  void generateDTSelfPattern(int iGroup, int iRefLayer, int iLayer);
  void getDTLayers();
  double getMeanSigma(TH2F* h);
  int getMenDistPhiLoadHisto(int iGroup, int iRefLayer, int iLayer, std::vector<TH2F*>* h);
  int getMenDistPhiLoadHisto(int iGroup, int iRefLayer, int iLayer, std::vector<TH1F*>* h);
  TH1F* getHistogram1D(int ipt, int charge, int iLayer, int iReflayer);
  TH2F* getHistogram2D(int ipt, int charge, int iLayer, int iReflayer);

  std::vector<int> layersDT; //logic number of DT layers
  std::vector<int> refLayersDT; //refLayer number of DT layers
  std::string path;
  int emptyHistogramCount = 1000;
};

#endif /* OMTF_PATTERNSGENERATORPROCESSOR_H_ */
