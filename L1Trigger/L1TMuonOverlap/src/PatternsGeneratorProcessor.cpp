/*
 * PatternsGeneratorProcessor.cpp
 *
 *  Created on: Aug 3, 2017
 *      Author: kbunkow
 */

#include "L1Trigger/L1TMuonOverlap/interface/PatternsGeneratorProcessor.h"
#include <cmath>
#include <vector>
#include <algorithm>
#include "TROOT.h"
#include "TString.h"
#include "TStyle.h"
#include "TH2D.h"
#include "TH1F.h"
#include "TFile.h"
#include "TF1.h"
#include "TMath.h"
#include "TTree.h"
#include "TSystem.h"
#include <sys/stat.h>
#include <unistd.h>
#include <string>
#include <sstream>
#include "L1Trigger/RPCTrigger/interface/RPCConst.h"

PatternsGeneratorProcessor::PatternsGeneratorProcessor(): ProcessorBase<GoldenPattern>() {
  // TODO Auto-generated constructor stub

}

PatternsGeneratorProcessor::~PatternsGeneratorProcessor() {
  // TODO Auto-generated destructor stub
}

void PatternsGeneratorProcessor::generatePatterns() {
	path = "/afs/cern.ch/work/k/kpijanow/public/CMSSW_9_2_0/src/L1Trigger/L1TMuonOverlap/test/expert/";
	OMTFConfiguration::vector2D mergedPartters = myOmtfConfig->getMergedPartters();
	std::cout<<"#####   gneratePatterns   #####"<<std::endl;

	for(auto& it : theGPs){
		it->reset();
	}

	//TH1::AddDirectory(kFALSE);
	getDTLayers();

	for(unsigned int iGroup = 0; iGroup < mergedPartters.size(); iGroup++)
	{
		///Mean dist phi data
		for(unsigned int iLayer=0;iLayer<myOmtfConfig->nLayers();++iLayer){
		  for(unsigned int iRefLayer=0;iRefLayer<myOmtfConfig->nRefLayers();++iRefLayer)
		  {
          //check if refLayer is a DT
          if(std::find(refLayersDT.begin(), refLayersDT.end(), iRefLayer) != refLayersDT.end()){
            //histograms of DT layers vs DT_bending
            if(myOmtfConfig->getRefToLogicNumber().at(iRefLayer) + 1 == (int)iLayer){
              generateDTSelfPattern(iGroup, iRefLayer, iLayer);
            }
            else{
              generateDTPattern(iGroup, iRefLayer, iLayer);
              //copyDTPattern(iGroup, iRefLayer, iLayer);
            }
          }
          else{
            generateNonDTPattern(iGroup, iRefLayer, iLayer);
          }
			  //}
		  }
		}
		std::cout<<"######\t\tGROUP"<<iGroup<<"\t\t######"<<std::endl;
	}
	std::cout<<"######\t\t END \t\t######"<<std::endl;
}

void PatternsGeneratorProcessor::generateDTPattern(int iGroup, int iRefLayer, int iLayer){
	OMTFConfiguration::vector2D mergedPartters = myOmtfConfig->getMergedPartters();
	int iMerge = mergedPartters[iGroup].size();
	std::vector<TH2F*> h;
	bool empty = true;

	for(int k = 0; k < iMerge; k++)
	{
    unsigned int ipt = theGPs.at(mergedPartters[iGroup][k])->key().thePt;
    int charge = theGPs.at(mergedPartters[iGroup][k])->key().theCharge;

		h.push_back(getHistogram2D(ipt, charge,
				iLayer, iRefLayer));
		TH2F* hRef = getHistogram2D(ipt, charge, myOmtfConfig->getRefToLogicNumber().at(iRefLayer), iRefLayer);
    if(h.at(k)->GetSum() < emptyHistogramCount){ //if this happens then histogram is empty
      hRef->Delete();
      continue;
    }

		double sigY = getMeanSigma(h.at(k)) * 0.8; 		//sigma in Y direction, multiplication by 0.8 is empirical

		const float logC = log(h.at(k)->GetSum()/hRef->GetSum() * 1.0 / sqrt(2 * TMath::Pi()) / sigY);
		const float minPlog =  log(myOmtfConfig->minPdfVal());
		const int nPdfValBits = myOmtfConfig->nPdfValBits();

		///Pdf data
		for(unsigned int iPdf=0;iPdf<(unsigned int)(1<<myOmtfConfig->nPdfAddrBits());++iPdf)
		{
			float pVal = logC - (iPdf - (int)(1<<myOmtfConfig->nPdfAddrBits())/2) * (iPdf - (int)(1<<myOmtfConfig->nPdfAddrBits())/2)/2.0/sigY/sigY;
			if(pVal<minPlog || hRef->GetSum() == 0){
				theGPs.at(mergedPartters[iGroup][k])->setPdfValue(0, iLayer, iRefLayer, iPdf);
			  continue;
			}

			int digitisedVal = rint((std::pow(2,nPdfValBits)-1) - (pVal/minPlog)*(std::pow(2,nPdfValBits)-1));
			theGPs.at(mergedPartters[iGroup][k])->setPdfValue(digitisedVal, iLayer, iRefLayer, iPdf);
		}
		hRef->Delete();
		empty = false;
	}

	if(!empty)
	{
    //merging histograms
	  TH2F* hMerge = new TH2F();
	  unsigned int i  = 0;
	  //find non empty histogram
	  while (h.size() > i && !hMerge->GetSum())
	  {
	    if(h.at(i)->GetSum()) {
	      hMerge->Delete();
	      hMerge = (TH2F*)h.at(0)->Clone("hMerge");
	    }
	    i++;
	  }
    hMerge->Reset();
    for(int i = 0; i < iMerge; i++)
    {
      if(h.at(i)->GetSum() && h.at(i)->GetXaxis()->GetXmax() == hMerge->GetXaxis()->GetXmax())
        hMerge->Add(h.at(i));
    }

    //fitting
    TF1 *f1 = new TF1("f1","pol1");
    hMerge->Fit("f1", "QN");

    //function: y = A * x + B
    double B = f1->GetParameter(0);
    double A = f1->GetParameter(1);

    //saving calculated values to theGPs
    for(int k = 0; k < iMerge; k++)
    {
      int value = rint(B) - h.at(k)->GetNbinsX()/2 + h.at(k)->GetNbinsX()/2 * A;
      theGPs.at(mergedPartters[iGroup][k])->setMeanDistPhiValue(value, iLayer, iRefLayer, 0);
      value = rint(A * pow(2, myOmtfConfig->nPdfAddrBits())); //multiply A by power of 2 to save it as integer
      theGPs.at(mergedPartters[iGroup][k])->setMeanDistPhiValue(value, iLayer, iRefLayer, 1);
    }
    hMerge->Delete();
    f1->Delete();
	}
	for(auto &a : h)
	  a->Delete();
	h.clear();
}

void PatternsGeneratorProcessor::generateNonDTPattern(int iGroup, int iRefLayer, int iLayer){
	OMTFConfiguration::vector2D mergedPartters = myOmtfConfig->getMergedPartters();
	int iMerge = mergedPartters[iGroup].size();
	std::vector<TH1F*> h;
	int meanDistPhi = getMenDistPhiLoadHisto(iGroup, iRefLayer, iLayer, &h);

  //save PDF's
	for(int k = 0; k < iMerge; k++)
	{
		unsigned int ipt = theGPs.at(mergedPartters[iGroup][k])->key().thePt;
		int charge = theGPs.at(mergedPartters[iGroup][k])->key().theCharge;

		TH1F* hRef = getHistogram1D(ipt, charge, myOmtfConfig->getRefToLogicNumber().at(iRefLayer), iRefLayer);
		if(h.at(k)->GetSum() < emptyHistogramCount) //if this happens then histogram is empty
    {
		  hRef->Delete();
      continue;
    }

		theGPs.at(mergedPartters[iGroup][k])->setMeanDistPhiValue(meanDistPhi - h.at(k)->GetNbinsX()/2, iLayer, iRefLayer, 0);

		const float minPlog =  log(myOmtfConfig->minPdfVal());
		const int nPdfValBits = myOmtfConfig->nPdfValBits();
		///Pdf data
		for(unsigned int iPdf=0;iPdf<(unsigned int)(1<<myOmtfConfig->nPdfAddrBits());++iPdf)
		{
			// +1 in GetBinContent is due to bin numbers in ROOT (bin = 0; underflow bin)
			float pVal = log(h.at(k)->GetBinContent(meanDistPhi - exp2(myOmtfConfig->nPdfAddrBits()-1) + iPdf + 1)/hRef->GetSum());
			if(pVal < minPlog || hRef->GetSum() == 0){
				theGPs.at(mergedPartters[iGroup][k])->setPdfValue(0, iLayer, iRefLayer, iPdf);
			  continue;
			}
			int digitisedVal = rint((std::pow(2,nPdfValBits)-1) - (pVal/minPlog)*(std::pow(2,nPdfValBits)-1));
			theGPs.at(mergedPartters[iGroup][k])->setPdfValue(digitisedVal, iLayer, iRefLayer, iPdf);
		}
		hRef->Delete();
	}

  for(auto &a : h)
    a->Delete();
  h.clear();
}

void PatternsGeneratorProcessor::getDTLayers(){
	const std::map<int, int> hwToLogicLayer = myOmtfConfig->getHwToLogicLayer();
	layersDT.clear();
	refLayersDT.clear();

	layersDT.push_back(hwToLogicLayer.at(101));
	layersDT.push_back(hwToLogicLayer.at(102));
	layersDT.push_back(hwToLogicLayer.at(103));

	for(auto& L : layersDT){
		for(unsigned int i = 0; i < myOmtfConfig->getRefToLogicNumber().size(); i++){
			if(myOmtfConfig->getRefToLogicNumber().at(i) == L)
				refLayersDT.push_back(i);
		}
	}
}

void PatternsGeneratorProcessor::generateDTSelfPattern(int iGroup, int iRefLayer, int iLayer){
  OMTFConfiguration::vector2D mergedPartters = myOmtfConfig->getMergedPartters();
  int iMerge = mergedPartters[iGroup].size();
  std::vector<TH2F*> h;
  int meanDistPhi = getMenDistPhiLoadHisto(iGroup, iRefLayer, iLayer, &h);

  for(int k = 0; k < iMerge; k++)
  {
    TH1D* hProj = (TH1D*)h.at(k)->ProjectionY("hProj");

    theGPs.at(mergedPartters[iGroup][k])->setMeanDistPhiValue(meanDistPhi - hProj->GetNbinsX()/2, iLayer, iRefLayer, 0);

    const float minPlog =  log(myOmtfConfig->minPdfVal());
    const int nPdfValBits = myOmtfConfig->nPdfValBits();

    ///Pdf data
    for(unsigned int iPdf=0;iPdf<(unsigned int)(1<<myOmtfConfig->nPdfAddrBits());++iPdf)
    {
      // +1 in GetBinContent is due to bin numbers in ROOT (bin = 0; underflow bin)
      float pVal = log(hProj->GetBinContent(meanDistPhi - exp2(myOmtfConfig->nPdfAddrBits()-1) + iPdf + 1)/h.at(k)->GetSum());
      if(pVal < minPlog || hProj->GetSum() == 0){
        theGPs.at(mergedPartters[iGroup][k])->setPdfValue(0, iLayer, iRefLayer, iPdf);
        continue;
      }

      int digitisedVal = rint((std::pow(2,nPdfValBits)-1) - (pVal/minPlog)*(std::pow(2,nPdfValBits)-1));
      theGPs.at(mergedPartters[iGroup][k])->setPdfValue(digitisedVal, iLayer, iRefLayer, iPdf);
    }

    delete hProj;
  }

  for(auto &a : h)
    delete a;
  h.clear();
}

void PatternsGeneratorProcessor::copyDTPattern(int iGroup, int iRefLayer, int iLayer){
  OMTFConfiguration::vector2D mergedPartters = myOmtfConfig->getMergedPartters();
  int iMerge = mergedPartters[iGroup].size();
  std::vector<TH2F*> h;
  int meanDistPhi = getMenDistPhiLoadHisto(iGroup, iRefLayer, iLayer, &h);

  for(int k = 0; k < iMerge; k++)
  {
    unsigned int ipt = theGPs.at(mergedPartters[iGroup][k])->key().thePt;
    int charge = theGPs.at(mergedPartters[iGroup][k])->key().theCharge;

    TH2F* hRef = getHistogram2D(ipt, charge, myOmtfConfig->getRefToLogicNumber().at(iRefLayer), iRefLayer);
    if(h.at(k)->GetSum() < emptyHistogramCount){ //if this happens then histogram is empty
      delete hRef;
      continue;
    }
    TH1D* hProj = (TH1D*)h.at(k)->ProjectionY("hProj");
    theGPs.at(mergedPartters[iGroup][k])->setMeanDistPhiValue(meanDistPhi - h.at(k)->GetNbinsX()/2, iLayer, iRefLayer, 0);
    const float minPlog =  log(myOmtfConfig->minPdfVal());
    const int nPdfValBits = myOmtfConfig->nPdfValBits();

    ///Pdf data
    for(unsigned int iPdf=0;iPdf<(unsigned int)(1<<myOmtfConfig->nPdfAddrBits());++iPdf)
    {
      // +1 in GetBinContent is due to bin numbers in ROOT (bin = 0; underflow bin)
      float pVal = log(hProj->GetBinContent(meanDistPhi - exp2(myOmtfConfig->nPdfAddrBits()-1) + iPdf + 1)/hRef->GetSum());
      if(pVal < minPlog || hRef->GetSum() == 0){
        theGPs.at(mergedPartters[iGroup][k])->setPdfValue(0, iLayer, iRefLayer, iPdf);
        continue;
      }
      int digitisedVal = rint((std::pow(2,nPdfValBits)-1) - (pVal/minPlog)*(std::pow(2,nPdfValBits)-1));
      theGPs.at(mergedPartters[iGroup][k])->setPdfValue(digitisedVal, iLayer, iRefLayer, iPdf);
    }
    delete hRef;
    delete hProj;
  }

  for(auto &a : h)
    delete a;
  h.clear();
}

int PatternsGeneratorProcessor::getMenDistPhiLoadHisto(int iGroup, int iRefLayer, int iLayer, std::vector<TH2F*>* h)
{
  int meanDistPhi = 0;
  int meanDistCount = 0;
  OMTFConfiguration::vector2D mergedPartters = myOmtfConfig->getMergedPartters();
  int iMerge = mergedPartters[iGroup].size();
  //read histograms form a group and calculate meanDistPhi
  for(int k = 0; k < iMerge; k++)
  {
    unsigned int ipt = theGPs.at(mergedPartters[iGroup][k])->key().thePt;
    int charge = theGPs.at(mergedPartters[iGroup][k])->key().theCharge;
    (*h).push_back(getHistogram2D(ipt, charge, iLayer, iRefLayer));
    if((*h).at(k)->GetSum() < emptyHistogramCount){ //if this happens then histogram is considered empty
      continue;
    }
    meanDistPhi += (*h).at(k)->GetMean()* (*h).at(k)->GetSum();
    meanDistCount += (*h).at(k)->GetSum();
  }
  if(meanDistCount)
    meanDistPhi = rint(meanDistPhi/meanDistCount);

  return meanDistPhi;
}

int PatternsGeneratorProcessor::getMenDistPhiLoadHisto(int iGroup, int iRefLayer, int iLayer, std::vector<TH1F*>* h)
{
  int meanDistPhi = 0;
  int meanDistCount = 0;
  OMTFConfiguration::vector2D mergedPartters = myOmtfConfig->getMergedPartters();
  int iMerge = mergedPartters[iGroup].size();
  //read histograms form a group and calculate meanDistPhi
  for(int k = 0; k < iMerge; k++)
  {
    unsigned int ipt = theGPs.at(mergedPartters[iGroup][k])->key().thePt;
    int charge = theGPs.at(mergedPartters[iGroup][k])->key().theCharge;
    (*h).push_back(getHistogram1D(ipt, charge, iLayer, iRefLayer));
    if((*h).at(k)->GetSum() < emptyHistogramCount){ //if this happens then histogram is considered empty
      continue;
    }
    meanDistPhi += (*h).at(k)->GetMean()* (*h).at(k)->GetSum();
    meanDistCount += (*h).at(k)->GetSum();
  }
  if(meanDistCount)
    meanDistPhi = rint(meanDistPhi/meanDistCount);

  return meanDistPhi;
}

TH1F* PatternsGeneratorProcessor::getHistogram1D(int ipt, int charge, int iLayer, int iRefLayer)
{
  std::ostringstream s1;
  s1 << path << "bendingDistr_ptCode_" << ipt << "_ch_" << charge << ".root";
  TH1F* h = new TH1F();

  if(access( (s1.str()).c_str(), F_OK ) == -1){
    std::cout<<"NO ACCESS to "<<s1.str().c_str()<<std::endl;
    return h;
  }

  TFile* file = TFile::Open((TString)s1.str());
  h->Delete();
  h = (TH1F*)file->Get(Form("ipt_%i_ch%i_layer_%i_refLayer_%i", ipt, charge,
      iLayer, iRefLayer));
  h->SetDirectory(0);
  delete file;
  return h;
}

TH2F* PatternsGeneratorProcessor::getHistogram2D(int ipt, int charge, int iLayer, int iRefLayer)
{
  std::ostringstream s1;
  s1 << path << "bendingDistr_ptCode_" << ipt << "_ch_" << charge << ".root";
  TH2F* h = new TH2F();

  if(access( (s1.str()).c_str(), F_OK ) == -1){
    std::cout<<"NO ACCESS to "<<s1.str().c_str()<<std::endl;
    return h;
  }

  TFile* file = TFile::Open((TString)s1.str());
  h->Delete();
  h = (TH2F*)file->Get(Form("ipt_%i_ch%i_layer_%i_refLayer_%i", ipt, charge,
      iLayer, iRefLayer));
  h->SetDirectory(0);
  delete file;
  return h;
}

double PatternsGeneratorProcessor::getMeanSigma(TH2F* h)
{
  //mean sigma from slices
  TObjArray aSlices;
  h->FitSlicesY(0, 0, -1, 0, "QNR", &aSlices);

  double meanSigma = 0;
  double denom = 0;
  TH1D * as = (TH1D*)aSlices[2]->Clone("as");
  for(int d = 0; d < as->GetNbinsX(); d++)
  {
    //error of the value of sigma cannot be grater than 20%
    if(as->GetBinError(d)/as->GetBinContent(d) < 0.2 && as->GetBinError(d) > 0 && as->GetBinContent(d) > 0)
    {
      meanSigma += as->GetBinContent(d) * 1/as->GetBinError(d);
      denom += 1/as->GetBinError(d);
    }
    else
    {
      as->SetBinContent(d, 0);
      as->SetBinError(d, 0);
    }
  }

  if(denom > 0)
    meanSigma = meanSigma / denom;

  //if sigma is lower than 1 bin than all values are in one bin ie. it is a histogram of DT vs itself
  //the value was chosen to match PDF values from Arthur's patterns
  if (meanSigma < 0.25)
    meanSigma = 0.25;

  delete as;
  return meanSigma;
}
