//
// Name:  SignalShapingDUNE35tTest.h
//
// Purpose: SignalShapingDUNE35tTest module.  Test convolution/deconvolution.
//
// Created:  28-Nov-2011  H. Greenlee

#include <vector>
#include <iostream>
#include <cmath>

#include "art/Framework/Core/ModuleMacros.h" 
#include "art/Framework/Core/EDAnalyzer.h"
#include "art_root_io/TFileService.h"

#include "dune/Utilities/SignalShapingServiceDUNE35t.h"
#include "lardata/Utilities/LArFFT.h"

#include "TComplex.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"

// Local functions.

namespace {

  // Fill histogram from vector (set underflow/overflow bins to zero).

  void vector_to_hist(const std::vector<double>& v, TH1D* h)
  {
    assert(h != 0);
    int nvec = v.size();
    int nbins = h->GetNbinsX();
    int nfill = std::min(nvec, nbins);
    h->SetBinContent(0, 0.);
    int zerobin = h->GetXaxis()->FindBin(0.);
    for(int i=1; i<=nfill; ++i) {
      if(i >= zerobin)
	h->SetBinContent(i, v[i - zerobin]);
      else
	h->SetBinContent(i, v[i - zerobin + nvec]);
    }
    for(int i=nfill+1; i<=nbins+1; ++i)
      h->SetBinContent(i, 0.);
  }

  // Fill vector with initial delta-function at bin d.

  void fill_delta(std::vector<double>& v, int d)
  {
    int n = v.size();
    assert(d >= 0 && d < n);
    for(int i=0; i<n; ++i)
      v[i] = 0.;
    v[d] = 1.;
  }
}

namespace util
{
  class SignalShapingDUNE35tTest : public art::EDAnalyzer
  {
  public:

    // Constructor, destructor.

    explicit SignalShapingDUNE35tTest(fhicl::ParameterSet const& pset);
    virtual ~SignalShapingDUNE35tTest();

    // Overrides.

    void beginJob();
    void analyze(const art::Event& evt);

  };

  DEFINE_ART_MODULE(SignalShapingDUNE35tTest)

  SignalShapingDUNE35tTest::SignalShapingDUNE35tTest(const fhicl::ParameterSet& pset)
  : EDAnalyzer(pset)
  {}

  void SignalShapingDUNE35tTest::beginJob()
  {
    // Get services.

    art::ServiceHandle<art::TFileService> tfs;
    art::ServiceHandle<util::LArFFT> fft;
    art::ServiceHandle<util::SignalShapingServiceDUNE35t> sss;

    int nticks = fft->FFTSize();
    std::cout << "Number of ticks = " << nticks << std::endl;

    // Make collection plane histograms.

    art::TFileDirectory dirc = tfs->mkdir("Collection", "Collection");
    int nhist = std::min(200, nticks);
    int nfilt = std::min(5000, nticks/2);
    int nkern = nticks/2;

    // Make input pulse.

    std::vector<double> tinc(nticks, 0.);
    fill_delta(tinc, nhist/2);
    TH1D* hinc = dirc.make<TH1D>("input", "Collection Input", nhist+1, -0.5, nhist+0.5);
    vector_to_hist(tinc, hinc);

    // Convoluted pulse.

    std::vector<double> tconvc(tinc);
    sss->Convolute(6000, tconvc);
    TH1D* hconvc = dirc.make<TH1D>("conv", "Collection Convoluted", nhist+1, -0.5, nhist+0.5);
    vector_to_hist(tconvc, hconvc);

    // Deconvoluted pulse.

    std::vector<double> tdeconvc(tconvc);
    sss->Deconvolute(6000, tdeconvc);
    TH1D* hdeconvc = dirc.make<TH1D>("deconv", "Collection Deconvoluted", nhist+1, -0.5, nhist+0.5);
    vector_to_hist(tdeconvc, hdeconvc);

    // Get collection response function and fill histogram.

    const std::vector<double>&  respc = sss->SignalShaping(6000).Response();
    TH1D* hrfc = dirc.make<TH1D>("resp", "Collection Response", nhist+1, -nhist/2-0.5, nhist/2+0.5);
    vector_to_hist(respc, hrfc);

    // Get collection convolution kernel and fill histogram.

    const std::vector<TComplex>&  kernc = sss->SignalShaping(6000).ConvKernel();
    std::vector<double> kernrc(kernc.size());
    for(unsigned int i=0; i<kernrc.size(); ++i)
      kernrc[i] = kernc[i].Rho();
    TH1D* hkernc = dirc.make<TH1D>("kern", "Collection Convolution Kernel", nkern+1, -0.5, nkern+0.5);
    hkernc->SetMinimum(0.);
    vector_to_hist(kernrc, hkernc);

    // Get collection filter function and fill histogram.

    const std::vector<TComplex>&  filtc = sss->SignalShaping(6000).Filter();
    std::vector<double> filtrc(filtc.size());
    for(unsigned int i=0; i<filtrc.size(); ++i)
      filtrc[i] = filtc[i].Re();
    TH1D* hffc = dirc.make<TH1D>("filt", "Collection Filter", nfilt+1, -0.5, nfilt+0.5);
    vector_to_hist(filtrc, hffc);

    // Make induction plane histograms.

    art::TFileDirectory diri = tfs->mkdir("Induction", "Induction");

    // Make input pulse.

    std::vector<double> tini(nticks, 0.);
    fill_delta(tini, nhist/2);
    TH1D* hini = diri.make<TH1D>("input", "Induction Input", nhist+1, -0.5, nhist+0.5);
    vector_to_hist(tini, hini);

    // Convoluted pulse.

    std::vector<double> tconvi(tini);
    sss->Convolute(0, tconvi);
    TH1D* hconvi = diri.make<TH1D>("conv", "Induction Convoluted", nhist+1, -0.5, nhist+0.5);
    vector_to_hist(tconvi, hconvi);

    // Deconvoluted pulse.

    std::vector<double> tdeconvi(tconvi);
    sss->Deconvolute(0, tdeconvi);
    TH1D* hdeconvi = diri.make<TH1D>("deconv", "Induction Deconvoluted", nhist+1, -0.5, nhist+0.5);
    vector_to_hist(tdeconvi, hdeconvi);

    // Get induction response function and fill histogram.

    const std::vector<double>&  respi = sss->SignalShaping(0).Response();
    TH1D* hrfi = diri.make<TH1D>("resp", "Induction Response", nhist+1, -nhist/2-0.5, nhist/2+0.5);
    vector_to_hist(respi, hrfi);

    // Get induction convolution kernel and fill histogram.

    const std::vector<TComplex>&  kerni = sss->SignalShaping(0).ConvKernel();
    std::vector<double> kernri(kerni.size());
    for(unsigned int i=0; i<kernri.size(); ++i)
      kernri[i] = kerni[i].Rho();
    TH1D* hkerni = diri.make<TH1D>("kern", "Induction Convolution Kernel", nkern+1, -0.5, nkern+0.5);
    hkerni->SetMinimum(0.);
    vector_to_hist(kernri, hkerni);

    // Get induction filter function and fill histogram.

    const std::vector<TComplex>&  filti = sss->SignalShaping(0).Filter();
    std::vector<double> filtri(filti.size());
    for(unsigned int i=0; i<filtri.size(); ++i)
      filtri[i] = filti[i].Re();
    TH1D* hffi = diri.make<TH1D>("filt", "Induction Filter", nfilt+1, -0.5, nfilt+0.5);
    vector_to_hist(filtri, hffi);
  }

  SignalShapingDUNE35tTest::~SignalShapingDUNE35tTest()
  {}

  void SignalShapingDUNE35tTest::analyze(const art::Event& evt)
 {}
}
