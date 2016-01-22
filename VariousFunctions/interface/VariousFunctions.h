#ifndef ggA_GenLevel_Analyzer_VariousFunctions_interface_VariousFunctions_h
#define ggA_GenLevel_Analyzer_VariousFunctions_interface_VariousFunctions_h

#include <vector>
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "TCanvas.h"
#include "TGraphAsymmErrors.h"
#include "TAxis.h"
#include <string>
#include "TH2F.h"

class VariousFunctions { 

  public:
    static reco::GenParticleRef findDaughterInDaughters(const reco::GenParticleRef& , const double, const bool);
    static bool findIfInDaughters(const reco::GenParticleRef& , const double, const bool);
    static int tauDecayMode(const reco::GenParticleRef &);
    static reco::LeafCandidate::LorentzVector sumTauP4(const reco::GenParticleRef&, const int, const bool);
    static double getDiTauDR(const reco::GenParticleRef&, const reco::GenParticleRef&, const bool);
    static double getDiThingDR_1(const reco::GenParticleRef&, const reco::GenParticleRef&);
    static void formatAndDrawCanvasAndHist2D(TCanvas&, TH2F*, 
                                           const Int_t, const Int_t, const Int_t, 
                                           const Color_t, const Size_t, const Style_t,
                                           const char*, const Float_t, const Float_t, const Float_t, 
				           const char*, const Float_t, const Float_t, const Float_t,
  					   const char*, const Float_t, const Float_t, const Float_t);
    static double getHigherPt(const reco::GenParticleRef&, const reco::GenParticleRef&);
    static double getLowerPt(const reco::GenParticleRef&, const reco::GenParticleRef&);
};
#endif
