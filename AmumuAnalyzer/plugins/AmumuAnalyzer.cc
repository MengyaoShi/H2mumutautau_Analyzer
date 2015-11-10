// -*- C++ -*-
//
// Package:    Amumu/AmumuAnalyzer
// Class:      AmumuAnalyzer
// 
/**\class AmumuAnalyzer AmumuAnalyzer.cc Amumu/AmumuAnalyzer/plugins/AmumuAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Mengyao Shi
//         Created:  Mon, 19 Oct 2015 10:59:36 GMT
//
//


// system include files
#include <memory>
#include <map>
#include <string>
#include <vector>
#include "TH1D.h"
#include "TH2D.h"
#include <math.h>
#include <sstream>
#include <typeinfo>
#include "TCanvas.h"
#include "TFile.h"
// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
 #include "DataFormats/Candidate/interface/Candidate.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "ggA_GenLevel_Analyzer/VariousFunctions/interface/VariousFunctions.h"
//
// class declaration
//
using namespace std;
using namespace edm;
using namespace reco;

class AmumuAnalyzer : public edm::EDAnalyzer {
   public:
      explicit AmumuAnalyzer(const edm::ParameterSet&);
      ~AmumuAnalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
      
   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;
      void reset(const bool);
      double mass;
      TH2F* DeltaR_a_and_tau_VS_pt_of_a_;
      TFile* out_;
      std::string outFileName_;
      edm::InputTag genParticleTag_;
      std::map<std::string, TH1D*> histos1D_;
      std::map<std::string, TH2D*> histos2D_;

      //std::string jetOutputFileName_;
      //ofstream jetOutput_;
      //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
      //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

      // ----------member data ---------------------------
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
//
AmumuAnalyzer::AmumuAnalyzer(const edm::ParameterSet& iConfig):
 outFileName_(iConfig.getParameter<std::string>("outFileName")),
 //jetOutputFileName_(iConfig.getParameter<std::string>("jetOutputFileName")),
 genParticleTag_(iConfig.getParameter<edm::InputTag>("genParticleTag")),
 histos1D_(),
 histos2D_()
{
  reset(false);
}
AmumuAnalyzer::~AmumuAnalyzer()
{
   reset(true);
   // do anything here that needs to be done at desctruction time
   //    // (e.g. close files, deallocate resources etc.)
   //
}

void AmumuAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  //now do what ever initialization is needed
  edm::Handle<reco::GenParticleCollection> pGenParticles;
  iEvent.getByLabel(genParticleTag_, pGenParticles);
  for(reco::GenParticleCollection::const_iterator iGenParticle = pGenParticles->begin(); iGenParticle != pGenParticles->end(); ++iGenParticle){
    if((*iGenParticle).pdgId() == 35 && ((*iGenParticle).numberOfDaughters()==2) )
    {
      reco::GenParticleRef child0 = iGenParticle->daughterRef(0);
      reco::GenParticleRef child1 = iGenParticle->daughterRef(1);
      //std::cout << "\tmom pdgid=" << iGenParticle->pdgId() << "\tchild0 pdgid=" << child0->pdgId() <<"\tchild1 pdgid=" << child1->pdgId() << std::endl;
      if(child0->pdgId() == 36 && child1->pdgId()==36 && (child0->numberOfDaughters()==2)&&(child1->numberOfDaughters()==2) )
      { 
        reco::GenParticleRef childchild00=child0->daughterRef(0);
        reco::GenParticleRef childchild01=child0->daughterRef(1);
        reco::GenParticleRef childchild10=child1->daughterRef(0);
        reco::GenParticleRef childchild11=child1->daughterRef(1);
        if(fabs(childchild00->pdgId())==13 &&( fabs(childchild10->pdgId())==15))
        {
         // std::cout << "\tthen mom pdgid=" << child0->pdgId() << "\tchildchild00" << childchild00->pdgId() << "\tchildchild01" << childchild01->pdgId() <<std::endl;
          //std::cout << "\tAnd mom pdgid=" << child1->pdgId() << "\tchildchild10" << childchild10->pdgId() << "\tchildchild11" << childchild11->pdgId() <<std::endl;
          //std::cout << "\tStatus of Tau" << childchild10->status() << std::endl;
          //std::cout << "\tStatus of Mu" << childchild00->status() << std::endl;
          mass=0.0;
 	  mass=child0->mass();     
          double DR_tau=VariousFunctions::getDiThingDR_1(childchild10, childchild11);
          double DR_mu=VariousFunctions::getDiThingDR_1(childchild00, childchild01);
          double DR_a_and_tau=VariousFunctions::getDiThingDR_1(child1, childchild11);
          double p_of_a=child1->pt();
          if(DR_tau > 2.0)
          {
            std::cout<<"eta" << childchild10->eta() << "phi"<< childchild10->phi()<< "pt of a"<< p_of_a <<std::endl;
            std::cout<<"eta" << childchild11->eta() << "phi"<< childchild11->phi() <<std::endl;
          }  
          //mass=sqrt((child1->energy()+child0->energy())*(child1->energy()+child0->energy())
	  //	-(child1->px()+child0->px())*(child1->px()+child0->px())
	  //	-(child1->py()+child0->py())*(child1->py()+child0->py())
	  //	-(child1->pz()+child0->pz())*(child1->pz()+child0->pz()));
          //reco::LeafCandidate::LorentzVector sum;
          //sum=child0->p4() + child1->p4();
          histos1D_[ "test" ]->Fill(child0->mass());
          histos1D_[ "DeltaR" ]->Fill(DR_tau);
          histos1D_[ "DeltaR_mu" ]->Fill(DR_mu);
          histos2D_[ "DeltaR_a_and_tau_VS_p_of_a" ]->Fill(DR_a_and_tau, p_of_a);
          DeltaR_a_and_tau_VS_pt_of_a_->Fill(DR_a_and_tau, p_of_a);
        }
       if(fabs(childchild10->pdgId())==13 &&( fabs(childchild00->pdgId())==15))
        {
         // std::cout << "\tthen mom pdgid=" << child0->pdgId() << "\tchildchild00" << childchild00->pdgId() << "\tchildchild01" << childchild01->pdgId() <<std::endl;
         // std::cout << "\tAnd mom pdgid=" << child1->pdgId() << "\tchildchild10" << childchild10->pdgId() << "\tchildchild11" << childchild11->pdgId() <<std::endl;
         // std::cout << "\tStatus of Tau" << childchild00->status() << std::endl;
         // std::cout << "\tStatus of Mu" << childchild10->status() << std::endl;
          mass=0.0;
          mass=child0->mass();
          double DR_tau=VariousFunctions::getDiThingDR_1(childchild00, childchild01);
    
          double DR_mu=VariousFunctions::getDiThingDR_1(childchild10, childchild11);
          double DR_a_and_tau=VariousFunctions::getDiThingDR_1(child0, childchild01);
          double p_of_a=child1->pt();
 
          histos1D_[ "test" ]->Fill(child0->mass());
          histos1D_[ "DeltaR" ]->Fill(DR_tau);
          histos1D_[ "DeltaR_mu"] ->Fill(DR_mu);
          histos2D_[ "DeltaR_a_and_tau_VS_p_of_a" ]->Fill(DR_a_and_tau, p_of_a);
          DeltaR_a_and_tau_VS_pt_of_a_->Fill(DR_a_and_tau, p_of_a);
         // std::cout<< "mass is" << mass <<std::endl;
         // std::cout<< "DR_tau" << DR_tau<<std::endl;
        }
      }    
    }
  }  
}


void
AmumuAnalyzer::reset(const bool doDelete)
{
  if ((doDelete) && (DeltaR_a_and_tau_VS_pt_of_a_ != NULL)) delete DeltaR_a_and_tau_VS_pt_of_a_;
  DeltaR_a_and_tau_VS_pt_of_a_ = NULL;
}
//
// member functions
//

// ------------ method called for each event  ------------


// ------------ method called once each job just before starting event loop  ------------
void 
AmumuAnalyzer::beginJob()
{
  out_=new TFile(outFileName_.c_str(), "RECREATE");
  //jetOutput_.open(jetOutputFileName_.c_str());
  edm::Service< TFileService > fileService;
  histos1D_[ "test" ]=fileService->make< TH1D >("test", "invariant mass of a from DiMuon", 30, 0, 10);
  histos1D_[ "DeltaR" ]=fileService->make<TH1D>("DeltaR", "DeltaR", 100, 0.0, 5.0);
  histos1D_[ "DeltaR_mu"]=fileService->make<TH1D>("DeltaR_mu","DeltaR_mu", 100, 0, 5.0);
  histos2D_[ "DeltaR_a_and_tau_VS_p_of_a" ]=fileService->make<TH2D>("DeltaR_a_and_tau_VS_p_of_a", "DeltaR_a_and_tau_VS_p_of_a", 100, 0, 1.0, 100, 0, 120); 
  DeltaR_a_and_tau_VS_pt_of_a_ = new TH2F("DeltaR_a_and_tau_VS_pt_of_a", "", 100, 0, 1.0, 100, 0, 120);
}

// ------------ method called once each job just after ending the event loop  ------------
void 
AmumuAnalyzer::endJob() 
{
 // jetOutput_.close();
  TCanvas DeltaR_a_and_tau_VS_pt_of_a_Canvas("CanvasName","",600,600);
  VariousFunctions::formatAndDrawCanvasAndHist2D(DeltaR_a_and_tau_VS_pt_of_a_Canvas, DeltaR_a_and_tau_VS_pt_of_a_, 0, 0, 0, kBlack, 7, 20, "Delta_R_of_a_and_tau", .04, .04, 1.1, "pt of a", .04, .04, 1.6, "", .04, .04, 1.0);
  out_->cd();
  DeltaR_a_and_tau_VS_pt_of_a_Canvas.Write();
  out_->Write();
  out_->Close();
}

// ------------ method called when starting to processes a run  ------------
/*
void 
AmumuAnalyzer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a run  ------------
/*
void 
AmumuAnalyzer::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void 
AmumuAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void 
AmumuAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
AmumuAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(AmumuAnalyzer);
