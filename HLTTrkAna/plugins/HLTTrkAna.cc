// -*- C++ -*-
//
// Package:    UserCode/HLTTrkAna
// Class:      HLTTrkAna
// 
/**\class HLTTrkAna HLTTrkAna.cc UserCode/HLTTrkAna/plugins/HLTTrkAna.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  dylan rankin
//         Created:  Fri, 31 Mar 2017 04:30:28 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "DataFormats/Math/interface/deltaR.h"

#include "DataFormats/Common/interface/AssociationMap.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "DataFormats/RecoCandidate/interface/RecoEcalCandidate.h"

#include "DataFormats/CaloRecHit/interface/CaloCluster.h"
#include "DataFormats/EgammaReco/interface/ElectronSeed.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"

#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"

#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

#include "HLTrigger/HLTcore/interface/TriggerSummaryAnalyzerAOD.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "HLTrigger/HLTcore/interface/HLTPrescaleProvider.h"

#include "TTree.h"
#include "TLorentzVector.h"

#include <tuple>

namespace edm {
  class ConfigurationDescriptions;
}

//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

//class HLTTrkAna : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
class HLTTrkAna : public edm::stream::EDAnalyzer< >  {
   public:
      explicit HLTTrkAna(const edm::ParameterSet&);
      ~HLTTrkAna();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      double detEtaFromEvnt(double evntEta, double z0);
      double etaToTheta(double eta);
      double thetaToEta(double theta);

      virtual void beginRun(edm::Run const &, edm::EventSetup const&) override;
      virtual void endRun(edm::Run const &, edm::EventSetup const&) override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;

      const edm::EDGetTokenT<edm::TriggerResults> trigResToken_;
      const edm::InputTag inputTag_;
      const edm::EDGetTokenT<trigger::TriggerEvent> inputToken_;
      HLTConfigProvider hltConfig_;

      std::string processName_;

      //edm::EDGetTokenT<std::vector<reco::GsfElectron>> electronToken_;
      edm::EDGetTokenT<std::vector<reco::GenParticle>> genToken_;
      edm::EDGetTokenT<std::vector<PileupSummaryInfo>> puToken_;
      std::string triggerName_;
      std::vector<std::string> filterNames0_;
      std::vector<std::string> filterNames1_;
      std::vector<std::string> filterNames2_;

      TTree *Tree0;
      TTree *Tree1;
      TTree *Tree2;

      std::vector<float> genpt_;
      std::vector<float> geneta_;
      std::vector<float> genphi_;
      std::vector<float> genet_;
      std::vector<float> hltpt_;
      std::vector<float> hlteta_;
      std::vector<float> hltphi_;
      std::vector<float> hltet_;
      std::vector<float> vz_;
      std::vector<int> nPU_;

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
HLTTrkAna::HLTTrkAna(const edm::ParameterSet& iConfig):
   trigResToken_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("trigresults"))),
   inputTag_(iConfig.getParameter<edm::InputTag>("triggers")),
   inputToken_(consumes<trigger::TriggerEvent>(inputTag_)),
   processName_(iConfig.getParameter<std::string>("processname")),
   genToken_(consumes<std::vector<reco::GenParticle>>(iConfig.getParameter<edm::InputTag>("genparticles"))),
   puToken_(consumes<std::vector<PileupSummaryInfo>>(iConfig.getParameter<edm::InputTag>("pusummary"))),
   triggerName_(iConfig.getParameter<std::string>("triggername")),
   filterNames0_(iConfig.getParameter<std::vector<std::string>>("filternames0")),
   filterNames1_(iConfig.getParameter<std::vector<std::string>>("filternames1")),
   filterNames2_(iConfig.getParameter<std::vector<std::string>>("filternames2"))
{
   //now do what ever initialization is needed
   edm::Service<TFileService> fs;

   Tree0 = fs->make<TTree>( "Tree0", "Tree0");
   Tree1 = fs->make<TTree>( "Tree1", "Tree1");
   Tree2 = fs->make<TTree>( "Tree2", "Tree2");

   Tree0->Branch("nPU", &nPU_);
   Tree0->Branch("genpt", &genpt_);
   Tree0->Branch("geneta", &geneta_);
   Tree0->Branch("genphi", &genphi_);
   Tree0->Branch("genet", &genet_);
   Tree0->Branch("hltpt", &hltpt_);
   Tree0->Branch("hlteta", &hlteta_);
   Tree0->Branch("hltphi", &hltphi_);
   Tree0->Branch("hltet", &hltet_);
   Tree0->Branch("vz", &vz_);

   Tree1->Branch("nPU", &nPU_);
   Tree1->Branch("genpt", &genpt_);
   Tree1->Branch("geneta", &geneta_);
   Tree1->Branch("genphi", &genphi_);
   Tree1->Branch("genet", &genet_);
   Tree1->Branch("hltpt", &hltpt_);
   Tree1->Branch("hlteta", &hlteta_);
   Tree1->Branch("hltphi", &hltphi_);
   Tree1->Branch("hltet", &hltet_);
   Tree1->Branch("vz", &vz_);

   Tree2->Branch("nPU", &nPU_);
   Tree2->Branch("genpt", &genpt_);
   Tree2->Branch("geneta", &geneta_);
   Tree2->Branch("genphi", &genphi_);
   Tree2->Branch("genet", &genet_);
   Tree2->Branch("hltpt", &hltpt_);
   Tree2->Branch("hlteta", &hlteta_);
   Tree2->Branch("hltphi", &hltphi_);
   Tree2->Branch("hltet", &hltet_);
   Tree2->Branch("vz", &vz_);

}


HLTTrkAna::~HLTTrkAna()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

void
HLTTrkAna::endRun(edm::Run const & iRun, edm::EventSetup const& iSetup) {}

void
HLTTrkAna::beginRun(edm::Run const & iRun, edm::EventSetup const& iSetup)
{
   using namespace std;
   using namespace edm;

   bool changed(true);
   if (hltConfig_.init(iRun,iSetup,processName_,changed)) {
     //std::cout<<"Got hltConfig_"<<std::endl;
     hltConfig_.dump("Triggers");
   } else {
     std::cout << "HLTTrkAna::analyze:"
 	 << " config extraction failure with process name "
 	 << processName_ << std::endl;
   }
}

// ------------ method called for each event  ------------
void
HLTTrkAna::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace std;
   using namespace edm;
   using namespace reco;
   using namespace trigger;

   edm::Handle<std::vector<PileupSummaryInfo>> pusumm;
   iEvent.getByToken(puToken_, pusumm);

   int npu = -1;

   std::vector<PileupSummaryInfo>::const_iterator PVI;                       
   for (const PileupSummaryInfo &psi : *pusumm) {
       if (psi.getBunchCrossing() == 0){//Only care about in time PU for now 
           npu = psi.getPU_NumInteractions();           
           break;
       }
   }//Looping over different Bunch Crossings

   std::vector<TLorentzVector> hlt0p4;
   std::vector<TLorentzVector> hlt1p4;
   std::vector<TLorentzVector> hlt2p4;
   TLorentzVector tmpvec;

   //std::cout<<inputTag_.encode()<<std::endl;
   edm::Handle<edm::TriggerResults> trigreshand_;
   iEvent.getByToken(trigResToken_,trigreshand_);

   assert(trigreshand_->size()==hltConfig_.size());

   const unsigned int triggerIndex(hltConfig_.triggerIndex(triggerName_));
   const unsigned int m(hltConfig_.size(triggerIndex));
   const vector<string>& moduleLabels(hltConfig_.moduleLabels(triggerIndex));
   const unsigned int moduleIndex(trigreshand_->index(triggerIndex));
   int runmodule = -1;
   for (unsigned int mm = 0; mm<=moduleIndex; mm++) {
     if (runmodule==2 or runmodule==-2) break;
     if (runmodule==-1) {
       for (unsigned int fi = 0; fi < filterNames0_.size(); fi++) {
         if (moduleLabels[mm]==filterNames0_[fi].substr(0,filterNames0_[fi].find(":"))) {
           runmodule = 0;
           if (mm==moduleIndex) runmodule = -2;
           break;
         }
       }
     }
     if (runmodule==0) {
       for (unsigned int fi = 0; fi < filterNames1_.size(); fi++) {
         if (moduleLabels[mm]==filterNames1_[fi].substr(0,filterNames1_[fi].find(":"))) {
           runmodule = 1;
           break;
         }
       }
     }
     if (runmodule==1) {
       for (unsigned int fi = 0; fi < filterNames2_.size(); fi++) {
         if (moduleLabels[mm]==filterNames2_[fi].substr(0,filterNames2_[fi].find(":"))) {
           runmodule = 2;
           break;
         }
       }
     }
   }
   /*if (moduleLabels[moduleIndex]==filterName0_.substr(0,filterName0_.find(":"))) lastmodule = 0;
   else {
     if (moduleLabels[moduleIndex]==filterName1_.substr(0,filterName1_.find(":"))) lastmodule = 1;
     else {
       if (moduleLabels[moduleIndex]==filterName2_.substr(0,filterName2_.find(":"))) lastmodule = 2;
     }
   }*/

   assert (moduleIndex<m);

   edm::Handle<TriggerEvent> handle;
   iEvent.getByToken(inputToken_,handle);
   if (handle.isValid()) {
     const size_type nF(handle->sizeFilters());
     for (size_type iF=0; iF!=nF; ++iF) {
       bool filterMatch0 = false;
       bool filterMatch1 = false;
       bool filterMatch2 = false;
       //std::cout<<handle->filterTag(iF).encode()<<std::endl;
       for (unsigned int fi = 0; fi < filterNames0_.size(); fi++) {
         if (runmodule>=0 && handle->filterTag(iF).encode()==filterNames0_[fi]) {
           filterMatch0 = true;
           break;
         }
       }
       for (unsigned int fi = 0; fi < filterNames1_.size(); fi++) {
         if (runmodule>=1 && handle->filterTag(iF).encode()==filterNames1_[fi]) {
           filterMatch1 = true;
         }
       }
       for (unsigned int fi = 0; fi < filterNames2_.size(); fi++) {
         if (runmodule>=2 && handle->filterTag(iF).encode()==filterNames2_[fi]) {
           filterMatch2 = true;
         }
       }
       if (filterMatch0 || filterMatch1 || filterMatch2) {}
       else continue;
       const Vids& VIDS (handle->filterIds(iF));
       const Keys& KEYS(handle->filterKeys(iF));
       const size_type nI(VIDS.size());
       const size_type nK(KEYS.size());
       const TriggerObjectCollection& TOC(handle->getObjects());
       const size_type n(max(nI,nK));
       //std::cout<<"Filter "<<handle->filterTag(iF).encode()<<std::endl;
       for (size_type i=0; i!=n; ++i) {
         const TriggerObject& TO(TOC[KEYS[i]]);
         //std::cout<<"\tPt= "<<TO.pt()<<" , Eta= "<<TO.eta()<<" , Phi= "<<TO.phi()<<std::endl;
         tmpvec.SetPtEtaPhiM(TO.pt(),TO.eta(),TO.phi(),TO.mass());
         if (filterMatch0) {
           bool newcand = true;
           for (unsigned int ic = 0; ic<hlt0p4.size(); ic++) {
             if (tmpvec.DeltaR(hlt0p4[ic])<0.01) {newcand = false; break;}
           }
           if (newcand) hlt0p4.push_back(tmpvec);
         }
         if (filterMatch1) {
           bool newcand = true;
           for (unsigned int ic = 0; ic<hlt1p4.size(); ic++) {
             if (tmpvec.DeltaR(hlt1p4[ic])<0.01) {newcand = false; break;}
           }
           if (newcand) hlt1p4.push_back(tmpvec);
         }
         if (filterMatch2) {
           bool newcand = true;
           for (unsigned int ic = 0; ic<hlt2p4.size(); ic++) {
             if (tmpvec.DeltaR(hlt2p4[ic])<0.01) {newcand = false; break;}
           }
           if (newcand) hlt2p4.push_back(tmpvec);
         }
       }
     }
   } else {
     LogVerbatim("TriggerSummaryAnalyzerAOD") << "Handle invalid! Check InputTag provided." << endl;
   }

   if (hlt0p4.size()>0 || hlt1p4.size()>0 || hlt2p4.size()>0) {}
   else return;
   if ( hlt0p4.size()<hlt1p4.size() || hlt0p4.size()<hlt2p4.size() || hlt1p4.size()<hlt2p4.size() ) std::cout<<"Issue with size of trigger object vectors..."<<std::endl;

   edm::Handle<std::vector<reco::GenParticle>> genparts;
   iEvent.getByToken(genToken_, genparts);

   std::vector<std::pair<TLorentzVector,float>> genvec;

   for (const reco::GenParticle &g : *genparts) {
       if (abs(g.pdgId())!=11) continue;
       //std::cout<<"\nGenPart ID"<<g.pdgId()<<" Stat"<<g.status()<<" Pt"<<g.pt()<<" Eta"<<g.eta()<<" Phi"<<g.phi();
       //if (g.mother()) std::cout<<" IDMom"<<g.mother()->pdgId()<<" StatMom"<<g.mother()->status();
       //std::cout<<std::endl;
       if (g.pt() < 5 or g.status()!=23) continue;
       //tmpvec.SetPtEtaPhiE(g.pt(),g.eta(),g.phi(),g.energy());
       tmpvec.SetPtEtaPhiE(g.pt(),HLTTrkAna::detEtaFromEvnt(g.eta(),g.vz()),g.phi(),g.energy());
       std::pair<TLorentzVector,float> tmppair;
       tmppair.first = tmpvec;
       tmppair.second = g.vz();
       genvec.push_back(tmppair);
   }


   nPU_.clear();
   genpt_.clear();
   geneta_.clear();
   genphi_.clear();
   genet_.clear();
   hltpt_.clear();
   hlteta_.clear();
   hltphi_.clear();
   hltet_.clear();
   vz_.clear();
   for (unsigned int i=0; i<hlt0p4.size(); i++) {
       for (unsigned int j=0; j<genvec.size(); j++) {
           if (hlt0p4[i].DeltaR(genvec[j].first)<0.3 && hlt0p4[i].Pt()/genvec[j].first.Pt() > 0.5 && hlt0p4[i].Pt()/genvec[j].first.Pt() < 2.0 ) {
               nPU_.push_back(npu);
               genpt_.push_back(genvec[j].first.Pt());
               geneta_.push_back(genvec[j].first.Eta());
               genphi_.push_back(genvec[j].first.Phi());
               genet_.push_back(genvec[j].first.Et());
               hltpt_.push_back(hlt0p4[i].Pt());
               hlteta_.push_back(hlt0p4[i].Eta());
               hltphi_.push_back(hlt0p4[i].Phi());
               hltet_.push_back(hlt0p4[i].Et());
               vz_.push_back(genvec[j].second);
               break;
           }
       }
   }
   Tree0->Fill();

   nPU_.clear();
   genpt_.clear();
   geneta_.clear();
   genphi_.clear();
   genet_.clear();
   hltpt_.clear();
   hlteta_.clear();
   hltphi_.clear();
   hltet_.clear();
   vz_.clear();
   for (unsigned int i=0; i<hlt1p4.size(); i++) {
       for (unsigned int j=0; j<genvec.size(); j++) {
           if (hlt1p4[i].DeltaR(genvec[j].first)<0.3 && hlt1p4[i].Pt()/genvec[j].first.Pt() > 0.5 && hlt1p4[i].Pt()/genvec[j].first.Pt() < 2.0 ) {
               nPU_.push_back(npu);
               genpt_.push_back(genvec[j].first.Pt());
               geneta_.push_back(genvec[j].first.Eta());
               genphi_.push_back(genvec[j].first.Phi());
               genet_.push_back(genvec[j].first.Et());
               hltpt_.push_back(hlt1p4[i].Pt());
               hlteta_.push_back(hlt1p4[i].Eta());
               hltphi_.push_back(hlt1p4[i].Phi());
               hltet_.push_back(hlt1p4[i].Et());
               vz_.push_back(genvec[j].second);
               break;
           }
       }
   }
   Tree1->Fill();

   nPU_.clear();
   genpt_.clear();
   geneta_.clear();
   genphi_.clear();
   genet_.clear();
   hltpt_.clear();
   hlteta_.clear();
   hltphi_.clear();
   hltet_.clear();
   vz_.clear();
   for (unsigned int i=0; i<hlt2p4.size(); i++) {
       for (unsigned int j=0; j<genvec.size(); j++) {
           if (hlt2p4[i].DeltaR(genvec[j].first)<0.3 && hlt2p4[i].Pt()/genvec[j].first.Pt() > 0.5 && hlt2p4[i].Pt()/genvec[j].first.Pt() < 2.0 ) {
               nPU_.push_back(npu);
               genpt_.push_back(genvec[j].first.Pt());
               geneta_.push_back(genvec[j].first.Eta());
               genphi_.push_back(genvec[j].first.Phi());
               genet_.push_back(genvec[j].first.Et());
               hltpt_.push_back(hlt2p4[i].Pt());
               hlteta_.push_back(hlt2p4[i].Eta());
               hltphi_.push_back(hlt2p4[i].Phi());
               hltet_.push_back(hlt2p4[i].Et());
               vz_.push_back(genvec[j].second);
               break;
           }
       }
   }
   Tree2->Fill();

}

double HLTTrkAna::etaToTheta(double eta)
{
  //  if(eta<0) return -2*atan(exp(eta));
  //  else return 2*atan(exp(-1*eta));
  return 2*atan(exp(-1*eta));
  //else return 2*atan(exp(-1*eta));

}

double HLTTrkAna::thetaToEta(double theta)
{
  //first bounds check theta to get into -pi/2 - pi/2 range
  while( fabs(theta) > asin(1.)){
    if(theta>0) theta-=2*asin(1.);
    else theta+=2*asin(1.);
  }
  //now check sign
  if(theta<0) return log(tan(fabs(theta/2.)));
  else return -1.*log(tan(theta/2.));
}

double HLTTrkAna::detEtaFromEvnt(double evntEta,double z0)
{
 
  double thetaEvt = HLTTrkAna::etaToTheta(evntEta);
  double z = 129.4 / tan(thetaEvt); //129.4 is the average barrel radius
  double zTot = z+z0;

  if(fabs(zTot)<269){ //269 is an emperically derived number which means that < its likely in th barrel
    return zTot !=0 ? HLTTrkAna::thetaToEta(atan(129.4/zTot)) : 0.; //otherwise endcap time
  }
  double endcapZ = 319.2; //average z position of endcap
  if(evntEta<0) endcapZ*=-1;
  double rxy = tan(thetaEvt) * (endcapZ-z0);
  return HLTTrkAna::thetaToEta(atan(rxy/endcapZ));

}


// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
HLTTrkAna::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(HLTTrkAna);
