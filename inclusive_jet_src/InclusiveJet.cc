//module for producing a TTree with jet information for doing jet validation studies
// for questions/bugs please contact Virginia Bailey vbailey13@gsu.edu
#include <fun4all/Fun4AllBase.h>
#include <InclusiveJet.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/PHTFileServer.h>
#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>
#include <jetbase/JetMap.h>
#include <jetbase/JetContainer.h>
#include <jetbase/Jetv2.h>
#include <jetbase/Jetv1.h>
#include <centrality/CentralityInfo.h>
#include <globalvertex/GlobalVertex.h>
#include <globalvertex/GlobalVertexMap.h>
#include <calobase/RawTower.h>
#include <calobase/RawTowerContainer.h>
#include <calobase/RawTowerGeom.h>
#include <calobase/RawTowerGeomContainer.h>
#include <calobase/TowerInfoContainer.h>
#include <calobase/TowerInfo.h>
#include <ffarawobjects/Gl1Packet.h>
#include <jetbackground/TowerBackground.h>
#include <calobase/RawTowerDefs.h>
#include <g4main/PHG4TruthInfoContainer.h>
#include <g4main/PHG4Particle.h>
#include <ffaobjects/EventHeaderv1.h>

#include <calobase/RawCluster.h>
#include <calobase/RawClusterContainer.h>
#include <calobase/RawClusterUtility.h>
#include <calobase/RawTowerGeomContainer.h>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#include <HepMC/GenEvent.h>
#pragma GCC diagnostic pop
#include <HepMC/GenVertex.h>
#include <phhepmc/PHHepMCGenEvent.h>
#include <phhepmc/PHHepMCGenEventMap.h>

#include <TTree.h>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <vector>

#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include <TTree.h>
#include <TVector3.h>

//____________________________________________________________________________..
InclusiveJet::InclusiveJet(const std::string& recojetname, const std::string& truthjetname, const std::string& outputfilename):
  SubsysReco("InclusiveJet_" + recojetname + "_" + truthjetname)
  , m_recoJetName(recojetname)
  , m_truthJetName(truthjetname)
  , m_outputFileName(outputfilename)
  , m_etaRange(-1, 1)
  , m_ptRange(5, 100)
  , m_doTruthJets(0)
  , m_doTruth(0)
  , m_doTowers(0)
  , m_doTopoclusters(0)
  , m_doSeeds(0)
  , m_T(nullptr)
  , m_event(-1)
  , m_nTruthJet(-1)
  , m_nJet(-1)
  , m_id()
  , m_nComponent()
  , m_eta()
  , m_phi()
  , m_e()
  , m_pt()
  , m_truthID()
  , m_truthNComponent()
  , m_truthEta()
  , m_truthPhi()
  , m_truthE()
  , m_truthPt()
  , m_eta_rawseed()
  , m_phi_rawseed()
  , m_pt_rawseed()
  , m_e_rawseed()
  , m_rawseed_cut()
  , m_eta_subseed()
  , m_phi_subseed()
  , m_pt_subseed()
  , m_e_subseed()
  , m_subseed_cut()
{
  std::cout << "InclusiveJet::InclusiveJet(const std::string &name) Calling ctor" << std::endl;
}

//____________________________________________________________________________..
InclusiveJet::~InclusiveJet()
{
  std::cout << "InclusiveJet::~InclusiveJet() Calling dtor" << std::endl;
}

//____________________________________________________________________________..
//____________________________________________________________________________..
int InclusiveJet::Init(PHCompositeNode *topNode)
{
  std::cout << "InclusiveJet::Init(PHCompositeNode *topNode) Initializing" << std::endl;
  PHTFileServer::get().open(m_outputFileName, "RECREATE");
  std::cout << "InclusiveJet::Init - Output to " << m_outputFileName << std::endl;

  // configure Tree
  m_T = new TTree("T", "MyJetAnalysis Tree");
  m_T->Branch("m_event", &m_event, "event/I");
  m_T->Branch("nJet", &m_nJet, "nJet/I");
  m_T->Branch("cent", &m_centrality);
  m_T->Branch("zvtx", &m_zvtx);
  m_T->Branch("b", &m_impactparam);
  m_T->Branch("id", &m_id);
  m_T->Branch("nComponent", &m_nComponent);
  m_T->Branch("triggerVector", &m_triggerVector);

  m_T->Branch("eta", &m_eta);
  m_T->Branch("phi", &m_phi);
  m_T->Branch("e", &m_e);
  m_T->Branch("pt", &m_pt);

  if(m_doTruthJets){
    m_T->Branch("nTruthJet", &m_nTruthJet);
    m_T->Branch("truthID", &m_truthID);
    m_T->Branch("truthNComponent", &m_truthNComponent);
    m_T->Branch("truthEta", &m_truthEta);
    m_T->Branch("truthPhi", &m_truthPhi);
    m_T->Branch("truthE", &m_truthE);
    m_T->Branch("truthPt", &m_truthPt);
  }

  if(m_doSeeds){
    m_T->Branch("rawseedEta", &m_eta_rawseed);
    m_T->Branch("rawseedPhi", &m_phi_rawseed);
    m_T->Branch("rawseedPt", &m_pt_rawseed);
    m_T->Branch("rawseedE", &m_e_rawseed);
    m_T->Branch("rawseedCut", &m_rawseed_cut);
    m_T->Branch("subseedEta", &m_eta_subseed);
    m_T->Branch("subseedPhi", &m_phi_subseed);
    m_T->Branch("subseedPt", &m_pt_subseed);
    m_T->Branch("subseedE", &m_e_subseed);
    m_T->Branch("subseedCut", &m_subseed_cut);
  }

  if(m_doTowers) {
    m_T->Branch("emcaln",&m_emcaln,"emcaln/I");
    m_T->Branch("emcale",m_emcale,"emcale[emcaln]/F");
    m_T->Branch("emcalchi2", m_emcalchi2,"emcalchi2[emcaln]/F");
    m_T->Branch("emcaleta",m_emcaleta,"emcaleta[emcaln]/F");
    m_T->Branch("emcalphi",m_emcalphi,"emcalphi[emcaln]/F");
    m_T->Branch("emcalieta",m_emcalieta,"emcalieta[emcaln]/I");
    m_T->Branch("emcaliphi",m_emcaliphi,"emcaliphi[emcaln]/I");
    m_T->Branch("emcalstatus",m_emcalstatus,"emcalstatus[emcaln]/I");
    m_T->Branch("emcaltime",m_emcaltime,"emcaltime[emcaln]/F");

    m_T->Branch("ihcaln",&m_ihcaln,"ihcaln/I");
    m_T->Branch("ihcale",m_ihcale,"ihcale[ihcaln]/F");
    m_T->Branch("ihcalchi2", m_ihcalchi2,"ihcalchi2[ihcaln]/F");
    m_T->Branch("ihcaleta",m_ihcaleta,"ihcaleta[ihcaln]/F");
    m_T->Branch("ihcalphi",m_ihcalphi,"ihcalphi[ihcaln]/F");
    m_T->Branch("ihcalieta",m_ihcalieta,"ihcalieta[ihcaln]/I");
    m_T->Branch("ihcaliphi",m_ihcaliphi,"ihcaliphi[ihcaln]/I");
    m_T->Branch("ihcalstatus",m_ihcalstatus,"ihcalstatus[ihcaln]/I");
    m_T->Branch("ihcaltime",m_ihcaltime,"ihcaltime[ihcaln]/F");

    m_T->Branch("ohcaln",&m_ohcaln,"ohcaln/I");
    m_T->Branch("ohcale",m_ohcale,"ohcale[ohcaln]/F");
    m_T->Branch("ohcalchi2", m_ohcalchi2,"ohcalchi2[ohcaln]/F");
    m_T->Branch("ohcaleta",m_ohcaleta,"ohcaleta[ohcaln]/F");
    m_T->Branch("ohcalphi",m_ohcalphi,"ohcalphi[ohcaln]/F");
    m_T->Branch("ohcalieta",m_ohcalieta,"ohcalieta[ohcaln]/I");
    m_T->Branch("ohcaliphi",m_ohcaliphi,"ohcaliphi[ohcaln]/I");
    m_T->Branch("ohcalstatus",m_ohcalstatus,"ohcalstatus[ohcaln]/I");
    m_T->Branch("ohcaltime",m_ohcaltime,"ohcaltime[ohcaln]/F");
  }

  if(m_doTopoclusters) {
    m_T->Branch("clsmult",&m_clsmult,"clsmult/I");
    m_T->Branch("cluster_e",m_cluster_e,"cluster_e[clsmult]/F");
    m_T->Branch("cluster_eta",m_cluster_eta,"cluster_eta[clsmult]/F");
    m_T->Branch("cluster_phi",m_cluster_phi,"cluster_phi[clsmult]/F");

  }

  if(m_doTruth) {
    m_T->Branch("truthpar_n",&truthpar_n,"truthpar_n/I");
    m_T->Branch("truthpar_pz",truthpar_pz,"truthpar_pz[truthpar_n]/F");
    m_T->Branch("truthpar_pt",truthpar_pt,"truthpar_pt[truthpar_n]/F");
    m_T->Branch("truthpar_e",truthpar_e,"truthpar_e[truthpar_n]/F");
    m_T->Branch("truthpar_eta",truthpar_eta,"truthpar_eta[truthpar_n]/F");
    m_T->Branch("truthpar_phi",truthpar_phi,"truthpar_phi[truthpar_n]/F");
    m_T->Branch("truthpar_pid",truthpar_pid,"truthpar_pid[truthpar_n]/I");
  }
  
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int InclusiveJet::InitRun(PHCompositeNode *topNode)
{
  std::cout << "InclusiveJet::InitRun(PHCompositeNode *topNode) Initializing for Run XXX" << std::endl;
  
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int InclusiveJet::process_event(PHCompositeNode *topNode)
{
  //  std::cout << "InclusiveJet::process_event(PHCompositeNode *topNode) Processing Event" << std::endl;
  ++m_event;

  // interface to reco jets
  JetContainer* jets = findNode::getClass<JetContainer>(topNode, m_recoJetName);
  if (!jets)
    {
      std::cout
	<< "MyJetAnalysis::process_event - Error can not find DST Reco JetContainer node "
	<< m_recoJetName << std::endl;
      exit(-1);
    }

  //interface to truth jets
  //JetMap* jetsMC = findNode::getClass<JetMap>(topNode, m_truthJetName);
  JetContainer* jetsMC = findNode::getClass<JetContainer>(topNode, m_truthJetName);
  if (!jetsMC && m_doTruthJets)
    {
      std::cout
	<< "MyJetAnalysis::process_event - Error can not find DST Truth JetMap node "
	<< m_truthJetName << std::endl;
      exit(-1);
    }
  
  // interface to jet seeds
  JetContainer* seedjetsraw = findNode::getClass<JetContainer>(topNode, "AntiKt_TowerInfo_HIRecoSeedsRaw_r02");
  if (!seedjetsraw && m_doSeeds)
    {
      std::cout
	<< "MyJetAnalysis::process_event - Error can not find DST raw seed jets "
<< std::endl;
      exit(-1);
    }

  JetContainer* seedjetssub = findNode::getClass<JetContainer>(topNode, "AntiKt_TowerInfo_HIRecoSeedsSub_r02");
  if (!seedjetssub && m_doSeeds)
    {
      std::cout
<< "MyJetAnalysis::process_event - Error can not find DST subtracted seed jets "
<< std::endl;
      exit(-1);
    }

  //centrality
  CentralityInfo* cent_node = findNode::getClass<CentralityInfo>(topNode, "CentralityInfo");
  if (!cent_node)
    {
      std::cout
        << "MyJetAnalysis::process_event - Error can not find centrality node "
        << std::endl;
      exit(-1);
    }
  
  //zvertex
  GlobalVertexMap *vertexmap = findNode::getClass<GlobalVertexMap>(topNode, "GlobalVertexMap");
  if (!vertexmap)
    {
      std::cout
        << "MyJetAnalysis::process_event - Error can not find global vertex  node "
        << std::endl;
      exit(-1);
    }
  if (vertexmap->empty())
    {
      std::cout
        << "MyJetAnalysis::process_event - global vertex node is empty "
        << std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }
  else
    {
      GlobalVertex *vtx = vertexmap->begin()->second;
      m_zvtx = vtx->get_z();
    }

    if (fabs(m_zvtx) > 30) {
      return Fun4AllReturnCodes::EVENT_OK;
    }

  //calorimeter towers
  TowerInfoContainer *towersEM3 = findNode::getClass<TowerInfoContainer>(topNode, "TOWERINFO_CALIB_CEMC_RETOWER");
  TowerInfoContainer *towersIH3 = findNode::getClass<TowerInfoContainer>(topNode, "TOWERINFO_CALIB_HCALIN");
  TowerInfoContainer *towersOH3 = findNode::getClass<TowerInfoContainer>(topNode, "TOWERINFO_CALIB_HCALOUT");
  RawTowerGeomContainer *tower_geomEM = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_CEMC");
  RawTowerGeomContainer *tower_geom = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_HCALIN");
  RawTowerGeomContainer *tower_geomOH = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_HCALOUT");
  if(!towersEM3 || !towersIH3 || !towersOH3){
    std::cout
      <<"MyJetAnalysis::process_event - Error can not find raw tower node "
      << std::endl;
    exit(-1);
  }

  TowerInfoContainer *EMtowers = findNode::getClass<TowerInfoContainer>(topNode, "TOWERINFO_CALIB_CEMC");
  TowerInfoContainer *IHtowers = findNode::getClass<TowerInfoContainer>(topNode, "TOWERINFO_CALIB_HCALIN");
  TowerInfoContainer *OHtowers = findNode::getClass<TowerInfoContainer>(topNode, "TOWERINFO_CALIB_HCALOUT");

  if(!tower_geom || !tower_geomOH){
    std::cout
      <<"MyJetAnalysis::process_event - Error can not find raw tower geometry "
      << std::endl;
    exit(-1);
  }

  RawClusterContainer *topoclusters = findNode::getClass<RawClusterContainer>(topNode,"TOPOCLUSTER_ALLCALO"); 
  if (m_doTopoclusters && !topoclusters) {
        std::cout
    <<"MyJetAnalysis::process_event - Error can not find topoclusters "
    << std::endl;
  exit(-1);
  }

  //get the event centrality/impact parameter from HIJING
  //m_centrality =  cent_node->get_centile(CentralityInfo::PROP::mbd_NS);
  m_centrality = (int)(100*cent_node->get_centile(CentralityInfo::PROP::mbd_NS));
  m_impactparam =  cent_node->get_quantity(CentralityInfo::PROP::bimp);

  //get reco jets
  m_nJet = 0;

  for (auto jet : *jets)
    {

      if(jet->get_pt() < 5) continue; // to remove noise jets
      if(jet->get_e() < 0) {
        return Fun4AllReturnCodes::EVENT_OK; // currently applied to deal with cold EMCal IB
      }
      if (fabs(jet->get_eta()) > 0.7) continue;

      m_id.push_back(jet->get_id());
      m_nComponent.push_back(jet->size_comp());
      m_eta.push_back(jet->get_eta());
      m_phi.push_back(jet->get_phi());
      m_e.push_back(jet->get_e());
      m_pt.push_back(jet->get_pt());

      float emcalE = 0;
      float ihcalE = 0;
      float ohcalE = 0;
      for (auto comp: jet->get_comp_vec())
      {
        TowerInfo *tower;
        unsigned int channel = comp.second;
        if (comp.first == 14 || comp.first == 29 || comp.first == 28) {
          tower = towersEM3->get_tower_at_channel(channel);
          if(!tower) { continue; }
          emcalE += tower->get_energy();
        }    
        if (comp.first == 15 ||  comp.first == 30 || comp.first == 26)
        {
          tower = towersIH3->get_tower_at_channel(channel);
          if(!tower) { continue; }
          ihcalE += tower->get_energy();
        }

        if (comp.first == 16 ||  comp.first == 31 || comp.first == 27)
        {
          tower = towersOH3->get_tower_at_channel(channel);
          if(!tower) { continue; }
          ohcalE += tower->get_energy();
        }
      }

      if (ohcalE/m_e.back() > 0.9) {
        std::cout << "event " << m_event << " emcalE " << emcalE << " ohcalE " << ohcalE << " totalE " << emcalE + ihcalE + ohcalE << " Jet E " << m_e.back() << std::endl;
        return Fun4AllReturnCodes::EVENT_OK;
      }
      m_nJet++;
    }

  //get truth jets
  if(m_doTruthJets)
    {
      m_nTruthJet = 0;
      //for (JetMap::Iter iter = jetsMC->begin(); iter != jetsMC->end(); ++iter)
      for (auto truthjet : *jetsMC)
	{
	  //Jet* truthjet = iter->second;
	    
	  bool eta_cut = (truthjet->get_eta() >= m_etaRange.first) and (truthjet->get_eta() <= m_etaRange.second);
	  bool pt_cut = (truthjet->get_pt() >= m_ptRange.first) and (truthjet->get_pt() <= m_ptRange.second);
	  if ((not eta_cut) or (not pt_cut)) continue;
	  m_truthID.push_back(truthjet->get_id());
	  m_truthNComponent.push_back(truthjet->size_comp());
	  m_truthEta.push_back(truthjet->get_eta());
	  m_truthPhi.push_back(truthjet->get_phi());
	  m_truthE.push_back(truthjet->get_e());
	  m_truthPt.push_back(truthjet->get_pt());
	  m_nTruthJet++;
	}
    }
  
  //get seed jets
  if(m_doSeeds)
    {
      for (auto jet : *seedjetsraw)
	{
	  int passesCut = jet->get_property(seedjetsraw->property_index(Jet::PROPERTY::prop_SeedItr));
	  m_eta_rawseed.push_back(jet->get_eta());
	  m_phi_rawseed.push_back(jet->get_phi());
	  m_e_rawseed.push_back(jet->get_e());
	  m_pt_rawseed.push_back(jet->get_pt());
	  m_rawseed_cut.push_back(passesCut);
	}
      
      for (auto jet : *seedjetssub)
	{
	  int passesCut = jet->get_property(seedjetssub->property_index(Jet::PROPERTY::prop_SeedItr));
	  m_eta_subseed.push_back(jet->get_eta());
	  m_phi_subseed.push_back(jet->get_phi());
	  m_e_subseed.push_back(jet->get_e());
	  m_pt_subseed.push_back(jet->get_pt());
	  m_subseed_cut.push_back(passesCut);
	}
    }

  //grab the gl1 data
  Gl1Packet *gl1PacketInfo = findNode::getClass<Gl1Packet>(topNode, "GL1Packet");
  if (!gl1PacketInfo)
    {
      std::cout << PHWHERE << "caloTreeGen::process_event: GL1Packet node is missing. Output related to this node will be empty" << std::endl;
    }

  if (gl1PacketInfo) {
    bool jettrig = false;
    uint64_t triggervec = gl1PacketInfo->getScaledVector();
    for (int i = 0; i < 64; i++) {
	    bool trig_decision = ((triggervec & 0x1U) == 0x1U);
      if (trig_decision) {
        m_triggerVector.push_back(i);
        if (i >= 16 && i <= 23) { jettrig = true; }
      }
	    triggervec = (triggervec >> 1U) & 0xffffffffU;
	  }
    if (!jettrig) {
      return Fun4AllReturnCodes::EVENT_OK;
    }
  }

  if (m_doTowers) {
    if(EMtowers) { 
      m_emcaln = 0;
      int nchannels = 24576; 
      for(int i=0; i<nchannels; ++i) {
        TowerInfo *tower = EMtowers->get_tower_at_channel(i); //get EMCal tower
        if(tower->get_isHot() || tower->get_isNoCalib() || tower->get_isNotInstr()) { continue; }
        int key = EMtowers->encode_key(i);
        int etabin = EMtowers->getTowerEtaBin(key);
        int phibin = EMtowers->getTowerPhiBin(key);
        m_emcale[m_emcaln] = tower->get_energy(); 
        m_emcalchi2[m_emcaln] = tower->get_chi2();
        if (!tower->get_isBadChi2() && tower->get_chi2() > 10000) {
          std::cout << "event " << m_event << " ieta " << etabin << " iphi " << phibin << " chi2 " << tower->get_chi2() << " e " << tower->get_energy() << std::endl;
        }
        const RawTowerDefs::keytype geomkey = RawTowerDefs::encode_towerid(RawTowerDefs::CalorimeterId::CEMC, etabin, phibin);
        RawTowerGeom *geom = tower_geomEM->get_tower_geometry(geomkey); 
        TVector3 tower_pos;
        tower_pos.SetXYZ(geom->get_center_x(),geom->get_center_y(),geom->get_center_z() - m_zvtx);
        m_emcaleta[m_emcaln] = tower_pos.Eta();
        m_emcalphi[m_emcaln] = tower_pos.Phi();
        m_emcalieta[m_emcaln] = etabin;
        m_emcaliphi[m_emcaln] = phibin;
        m_emcalstatus[m_emcaln] = tower->get_status();
        m_emcaltime[m_emcaln] = tower->get_time_float();
        m_emcaln++;
      }
    }

    if(IHtowers) { 
      m_ihcaln = 0;
      int nchannels = 1536; 
      for(int i=0; i<nchannels; ++i) {
        TowerInfo *tower = IHtowers->get_tower_at_channel(i); //get IHCal tower
        if(tower->get_isHot() || tower->get_isNoCalib() || tower->get_isNotInstr()) { continue; }
        int key = IHtowers->encode_key(i);
        int etabin = IHtowers->getTowerEtaBin(key);
        int phibin = IHtowers->getTowerPhiBin(key);
        m_ihcale[m_ihcaln] = tower->get_energy(); 
        m_ihcalchi2[m_ihcaln] = tower->get_chi2();
        const RawTowerDefs::keytype geomkey = RawTowerDefs::encode_towerid(RawTowerDefs::CalorimeterId::HCALIN, etabin, phibin);
        RawTowerGeom *geom = tower_geom->get_tower_geometry(geomkey); 
        TVector3 tower_pos;
        tower_pos.SetXYZ(geom->get_center_x(),geom->get_center_y(),geom->get_center_z() - m_zvtx);
        m_ihcaleta[m_ihcaln] = tower_pos.Eta();
        m_ihcalphi[m_ihcaln] = tower_pos.Phi();
        m_ihcalieta[m_ihcaln] = etabin;
        m_ihcaliphi[m_ihcaln] = phibin;
        m_ihcalstatus[m_ihcaln] = tower->get_status();
        m_ihcaltime[m_ihcaln] = tower->get_time_float();
        m_ihcaln++;
      }
    }

    if(OHtowers) { 
      m_ohcaln = 0;
      int nchannels = 1536; 
      for(int i=0; i<nchannels; ++i) {
        TowerInfo *tower = OHtowers->get_tower_at_channel(i); //get OHCal tower
        if(tower->get_isHot() || tower->get_isNoCalib() || tower->get_isNotInstr()) { continue; }
        int key = OHtowers->encode_key(i);
        int etabin = OHtowers->getTowerEtaBin(key);
        int phibin = OHtowers->getTowerPhiBin(key);
        m_ohcale[m_ohcaln] = tower->get_energy(); 
        m_ohcalchi2[m_ohcaln] = tower->get_chi2();
        const RawTowerDefs::keytype geomkey = RawTowerDefs::encode_towerid(RawTowerDefs::CalorimeterId::HCALOUT, etabin, phibin);
        RawTowerGeom *geom = tower_geomOH->get_tower_geometry(geomkey); 
        TVector3 tower_pos;
        tower_pos.SetXYZ(geom->get_center_x(),geom->get_center_y(),geom->get_center_z() - m_zvtx);
        m_ohcaleta[m_ohcaln] = tower_pos.Eta();
        m_ohcalphi[m_ohcaln] = tower_pos.Phi();
        m_ohcalieta[m_ohcaln] = etabin;
        m_ohcaliphi[m_ohcaln] = phibin;
        m_ohcalstatus[m_ohcaln] = tower->get_status();
        m_ohcaltime[m_ohcaln] = tower->get_time_float();
        m_ohcaln++;
      }
    }
  }
  if (m_doTopoclusters) {
    if (topoclusters) {
      RawClusterContainer::Map clusterMap = topoclusters->getClustersMap();
      m_clsmult = 0;
      for(auto entry : clusterMap){
        RawCluster* cluster = entry.second;
        CLHEP::Hep3Vector origin(0, 0, m_zvtx);
        m_cluster_e[m_clsmult] = cluster->get_energy();
        m_cluster_eta[m_clsmult] = RawClusterUtility::GetPseudorapidity(*cluster, origin);
        m_cluster_phi[m_clsmult] = RawClusterUtility::GetAzimuthAngle(*cluster, origin);
        m_clsmult++;
        if (m_clsmult == 10000) { break; }
      } 
    }
  }

  if (m_doTruth) {

    PHHepMCGenEventMap *hepmceventmap = findNode::getClass<PHHepMCGenEventMap>(topNode, "PHHepMCGenEventMap");
    if (!hepmceventmap) {
      std::cout << PHWHERE << "HEPMC event map node is missing, can't collected HEPMC truth particles"<< std::endl;
      return Fun4AllReturnCodes::EVENT_OK;
    }
    for (PHHepMCGenEventMap::ConstIter eventIter = hepmceventmap->begin(); eventIter != hepmceventmap->end(); ++eventIter) {
    /// Get the event
      PHHepMCGenEvent *hepmcevent = eventIter->second; 
      // To fill TTree, require that the event be the primary event (embedding_id > 0)
      if (hepmcevent && hepmcevent->get_embedding_id() == 1) {
        /// Get the event characteristics, inherited from HepMC classes
        HepMC::GenEvent *truthevent = hepmcevent->getEvent();
        if (!truthevent) {
          std::cout << PHWHERE << "no evt pointer under phhepmvgeneventmap found " << std::endl;
          return Fun4AllReturnCodes::EVENT_OK;
        }
        //int process_id = truthevent->signal_process_id();
        //std::cout << process_id << std::endl;
        /// Loop over all the truth particles and get their information
        truthpar_n = 0;
        for (HepMC::GenEvent::particle_const_iterator iter = truthevent->particles_begin(); iter != truthevent->particles_end(); ++iter) {
          if (!(*iter)->end_vertex() && (*iter)->status() == 1) {
            truthpar_e[truthpar_n] = (*iter)->momentum().e();
            double px = (*iter)->momentum().px();
            double py = (*iter)->momentum().py();
            double pz = (*iter)->momentum().pz();
            truthpar_pt[truthpar_n] = sqrt(px*px+py*py);
            truthpar_pz[truthpar_n] = pz;
            truthpar_phi[truthpar_n]=atan2(py,px);
            truthpar_eta[truthpar_n]=atanh(pz/sqrt(px*px+py*py+pz*pz));
            truthpar_pid[truthpar_n]=(*iter)->pdg_id();
            /// Fill the truth tree
            truthpar_n++;
          }
        }
        break;
      } 
    }
  }
  
  //fill the tree
  m_T->Fill();

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int InclusiveJet::ResetEvent(PHCompositeNode *topNode)
{
  //std::cout << "InclusiveJet::ResetEvent(PHCompositeNode *topNode) Resetting internal structures, prepare for next event" << std::endl;
  m_id.clear();
  m_nComponent.clear();
  m_eta.clear();
  m_phi.clear();
  m_e.clear();
  m_pt.clear();

  m_truthID.clear();
  m_truthNComponent.clear();
  m_truthEta.clear();
  m_truthPhi.clear();
  m_truthE.clear();
  m_truthPt.clear();
  m_truthdR.clear();

  m_eta_subseed.clear();
  m_phi_subseed.clear();
  m_e_subseed.clear();
  m_pt_subseed.clear();
  m_subseed_cut.clear();

  m_eta_rawseed.clear();
  m_phi_rawseed.clear();
  m_e_rawseed.clear();
  m_pt_rawseed.clear();
  m_rawseed_cut.clear();
  
  m_triggerVector.clear();

  m_emcaln = 0;
  m_ihcaln = 0;
  m_ohcaln = 0;
  m_clsmult = 0;

  for (int i = 0; i < 24576; i++) {
    m_emcale[i] = 0;
    m_emcalieta[i] = 0;
    m_emcaliphi[i] = 0;
    m_emcaleta[i] = 0;
    m_emcalphi[i] = 0;
    m_emcalchi2[i] = 0;
    m_emcalstatus[i] = 0;
    m_emcaltime[i] = 0;
  }

  for (int i = 0; i < 1536; i++) {
    m_ihcale[i] = 0;
    m_ihcalieta[i] = 0;
    m_ihcaliphi[i] = 0;
    m_ihcaleta[i] = 0;
    m_ihcalphi[i] = 0;
    m_ihcalchi2[i] = 0;
    m_ihcalstatus[i] = 0;
    m_ihcaltime[i] = 0;
    m_ohcale[i] = 0;
    m_ohcalieta[i] = 0;
    m_ohcaliphi[i] = 0;
    m_ohcaleta[i] = 0;
    m_ohcalphi[i] = 0;
    m_ohcalchi2[i] = 0;
    m_ohcalstatus[i] = 0;
    m_ohcaltime[i] = 0;
  }

  for (int i = 0; i < 100000; i++) {
    truthpar_pt[i] = 0;
    truthpar_pz[i] = 0;
    truthpar_e[i] = 0;
    truthpar_phi[i] = 0;
    truthpar_eta[i] = 0;
    truthpar_pid[i] = 0;
  }

  for (int i = 0; i < 10000; i++) {
    m_cluster_e[i] = 0;
    m_cluster_eta[i] = 0;
    m_cluster_phi[i] = 0;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int InclusiveJet::EndRun(const int runnumber)
{
  std::cout << "InclusiveJet::EndRun(const int runnumber) Ending Run for Run " << runnumber << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int InclusiveJet::End(PHCompositeNode *topNode)
{
  std::cout << "InclusiveJet::End - Output to " << m_outputFileName << std::endl;
  PHTFileServer::get().cd(m_outputFileName);

  m_T->Write();
  std::cout << "InclusiveJet::End(PHCompositeNode *topNode) This is the End..." << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int InclusiveJet::Reset(PHCompositeNode *topNode)
{
  std::cout << "InclusiveJet::Reset(PHCompositeNode *topNode) being Reset" << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
void InclusiveJet::Print(const std::string &what) const
{
  std::cout << "InclusiveJet::Print(const std::string &what) const Printing info for " << what << std::endl;
}