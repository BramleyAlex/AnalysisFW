// Fill and Fill tree functions to produce BDT training set
#include <vector>
#include <algorithm>
#include <iostream>
#include "CLoop.h"
#include "OutputTree.h"

double del_phi(double phi_1, double phi_2);
double min_deltaR(const TLorentzVector& test_particle, const TLorentzVector& jet1, const TLorentzVector& jet2);
TLorentzVector& toGeV(TLorentzVector &v);
double CalculatePtBalance(const std::vector<TLorentzVector>& particles);
double CalculateOmega(const TLorentzVector& tau_0_p4, const TLorentzVector& tau_1_p4, const TLorentzVector& met_p4);
std::pair<TLorentzVector,TLorentzVector> CalculateNeutrinoVector(const TLorentzVector& met_p4, const TLorentzVector& tau_0_p4, 
    const TLorentzVector& tau_1_p4);
void PrintTLorentzVectorPtEtaPhi(const TLorentzVector& input);
void PrintTLorentzVectorXYZT(const TLorentzVector& input);
bool Region(std::string region_id, const bool& N_gap_jets, const double& Z_centrality);
bool CalculateNGapJets(const double &ljet_0_rapidity, const double &ljet_1_rapidity, const std::vector<float> *JetEta);
bool PassTauID(const std::string& working_point, const double& RNN_score, const int& n_tracks);

void BDT_Sample::Fill(double weight, int z_sample, const std::string& sampleName) {
  double pi=TMath::Pi();

  // Jet vectors
  TLorentzVector ljet_0_p4;
  TLorentzVector ljet_1_p4;
  ljet_1_p4.SetPtEtaPhiE(JetPt->at(1),JetEta->at(1),JetPhi->at(1),JetE->at(1));
  ljet_0_p4.SetPtEtaPhiE(JetPt->at(0),JetEta->at(0),JetPhi->at(0),JetE->at(0));
  ljet_0_p4 = toGeV(ljet_0_p4);
  ljet_1_p4 = toGeV(ljet_1_p4);
  // Tau vectors
  TLorentzVector tau_0_p4;
  TLorentzVector tau_1_p4;
  tau_0_p4.SetPtEtaPhiE(TauPt->at(0),TauEta->at(0),TauPhi->at(0),TauE->at(0));
  tau_1_p4.SetPtEtaPhiE(TauPt->at(1),TauEta->at(1),TauPhi->at(1),TauE->at(1));
  tau_0_p4 = toGeV(tau_0_p4);
  tau_1_p4 = toGeV(tau_1_p4);

  // MET vector
  TLorentzVector met_reco_p4;
  met_reco_p4.SetPtEtaPhiE(MET_met,0,MET_phi,MET_met);
  met_reco_p4 = toGeV(met_reco_p4);
  
  //Charges and lepton ID
  float qtau0=TauCharge->at(0);
  float qtau1=TauCharge->at(1);
  // Tau medium working point
  bool tau_0_passed_medium_RNN = PassTauID("medium", TauRNNJetScore->at(0), TauNCoreTracks->at(0));
  bool tau_1_passed_medium_RNN = PassTauID("medium", TauRNNJetScore->at(1), TauNCoreTracks->at(1));

  std::size_t nTaus = TauPt->size();

  if (qtau0!=qtau1 && nTaus==2 && tau_0_passed_medium_RNN && tau_1_passed_medium_RNN && JetPt->size()<4 && n_bjets==0){
    // Angle between taus
    double angle=del_phi(tau_0_p4.Phi(),tau_1_p4.Phi());

    //trigger decision
    bool trigger_decision= passTrigger;
    // INVARIANT MASS 2-JETS
    double mjj=sqrt(2*(ljet_0_p4.Dot(ljet_1_p4)));

    if (mjj>=250 && trigger_decision) {

      // ZpT calculations
      double truth_z_pt=0.0;

      double Z_pt = (tau_1_p4 + tau_0_p4).Pt();
      if (z_sample==0) truth_z_pt=Z_pt;
    
      // Neutrinos
      std::pair<TLorentzVector,TLorentzVector> neus_p4;
      neus_p4 = CalculateNeutrinoVector(met_reco_p4, tau_0_p4, tau_1_p4);
      TLorentzVector tau_0_reco_p4 = tau_0_p4 + neus_p4.first;
      TLorentzVector tau_1_reco_p4 = tau_1_p4 + neus_p4.second;

      // TAU-TAU INVARIANT MASS
      double inv_tautau=sqrt((2*tau_1_p4.Pt()*tau_0_p4.Pt())*(cosh(tau_1_p4.Eta()-tau_0_p4.Eta())-cos(tau_1_p4.Phi()-tau_0_p4.Phi())));

      // TAU-TAU invariant mass with neutrinos
      TLorentzVector ttvv_p4 = tau_0_p4 + tau_1_p4 + neus_p4.first + neus_p4.second;
      double inv_tt_reco = ttvv_p4.Mag();

      //TAU-Tau invariant mass with met
      TLorentzVector ttmet_p4 = tau_0_p4 + tau_1_p4 + met_reco_p4;
      double inv_ttmet = ttmet_p4.Mag();

      // Minimum DeltaR between lepton and jets
      double min_dR_tau = min_deltaR(tau_0_p4,ljet_0_p4,ljet_1_p4);
      double min_dR_lep = min_deltaR(tau_1_p4,ljet_0_p4,ljet_1_p4);

      // Transverse mass
      double transverseMassTau1 = sqrt(2*tau_1_p4.Pt()*met_reco_p4.Pt()*(1-cos(tau_1_p4.Phi()-met_reco_p4.Phi())));

      // Handling BDT
      //float bdt_transmasstau1 = inv_tautau > 200 ? transverseMassTau1/std::pow(inv_tautau,0.3) : transverseMassTau1/std::pow(200,0.3); // for transverse-reco mass ratio
      //m_vbfBDT.update(mjj, 0.0, 0.0, 0.0, 0.0, bdt_transmasstau1, eventNumber);
      //double VBFBDT_score = m_vbfBDT.evaluate();

      // Rapidity seperation jets
      double delta_yjj = abs(ljet_0_p4.Rapidity()-ljet_1_p4.Rapidity());

      //Jet angle and Pseudo-Rapidity difference
      double delta_phijj = del_phi(ljet_0_p4.Phi(),ljet_1_p4.Phi());
      double delta_etajj = abs(ljet_0_p4.Eta()-ljet_1_p4.Eta());

      // Rapidity seperation taus
      double delta_ytt = abs(tau_0_reco_p4.Rapidity()-tau_1_reco_p4.Rapidity());
      double tautau_rapidity = (tau_0_reco_p4+tau_1_reco_p4).Rapidity();

      // Z centrality
      double z_centrality = abs(tautau_rapidity-0.5*(ljet_0_p4.Rapidity()+ljet_1_p4.Rapidity())) / delta_yjj;

      // N gap jets
      bool N_gapjets = CalculateNGapJets(ljet_0_p4.Rapidity(), ljet_1_p4.Rapidity(), JetEta);

      bool region_cut = Region("SR", N_gapjets, z_centrality);

      // Calculating pT balance
      std::vector<TLorentzVector> particles{ljet_0_p4,ljet_1_p4,tau_0_reco_p4,tau_1_reco_p4};
      if (N_gapjets) 
      {
        TLorentzVector ljet_2_p4;
        ljet_2_p4.SetPtEtaPhiE(JetPt->at(2),JetEta->at(2),JetPhi->at(2),JetE->at(2));
        ljet_2_p4 = toGeV(ljet_2_p4);
        particles.push_back(ljet_2_p4);
      }
      double pt_balance = CalculatePtBalance(particles);

      // Omega
      double omega = CalculateOmega(tau_0_p4, tau_1_p4, met_reco_p4);

      // MET angle
      double MET_angle = std::min(del_phi(met_reco_p4.Phi(),tau_0_p4.Phi()),del_phi(met_reco_p4.Phi(),tau_1_p4.Phi()));
    
      // cuts

      bool diLeptonMassRequirement = inv_tt_reco > 135;
      // Cuts vector
      std::vector<int> cuts={0,0,0,0,0,0,0,0,0,0};
      // CUTS
      if(tau_0_p4.Pt()>=80){cuts[0]=1;} //80
      if(tau_1_p4.Pt()>=50){cuts[1]=1;} //50
      if(ljet_0_p4.Pt()>=75){cuts[2]=1;}
      if(ljet_1_p4.Pt()>=70){cuts[3]=1;}
      if (diLeptonMassRequirement){cuts[4]=1;}
      if (delta_yjj > 2){cuts[5]=1;}
      if (omega >= -0.4 && omega <= 1.4){cuts[6]=1;}
      if (pt_balance <= 0.15){cuts[7]=1;}
      if (region_cut){cuts[8]=1;}
      if(mjj>=1000){cuts[9]=1;} 

      // SUM OF THE VECTOR STORING IF CUTS PASS OR NOT
      size_t sum{0};
      for(auto &j : cuts){sum=sum+j;}

      std::vector<int> cutsVector{1};
      cutsVector.insert(cutsVector.end(),cuts.begin(),cuts.end());
      bool passedAllCuts = (sum+1==cutsVector.size());
      std::vector<int> notFullCutsVector{1,static_cast<int>(passedAllCuts)};

      // Test of neutrino energy splitting
      TLorentzVector total_p4 = neus_p4.first + neus_p4.second;
      TLorentzVector difference_p4 = met_reco_p4 - total_p4;
      
      // FILLING CONATINER HISTOGRAMS
    }
  }
}

void BDT_Sample::FillTree(double weight, int z_sample, const std::string& sampleName) {
  double pi=TMath::Pi();
  // Jet vectors
  TLorentzVector ljet_0_p4;
  TLorentzVector ljet_1_p4;
  ljet_1_p4.SetPtEtaPhiE(JetPt->at(1),JetEta->at(1),JetPhi->at(1),JetE->at(1));
  ljet_0_p4.SetPtEtaPhiE(JetPt->at(0),JetEta->at(0),JetPhi->at(0),JetE->at(0));
  ljet_0_p4 = toGeV(ljet_0_p4);
  ljet_1_p4 = toGeV(ljet_1_p4);
  // Tau vectors
  TLorentzVector tau_0_p4;
  TLorentzVector tau_1_p4;
  tau_0_p4.SetPtEtaPhiE(TauPt->at(0),TauEta->at(0),TauPhi->at(0),TauE->at(0));
  tau_1_p4.SetPtEtaPhiE(TauPt->at(1),TauEta->at(1),TauPhi->at(1),TauE->at(1));
  tau_0_p4 = toGeV(tau_0_p4);
  tau_1_p4 = toGeV(tau_1_p4);

  // MET vector
  TLorentzVector met_reco_p4;
  met_reco_p4.SetPtEtaPhiE(MET_met,0,MET_phi,MET_met);
  met_reco_p4 = toGeV(met_reco_p4);
  
  //Charges and lepton ID
  float qtau0=TauCharge->at(0);
  float qtau1=TauCharge->at(1);
  // Tau medium working point
  bool tau_0_passed_medium_RNN = PassTauID("medium", TauRNNJetScore->at(0), TauNCoreTracks->at(0));
  bool tau_1_passed_medium_RNN = PassTauID("medium", TauRNNJetScore->at(1), TauNCoreTracks->at(1));

  std::size_t nTaus = TauPt->size();

  if (qtau0!=qtau1 && nTaus==2 && tau_0_passed_medium_RNN && tau_1_passed_medium_RNN && JetPt->size()<4 && n_bjets==0){
    // Angle between taus
    double angle=del_phi(tau_0_p4.Phi(),tau_1_p4.Phi());

    //trigger decision
    bool trigger_decision= passTrigger;
    // INVARIANT MASS 2-JETS
    double mjj=sqrt(2*(ljet_0_p4.Dot(ljet_1_p4)));

    if (mjj>=250 && trigger_decision) {

      // ZpT calculations
      double truth_z_pt=0.0;

      double Z_pt = (tau_1_p4 + tau_0_p4).Pt();
      if (z_sample==0) truth_z_pt=Z_pt;
    
      // Neutrinos
      std::pair<TLorentzVector,TLorentzVector> neus_p4;
      neus_p4 = CalculateNeutrinoVector(met_reco_p4, tau_0_p4, tau_1_p4);
      TLorentzVector tau_0_reco_p4 = tau_0_p4 + neus_p4.first;
      TLorentzVector tau_1_reco_p4 = tau_1_p4 + neus_p4.second;

      // TAU-TAU INVARIANT MASS
      double inv_tautau=sqrt((2*tau_1_p4.Pt()*tau_0_p4.Pt())*(cosh(tau_1_p4.Eta()-tau_0_p4.Eta())-cos(tau_1_p4.Phi()-tau_0_p4.Phi())));

      // TAU-TAU invariant mass with neutrinos
      TLorentzVector ttvv_p4 = tau_0_p4 + tau_1_p4 + neus_p4.first + neus_p4.second;
      double inv_tt_reco = ttvv_p4.Mag();

      //TAU-Tau invariant mass with met
      TLorentzVector ttmet_p4 = tau_0_p4 + tau_1_p4 + met_reco_p4;
      double inv_ttmet = ttmet_p4.Mag();

      // Minimum DeltaR between lepton and jets
      double min_dR_tau = min_deltaR(tau_0_p4,ljet_0_p4,ljet_1_p4);
      double min_dR_lep = min_deltaR(tau_1_p4,ljet_0_p4,ljet_1_p4);

      // Transverse mass
      double transverseMassTau1 = sqrt(2*tau_1_p4.Pt()*met_reco_p4.Pt()*(1-cos(tau_1_p4.Phi()-met_reco_p4.Phi())));
      double transMassRatio = transverseMassTau1 / inv_tt_reco;
      double transMassRatioFunc = transverseMassTau1 / std::pow(inv_tt_reco,0.3);

      // Rapidity seperation jets
      double delta_yjj = abs(ljet_0_p4.Rapidity()-ljet_1_p4.Rapidity());

      //Jet angle and Pseudo-Rapidity difference
      double delta_phijj = del_phi(ljet_0_p4.Phi(),ljet_1_p4.Phi());
      double delta_etajj = abs(ljet_0_p4.Eta()-ljet_1_p4.Eta());

      // Rapidity seperation taus
      double delta_ytt = abs(tau_0_reco_p4.Rapidity()-tau_1_reco_p4.Rapidity());
      double tautau_rapidity = (tau_0_reco_p4+tau_1_reco_p4).Rapidity();

      // Z centrality
      double z_centrality = abs(tautau_rapidity-0.5*(ljet_0_p4.Rapidity()+ljet_1_p4.Rapidity())) / delta_yjj;

      // N gap jets
      bool N_gapjets = CalculateNGapJets(ljet_0_p4.Rapidity(), ljet_1_p4.Rapidity(), JetEta);

      bool region_cut = Region("SR", N_gapjets, z_centrality);

      // Calculating pT balance
      std::vector<TLorentzVector> particles{ljet_0_p4,ljet_1_p4,tau_0_reco_p4,tau_1_reco_p4};
      if (N_gapjets) 
      {
        TLorentzVector ljet_2_p4;
        ljet_2_p4.SetPtEtaPhiE(JetPt->at(2),JetEta->at(2),JetPhi->at(2),JetE->at(2));
        ljet_2_p4 = toGeV(ljet_2_p4);
        particles.push_back(ljet_2_p4);
      }
      double pt_balance = CalculatePtBalance(particles);

      // Omega
      double omega = CalculateOmega(tau_0_p4, tau_1_p4, met_reco_p4);

      // MET angle
      double MET_angle = std::min(del_phi(met_reco_p4.Phi(),tau_0_p4.Phi()),del_phi(met_reco_p4.Phi(),tau_1_p4.Phi()));
    
      // cuts

      bool diLeptonMassRequirement = inv_tt_reco > 116 && inv_tt_reco < 15000;
      // Cuts vector
      std::vector<int> cuts={0,0,0,0,0,0,0,0,0,0};
      // CUTS
      if(tau_0_p4.Pt()>=80){cuts[0]=1;} //80
      if(tau_1_p4.Pt()>=50){cuts[1]=1;} //50
      if(ljet_0_p4.Pt()>=75){cuts[2]=1;}
      if(ljet_1_p4.Pt()>=70){cuts[3]=1;}
      if (diLeptonMassRequirement){cuts[4]=1;}
      if (delta_yjj > 0){cuts[5]=1;} //delta_yjj > 2
      if (omega >= -0.4 && omega <= 1.4){cuts[6]=1;}
      if (pt_balance <= 0.15){cuts[7]=1;} //pt_balance <= 0.15
      if (true){cuts[8]=1;} //region_cut
      if(mjj>=500){cuts[9]=1;}

      // SUM OF THE VECTOR STORING IF CUTS PASS OR NOT
      size_t sum{0};
      for(auto &j : cuts){sum=sum+j;}

      std::vector<int> cutsVector{1};
      cutsVector.insert(cutsVector.end(),cuts.begin(),cuts.end());
      bool passedAllCuts = (sum+1==cutsVector.size());
      bool failedOnlyRegionCut = cuts[8]==0 && (sum == cutsVector.size() - 2);
      std::vector<int> notFullCutsVector{1,static_cast<int>(passedAllCuts)};

      if (passedAllCuts){
      // FILLING TTree
      // Check if sample is VBF Ztautau or VBF Higgs
      bool isVBF = sampleName.find("VBF") != std::string::npos;
      bool isZtautau = sampleName.find("Ztautau") != std::string::npos;
      bool isHiggs = sampleName.find("Higgs") != std::string::npos;
      if (isVBF) 
      {
        if (isZtautau || isHiggs) 
        {
          m_signalTree.m_mcWeight = weight;
          m_signalTree.m_mjj = mjj;
          m_signalTree.m_deltaPhiLT = angle;
          m_signalTree.m_t0RNNScore = TauRNNJetScore->at(0);
          m_signalTree.m_t1RNNScore = TauRNNJetScore->at(1);
          m_signalTree.m_transverseMassLep = transverseMassTau1;
          m_signalTree.m_massTauTau = inv_tt_reco;
          m_signalTree.m_tau0_pT = tau_0_p4.Pt();
          m_signalTree.m_tau1_pT = tau_1_p4.Pt();
          m_signalTree.m_jet0_pT = ljet_0_p4.Pt();
          m_signalTree.m_jet1_pT = ljet_1_p4.Pt();
          m_signalTree.m_met_pT = met_reco_p4.Pt();
          m_signalTree.m_event_number = eventNumber;
          m_signalTree.m_pt_balance = pt_balance;
          m_signalTree.m_jet_gap = delta_yjj;
          m_signalTree.m_omega = omega;
          m_signalTree.m_centrality = z_centrality;
          m_signalTree.m_total_jet_pT = (ljet_0_p4 + ljet_1_p4).Pt();
          m_signalTree.m_total_tau_pT = (tau_0_p4 + tau_1_p4).Pt();
          m_signalTree.m_total_tau_mom = ttvv_p4.Pt();
          m_signalTree.m_transMassRatio = transMassRatio;
          m_signalTree.m_transMassRatioFunc = transMassRatioFunc;
          // Fill tree
          m_signalTree.FillTree();
        } 
      } else {
        m_backgroundTree.m_mcWeight = weight;
        m_backgroundTree.m_mjj = mjj;
        m_backgroundTree.m_deltaPhiLT = angle;
        m_backgroundTree.m_t0RNNScore = TauRNNJetScore->at(0);
        m_backgroundTree.m_t1RNNScore = TauRNNJetScore->at(1);
        m_backgroundTree.m_transverseMassLep = transverseMassTau1;
        m_backgroundTree.m_massTauTau = inv_tt_reco;
        m_backgroundTree.m_tau0_pT = tau_0_p4.Pt();
        m_backgroundTree.m_tau1_pT = tau_1_p4.Pt();
        m_backgroundTree.m_jet0_pT = ljet_0_p4.Pt();
        m_backgroundTree.m_jet1_pT = ljet_1_p4.Pt();
        m_backgroundTree.m_met_pT = met_reco_p4.Pt();
        m_backgroundTree.m_event_number = eventNumber;
        m_backgroundTree.m_pt_balance = pt_balance;
        m_backgroundTree.m_jet_gap = delta_yjj;
        m_backgroundTree.m_omega = omega;
        m_backgroundTree.m_centrality = z_centrality;
        m_backgroundTree.m_total_jet_pT = (ljet_0_p4 + ljet_1_p4).Pt();
        m_backgroundTree.m_total_tau_pT = (tau_0_p4 + tau_1_p4).Pt();
        m_backgroundTree.m_total_tau_mom = ttvv_p4.Pt();
        m_backgroundTree.m_transMassRatio = transMassRatio;
        m_backgroundTree.m_transMassRatioFunc = transMassRatioFunc;
        // Fill tree
        m_backgroundTree.FillTree();
      }
    }
  }
}
}

void BDTCuts::Fill(double weight, int z_sample, const std::string& sampleName) {
  double pi=TMath::Pi();

  // Jet vectors
  TLorentzVector ljet_0_p4;
  TLorentzVector ljet_1_p4;
  ljet_1_p4.SetPtEtaPhiE(JetPt->at(1),JetEta->at(1),JetPhi->at(1),JetE->at(1));
  ljet_0_p4.SetPtEtaPhiE(JetPt->at(0),JetEta->at(0),JetPhi->at(0),JetE->at(0));
  ljet_0_p4 = toGeV(ljet_0_p4);
  ljet_1_p4 = toGeV(ljet_1_p4);
  // Tau vectors
  TLorentzVector tau_0_p4;
  TLorentzVector tau_1_p4;
  tau_0_p4.SetPtEtaPhiE(TauPt->at(0),TauEta->at(0),TauPhi->at(0),TauE->at(0));
  tau_1_p4.SetPtEtaPhiE(TauPt->at(1),TauEta->at(1),TauPhi->at(1),TauE->at(1));
  tau_0_p4 = toGeV(tau_0_p4);
  tau_1_p4 = toGeV(tau_1_p4);

  // MET vector
  TLorentzVector met_reco_p4;
  met_reco_p4.SetPtEtaPhiE(MET_met,0,MET_phi,MET_met);
  met_reco_p4 = toGeV(met_reco_p4);
  
  //Charges and lepton ID
  float qtau0=TauCharge->at(0);
  float qtau1=TauCharge->at(1);
  bool same_sign = qtau0 == qtau1;

  // Medium working point
  bool tau_0_passed_medium_ID = PassTauID("medium", TauRNNJetScore->at(0), TauNCoreTracks->at(0));
  bool tau_1_passed_medium_ID = PassTauID("medium", TauRNNJetScore->at(1), TauNCoreTracks->at(1));

  // Tight working point
  bool tau_0_passed_tight_ID = PassTauID("tight", TauRNNJetScore->at(0), TauNCoreTracks->at(0));
  bool tau_1_passed_tight_ID = PassTauID("tight", TauRNNJetScore->at(1), TauNCoreTracks->at(1));

  std::size_t nTaus = TauPt->size();

  if (nTaus==2 && JetPt->size()<4 && n_bjets==0 && tau_0_passed_medium_ID && tau_1_passed_medium_ID && !same_sign){
    // Angle between taus
    double angle=del_phi(tau_0_p4.Phi(),tau_1_p4.Phi());

    //trigger decision
    bool trigger_decision= passTrigger;
    // INVARIANT MASS 2-JETS
    double mjj=sqrt(2*(ljet_0_p4.Dot(ljet_1_p4)));

    if (mjj>=250 && trigger_decision) {
      
      // ZpT calculations
      double truth_z_pt=0.0;

      double Z_pt = (tau_1_p4 + tau_0_p4).Pt();
      if (z_sample==0) truth_z_pt=Z_pt;
    
      // Neutrinos
      std::pair<TLorentzVector,TLorentzVector> neus_p4;
      neus_p4 = CalculateNeutrinoVector(met_reco_p4, tau_0_p4, tau_1_p4);
      TLorentzVector tau_0_reco_p4 = tau_0_p4 + neus_p4.first;
      TLorentzVector tau_1_reco_p4 = tau_1_p4 + neus_p4.second;

      // TAU-TAU VISIBLE INVARIANT MASS
      double inv_tautau=sqrt((2*tau_1_p4.Pt()*tau_0_p4.Pt())*(cosh(tau_1_p4.Eta()-tau_0_p4.Eta())-cos(tau_1_p4.Phi()-tau_0_p4.Phi())));

      // TAU-TAU reconstructed invariant
      TLorentzVector ttvv_p4 = tau_0_p4 + tau_1_p4 + neus_p4.first + neus_p4.second;
      double inv_tt_reco = ttvv_p4.Mag();
      double total_tau_pt = ttvv_p4.Pt();

      double mass_ratio = inv_tt_reco / inv_tautau;

      // Minimum DeltaR between lepton and jets
      double min_dR_tau = min_deltaR(tau_0_p4,ljet_0_p4,ljet_1_p4);
      double min_dR_lep = min_deltaR(tau_1_p4,ljet_0_p4,ljet_1_p4);

      // Transverse mass
      double transverseMassTau1 = sqrt(2*tau_1_p4.Pt()*met_reco_p4.Pt()*(1-cos(tau_1_p4.Phi()-met_reco_p4.Phi())));

      // Rapidity seperation jets
      double delta_yjj = abs(ljet_0_p4.Rapidity()-ljet_1_p4.Rapidity());

      // Z centrality
      double tautau_rapidity = (tau_0_reco_p4+tau_1_reco_p4).Rapidity();
      double z_centrality = abs(tautau_rapidity-0.5*(ljet_0_p4.Rapidity()+ljet_1_p4.Rapidity())) / delta_yjj;

      // N gap jets
      bool N_gapjets = CalculateNGapJets(ljet_0_p4.Rapidity(), ljet_1_p4.Rapidity(), JetEta);

      bool SR_cut = Region("SR", N_gapjets, z_centrality);

      // Calculating pT balance
      std::vector<TLorentzVector> particles{ljet_0_p4,ljet_1_p4,tau_0_reco_p4,tau_1_reco_p4};
      if (N_gapjets) 
      {
        TLorentzVector ljet_2_p4;
        ljet_2_p4.SetPtEtaPhiE(JetPt->at(2),JetEta->at(2),JetPhi->at(2),JetE->at(2));
        ljet_2_p4 = toGeV(ljet_2_p4);
        particles.push_back(ljet_2_p4);
      }
      double pt_balance = CalculatePtBalance(particles);

      // Omega
      double omega = CalculateOmega(tau_0_p4, tau_1_p4, met_reco_p4);

      // MET angle
      double MET_angle = std::min(del_phi(met_reco_p4.Phi(),tau_0_p4.Phi()),del_phi(met_reco_p4.Phi(),tau_1_p4.Phi()));

      // Handling BDT
      m_vbfBDT.update(mjj, delta_yjj, pt_balance, z_centrality, eventNumber);
      double VBFBDT_score = m_vbfBDT.evaluate();

      // cuts

      bool diLeptonMassRequirement;
      //diLeptonMassRequirement = inv_tt_reco > 160;
      diLeptonMassRequirement = inv_tt_reco < 116 && inv_tt_reco > 66;

      // Cuts vector
      std::vector<int> cuts={0,0,0,0,0,0,0,0,0,0,0,0,0};
      // CUTS
      if(tau_0_p4.Pt()>=80){cuts[0]=1;} //80
      if(tau_1_p4.Pt()>=60){cuts[1]=1;} //60
      if(ljet_0_p4.Pt()>=75){cuts[2]=1;}
      if(ljet_1_p4.Pt()>=70){cuts[3]=1;}
      if(diLeptonMassRequirement){cuts[4]=1;}
      if(delta_yjj>2){cuts[5]=1;}
      if(omega >= -0.4 && omega <= 1.4){cuts[6]=1;}
      if(pt_balance <= 0.15){cuts[7]=1;}
      if(mass_ratio < 4.0){cuts[8]=1;}
      if(mjj>=750){cuts[9]=1;}
      if(SR_cut){cuts[10]=1;}
      if(true){cuts[11]=1;} //tau_0_passed_tight_ID && tau_1_passed_tight_ID
      if(VBFBDT_score > 0.3){cuts[12]=1;}

      // SUM OF THE VECTOR STORING IF CUTS PASS OR NOT
      size_t sum{0};
      for(auto &j : cuts){sum=sum+j;}

      std::vector<int> cutsVector{1};
      cutsVector.insert(cutsVector.end(),cuts.begin(),cuts.end());
      bool passedAllCuts = (sum+1==cutsVector.size());
      std::vector<int> notFullCutsVector{1,static_cast<int>(passedAllCuts)};

      // Test of neutrino energy splitting
      TLorentzVector total_p4 = neus_p4.first + neus_p4.second;
      TLorentzVector difference_p4 = met_reco_p4 - total_p4;
      
      // FILL RAW HISTOGRAMS
      if (passedAllCuts){
        nJets->Fill(JetPt->size(),weight);
        tau0Eta->Fill(tau_0_p4.Eta(),weight);
      }
      // FILLING CONATINER HISTOGRAMS
      int tau0NTracks = TauNCoreTracks->at(0);
      int tau1NTracks = TauNCoreTracks->at(1);

      tau0_ptContainer.Fill(tau_0_p4.Pt(),weight,cutsVector);
      tau1_ptContainer.Fill(tau_1_p4.Pt(),weight,cutsVector);
      ljet0_ptContainer.Fill(ljet_0_p4.Pt(),weight,cutsVector);
      ljet1_ptContainer.Fill(ljet_1_p4.Pt(),weight,cutsVector);
      reco_mass_ttContainer.Fill(inv_tt_reco,weight,cutsVector);
      jet_gapContainer.Fill(delta_yjj,weight,cutsVector);
      omegaContainer.Fill(omega,weight,cutsVector);
      pt_balanceContainer.Fill(pt_balance,weight,cutsVector);
      mass_fracContainer.Fill(mass_ratio,weight,cutsVector);
      mass_jjContainer.Fill(mjj,weight,cutsVector);
      if (N_gapjets == 0) z_centralityContainer.Fill(z_centrality,weight,cutsVector);
      if (z_centrality < 0.5) N_gapjetsContainer.Fill(N_gapjets,weight,cutsVector);

      if (tau0NTracks==1) tau0_rnn_1pContainer.Fill(TauRNNJetScore->at(0),weight,cutsVector);
      else tau0_rnn_3pContainer.Fill(TauRNNJetScore->at(0),weight,cutsVector);
      if (tau1NTracks==1) tau1_rnn_1pContainer.Fill(TauRNNJetScore->at(1),weight,cutsVector);
      else tau1_rnn_3pContainer.Fill(TauRNNJetScore->at(1),weight,cutsVector);

      BDTScoreContainer.Fill(VBFBDT_score,weight,cutsVector);

      n_bjetsContainer.Fill(n_bjets,weight,notFullCutsVector);
      delta_phiContainer.Fill(angle,weight,cutsVector);
      
      nEventsContainer.Fill(1,weight,notFullCutsVector);


      // Only for MC samples
      if (sampleName.substr(0,4)!="data"){
        if (tau0NTracks==1) tau_matched_1pContainer.Fill(TauRNNJetScore->at(0),weight,notFullCutsVector);
        if (tau1NTracks==3) tau_matched_3pContainer.Fill(TauRNNJetScore->at(1),weight,notFullCutsVector);
      }
    }
  }
}