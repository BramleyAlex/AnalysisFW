// Include the file that lets the program know about the data
#include <vector>
#include <algorithm>
#include "CLoop.h"
#include "OutputTree.h"

double del_phi(double phi_1, double phi_2);
double min_deltaR(TLorentzVector* test_particle, std::vector<UInt_t>& bool_vector_container,const std::vector<TLorentzVector*>& jet_container);
// Tutorial: https://github.com/diegobaronm/AnalysisFW/wiki/Mini%E2%80%90tutorial
float metProjectionClosestLepton(TLorentzVector* met, TLorentzVector* muon, TLorentzVector* tau);

void CLoop::Fill(double weight, int z_sample, const std::string& sampleName) {
  double pi=TMath::Pi();
  //Charges and lepton ID
  float ql=muon_0_q;
  float qtau=tau_0_q;
  bool lepton_id=muon_0_id_medium;

  if (ql!=qtau && n_muons==1 && n_taus_rnn_loose>=1 && lepton_id){
    // Angle between tau and muon
    double angle=del_phi(tau_0_p4->Phi(),muon_0_p4->Phi());

    //trigger decision
    bool trigger_decision= false;
    bool trigger_match= false;
    if (run_number>= 276262 && run_number<=284484) {
      trigger_decision= bool(HLT_mu20_iloose_L1MU15 | HLT_mu50);
      trigger_match=bool(muTrigMatch_0_HLT_mu20_iloose_L1MU15 | muTrigMatch_0_HLT_mu50);
    } else {
      trigger_decision= bool(HLT_mu26_ivarmedium | HLT_mu50);
      trigger_match=bool(muTrigMatch_0_HLT_mu26_ivarmedium | muTrigMatch_0_HLT_mu50);
    }
    // INVARIANT MASS 2-JETS
    double mjj=sqrt(2*(ljet_0_p4->Dot(*ljet_1_p4)));

    if (mjj>=250 && trigger_decision  && trigger_match) {

      // ZpT calculations
      double truth_z_pt=0.0;

      // truth ZpT definition
      if (z_sample==1 || z_sample==2) truth_z_pt=truth_Z_p4->Pt()/1000;
      double Z_pt = (*muon_0_p4 + *tau_0_p4).Pt();
      if (z_sample==0) truth_z_pt=Z_pt;
    

      // LEP-TAU INVARIANT MASS
      double inv_taulep=sqrt((2*muon_0_p4->Pt()*tau_0_p4->Pt())*(cosh(muon_0_p4->Eta()-tau_0_p4->Eta())-cos(muon_0_p4->Phi()-tau_0_p4->Phi())));

      // Minimum DeltaR between lepton and jets
      std::vector<UInt_t> is_jet_present{ljet_0,ljet_1,ljet_2};
      std::vector<TLorentzVector*> jet_container{ljet_0_p4,ljet_1_p4,ljet_2_p4};
      double min_dR_tau = min_deltaR(tau_0_p4,is_jet_present,jet_container);
      double min_dR_lep = min_deltaR(muon_0_p4,is_jet_present,jet_container);

      // Transverse mass
      double transverseMassLep = sqrt(2*muon_0_p4->Pt()*met_reco_p4->Pt()*(1-cos(muon_0_p4->Phi()-met_reco_p4->Phi())));

      // Tutorial: https://github.com/diegobaronm/AnalysisFW/wiki/Mini%E2%80%90tutorial
      float metProjecClosestLep = metProjectionClosestLepton(met_reco_p4,muon_0_p4,tau_0_p4);

      // Handling BDT
      float bdt_transmasslep = inv_taulep > 200 ? transverseMassLep/std::pow(inv_taulep,0.3) : transverseMassLep/std::pow(200,0.3); // for transverse-reco mass ratio
      m_vbfBDT.update(mjj, 0.0, 0.0, 0.0, 0.0, bdt_transmasslep, event_number);
      double VBFBDT_score = m_vbfBDT.evaluate();
    
      // Cuts vector
      // Tutorial: https://github.com/diegobaronm/AnalysisFW/wiki/Mini%E2%80%90tutorial
      std::vector<int> cuts={0,0,0,0,0,0,0};
      // CUTS
      if (angle<=2.5){cuts[0]=1;}
      if(muon_0_p4->Pt()>=30){cuts[1]=1;}
      if (tau_0_p4->Pt()>=35){cuts[2]=1;}
      if(ljet_0_p4->Pt()>=75){cuts[3]=1;}
      if(mjj>=1000){cuts[4]=1;} 
      bool diLeptonMassRequirement =  inv_taulep >= 70;
      if (diLeptonMassRequirement){cuts[5]=1;}
      // Tutorial: https://github.com/diegobaronm/AnalysisFW/wiki/Mini%E2%80%90tutorial
      if (metProjecClosestLep>=30){cuts[6]=1;}

      

      // SUM OF THE VECTOR STORING IF CUTS PASS OR NOT
      size_t sum{0};
      for(auto &j : cuts){sum=sum+j;}

      std::vector<int> cutsVector{1};
      cutsVector.insert(cutsVector.end(),cuts.begin(),cuts.end());
      bool passedAllCuts = (sum+1==cutsVector.size());
      std::vector<int> notFullCutsVector{1,static_cast<int>(passedAllCuts)};

      // FILL RAW HISTOGRAMS
      if (passedAllCuts){
        nJets->Fill(n_jets,weight);
        tauEta->Fill(tau_0_p4->Eta(),weight);
      }

      // FILLING CONATINER HISTOGRAMS
      // Tutorial: https://github.com/diegobaronm/AnalysisFW/wiki/Mini%E2%80%90tutorial
      metProjecClosestLepContainer.Fill(metProjecClosestLep,weight,cutsVector);
      lep_ptContainer.Fill(muon_0_p4->Pt(),weight,cutsVector);
      tau_ptContainer.Fill(tau_0_p4->Pt(),weight,cutsVector);
      n_bjetsContainer.Fill(n_bjets_MV2c10_FixedCutBEff_85,weight,notFullCutsVector);
      lepisoContainer.Fill(muon_0_iso_TightTrackOnly_FixedRad,weight,notFullCutsVector);
      delta_phiContainer.Fill(angle,weight,cutsVector);
      mass_jjContainer.Fill(mjj,weight,cutsVector);
      ljet0_ptContainer.Fill(ljet_0_p4->Pt(),weight,cutsVector);
      ljet1_ptContainer.Fill(ljet_1_p4->Pt(),weight,cutsVector);
      bdtContainer.Fill(VBFBDT_score,weight,notFullCutsVector);
      visibleMassContainer.Fill(inv_taulep,weight,cutsVector);
      lepTransMassContainer.Fill(transverseMassLep,weight,cutsVector);
      delta_R_lepjetContainer.Fill(min_dR_lep,weight,cutsVector);
      delta_R_taujetContainer.Fill(min_dR_tau,weight,cutsVector);
      Z_ptContainer.Fill(truth_z_pt,weight,notFullCutsVector);

      if (tau_0_n_charged_tracks==1) rnn_score_1pContainer.Fill(tau_0_jet_rnn_score_trans,weight,notFullCutsVector);
      if (tau_0_n_charged_tracks==3) rnn_score_3pContainer.Fill(tau_0_jet_rnn_score_trans,weight,notFullCutsVector);

      // Only for MC samples
      if (sampleName.substr(0,4)!="data"){
        if (tau_0_n_charged_tracks==1) tau_matched_1pContainer.Fill(tau_0_truth_isHadTau,weight,notFullCutsVector);
        if (tau_0_n_charged_tracks==3) tau_matched_3pContainer.Fill(tau_0_truth_isHadTau,weight,notFullCutsVector);
      }
    }
  }
}

void CLoop::Style(double lumFactor) {
  nJets->Write();
  tauEta->Write();
}

void CLoop::FillTree(double weight, int z_sample, const std::string& sampleName) {
  double pi=TMath::Pi();
  //Charges and lepton ID
  float ql=muon_0_q;
  float qtau=tau_0_q;
  bool lepton_id=muon_0_id_medium;
  size_t n_ljets=n_jets-n_bjets_MV2c10_FixedCutBEff_85;

  if (ql!=qtau && n_muons==1 && n_taus_rnn_loose>=1 && lepton_id){
    
    //angle
    double angle=del_phi(tau_0_p4->Phi(),muon_0_p4->Phi());

    //trigger decision
    bool trigger_decision= false;
    bool trigger_match= false;
    if (run_number>= 276262 && run_number<=284484) {
      trigger_decision= bool(HLT_mu20_iloose_L1MU15 | HLT_mu50);
      trigger_match=bool(muTrigMatch_0_HLT_mu20_iloose_L1MU15 | muTrigMatch_0_HLT_mu50);
    } else {
      trigger_decision= bool(HLT_mu26_ivarmedium | HLT_mu50);
      trigger_match=bool(muTrigMatch_0_HLT_mu26_ivarmedium | muTrigMatch_0_HLT_mu50);
    }
    // INVARIANT MASS 2-JETS
    double mjj=sqrt(2*(ljet_0_p4->Dot(*ljet_1_p4)));
    if (mjj>=250 && trigger_decision  && trigger_match) {
      // ZpT calculations
      double truth_z_pt=0.0;

       // truth ZpT definition
      if (z_sample==1 || z_sample==2) truth_z_pt=truth_Z_p4->Pt()/1000;
      double Z_pt = (*muon_0_p4 + *tau_0_p4).Pt();
      if (z_sample==0) truth_z_pt=Z_pt;
    

      // LEP-TAU INVARIANT MASS
      double inv_taulep=sqrt((2*muon_0_p4->Pt()*tau_0_p4->Pt())*(cosh(muon_0_p4->Eta()-tau_0_p4->Eta())-cos(muon_0_p4->Phi()-tau_0_p4->Phi())));

      // Minimum DeltaR between lepton and jets
      std::vector<UInt_t> is_jet_present{ljet_0,ljet_1,ljet_2};
      std::vector<TLorentzVector*> jet_container{ljet_0_p4,ljet_1_p4,ljet_2_p4};
      double min_dR_tau = min_deltaR(tau_0_p4,is_jet_present,jet_container);
      double min_dR_lep = min_deltaR(muon_0_p4,is_jet_present,jet_container);

      // Transverse mass
      double transverseMassLep = sqrt(2*muon_0_p4->Pt()*met_reco_p4->Pt()*(1-cos(muon_0_p4->Phi()-met_reco_p4->Phi())));

      // Tutorial: https://github.com/diegobaronm/AnalysisFW/wiki/Mini%E2%80%90tutorial
      float metProjecClosestLep = metProjectionClosestLepton(met_reco_p4,muon_0_p4,tau_0_p4);

        // Cuts vector
      std::vector<int> cuts={0,0,0,0,0,0};
      // CUTS
      if (angle<=2.5){cuts[0]=1;}
      if(muon_0_p4->Pt()>=30){cuts[1]=1;}
      if (tau_0_p4->Pt()>=35){cuts[2]=1;}
      if(ljet_0_p4->Pt()>=75){cuts[3]=1;}
      if(mjj>=1000){cuts[4]=1;} 
      bool diLeptonMassRequirement =  inv_taulep >= 70;
      if (diLeptonMassRequirement){cuts[5]=1;}

      // SUM OF THE VECTOR STORING IF CUTS PASS OR NOT
      size_t sum{0};
      for(auto &j : cuts){sum=sum+j;}

      std::vector<int> cutsVector{1};
      cutsVector.insert(cutsVector.end(),cuts.begin(),cuts.end());
      bool passedAllCuts = (sum+1==cutsVector.size());

      if (passedAllCuts){
      // FILLING TTree
      // Check if sample is VBF Ztautau
      bool isZtautau = sampleName.find("Ztautau") != std::string::npos;
      if (isZtautau) {
        m_signalTree.m_mcWeight = weight;
        m_signalTree.m_mjj = mjj;
        m_signalTree.m_deltaPhiLT = angle;
        m_signalTree.m_jetRNNScore = tau_0_jet_rnn_score_trans;
        m_signalTree.m_transverseMassLep = transverseMassLep;
        m_signalTree.m_massTauLep = inv_taulep;
        m_signalTree.m_tau_pT = tau_0_p4->Pt();
        m_signalTree.m_lep_pT = muon_0_p4->Pt();
        m_signalTree.m_jet0_pT = ljet_0_p4->Pt();
        m_signalTree.m_jet1_pT = ljet_1_p4->Pt();
        m_signalTree.m_met_pT = met_reco_p4->Pt();
        m_signalTree.m_event_number = event_number;
        // Tutorial: https://github.com/diegobaronm/AnalysisFW/wiki/Mini%E2%80%90tutorial
        m_signalTree.m_metProjection = metProjecClosestLep;
        // Fill tree
        m_signalTree.FillTree();
      } else {
        m_backgroundTree.m_mcWeight = weight;
        m_backgroundTree.m_mjj = mjj;
        m_backgroundTree.m_deltaPhiLT = angle;
        m_backgroundTree.m_jetRNNScore = tau_0_jet_rnn_score_trans;
        m_backgroundTree.m_transverseMassLep = transverseMassLep;
        m_backgroundTree.m_massTauLep = inv_taulep;
        m_backgroundTree.m_tau_pT = tau_0_p4->Pt();
        m_backgroundTree.m_lep_pT = muon_0_p4->Pt();
        m_backgroundTree.m_jet0_pT = ljet_0_p4->Pt();
        m_backgroundTree.m_jet1_pT = ljet_1_p4->Pt();
        m_backgroundTree.m_met_pT = met_reco_p4->Pt();
        m_backgroundTree.m_event_number = event_number;
        // Tutorial: https://github.com/diegobaronm/AnalysisFW/wiki/Mini%E2%80%90tutorial
        m_backgroundTree.m_metProjection = metProjecClosestLep;
        // Fill tree
        m_backgroundTree.FillTree();
      }
    }
  }
}
}