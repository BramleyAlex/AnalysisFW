// Include the file that lets the program know about the data
#include <vector>
#include <algorithm>
#include "CLoop.h"
#include "OutputTree.h"

double del_phi(double phi_1, double phi_2);
double min_deltaR(const TLorentzVector& test_particle, const TLorentzVector& jet1, const TLorentzVector& jet2);
TLorentzVector& toGeV(TLorentzVector &v);
double CalculatePtBalance(const std::vector<TLorentzVector> particles);
double CalculateOmega(const TLorentzVector& tau_0_p4, const TLorentzVector& tau_1_p4, const TLorentzVector& met_p4);
std::pair<TLorentzVector,TLorentzVector> CalculateNeutrinoVector(const TLorentzVector& met_p4, const TLorentzVector& tau_0_p4, 
    const TLorentzVector& tau_1_p4);
void PrintTLorentzVectorPtEtaPhi(const TLorentzVector& input);
void PrintTLorentzVectorXYZT(const TLorentzVector& input);

void CLoop::Fill(double weight, int z_sample, const std::string& sampleName) {
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

  // Tau loose working point
  bool passed_medium_tau_RNN{false};
  if (TauRNNJetScore->at(0) > 0.25 && TauNCoreTracks->at(0) == 1) passed_medium_tau_RNN = true;
  else if (TauRNNJetScore->at(0) > 0.40 && TauNCoreTracks->at(0) == 3) passed_medium_tau_RNN = true;

  std::size_t nTaus = TauPt->size();

  if (qtau0!=qtau1 && nTaus==2 && passed_medium_tau_RNN && JetPt->size()<4 && n_bjets==0){
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

      if (inv_tt_reco <= 140) //Blinding  for high mass region
      {
        //TAU-Tau invariant mass with met (should be the same as met_ttvv)
        TLorentzVector ttmet_p4 = tau_0_p4 + tau_1_p4 + met_reco_p4;
        double inv_ttmet = ttmet_p4.Mag();

        // Minimum DeltaR between lepton and jets
        double min_dR_tau = min_deltaR(tau_0_p4,ljet_0_p4,ljet_1_p4);
        double min_dR_lep = min_deltaR(tau_1_p4,ljet_0_p4,ljet_1_p4);

        // Transverse mass
        double transverseMassTau1 = sqrt(2*tau_1_p4.Pt()*met_reco_p4.Pt()*(1-cos(tau_1_p4.Phi()-met_reco_p4.Phi())));

        // Handling BDT
        float bdt_transmasstau1 = inv_tautau > 200 ? transverseMassTau1/std::pow(inv_tautau,0.3) : transverseMassTau1/std::pow(200,0.3); // for transverse-reco mass ratio
        m_vbfBDT.update(mjj, 0.0, 0.0, 0.0, 0.0, bdt_transmasstau1, eventNumber);
        double VBFBDT_score = m_vbfBDT.evaluate();

        // Rapidity seperation jets
        double delta_yjj = abs(ljet_0_p4.Rapidity()-ljet_1_p4.Rapidity());

        // Rapidity seperation taus and centrality
        double delta_ytt = abs(tau_0_reco_p4.Rapidity()-tau_1_reco_p4.Rapidity());
        double tautau_rapidity = (tau_0_reco_p4+tau_1_reco_p4).Rapidity();

        double z_centrality = abs(tautau_rapidity-0.5*(ljet_0_p4.Rapidity()+ljet_1_p4.Rapidity())) / delta_yjj;

        // N gap jets
        std::vector<TLorentzVector> particles{ljet_0_p4,ljet_1_p4,tau_0_reco_p4,tau_1_reco_p4};

        bool N_gapjets{false};
        if (JetPt->size()>=3)
        {
          double lower_bound;
          double upper_bound;
          if (ljet_0_p4.Rapidity() < ljet_1_p4.Rapidity())
          {
            lower_bound = ljet_0_p4.Rapidity();
            upper_bound = ljet_1_p4.Rapidity();
          }
          else
          {
            lower_bound = ljet_1_p4.Rapidity();
            upper_bound = ljet_0_p4.Rapidity();
          }
          TLorentzVector ljet_2_p4;
          ljet_2_p4.SetPtEtaPhiE(JetPt->at(2),JetEta->at(2),JetPhi->at(2),JetE->at(2));
          ljet_2_p4 = toGeV(ljet_2_p4);
          if (ljet_2_p4.Rapidity() >= lower_bound && ljet_2_p4.Rapidity() <= upper_bound) 
          {
            N_gapjets = true;
            particles.push_back(ljet_2_p4);
          }
        }

        // Calculating pT balance
        double pt_balance = CalculatePtBalance(particles);

        // Omega
        double omega = CalculateOmega(tau_0_p4, tau_1_p4, met_reco_p4);

        // MET angle
        double MET_angle = std::min(del_phi(met_reco_p4.Phi(),tau_0_p4.Phi()),del_phi(met_reco_p4.Phi(),tau_1_p4.Phi()));
      
        // Cuts vector
        std::vector<int> cuts={0,0,0,0,0,0,0,0,0};
        // CUTS
        if(tau_0_reco_p4.Pt()>=80){cuts[0]=1;}
        if(tau_1_reco_p4.Pt()>=40){cuts[1]=1;}
        if(ljet_0_p4.Pt()>=75){cuts[2]=1;}
        if(ljet_1_p4.Pt()>=70){cuts[3]=1;}
        if(mjj>=1000){cuts[4]=1;} 
        bool diLeptonMassRequirement =  inv_tt_reco >= 66 && inv_tt_reco <=116;
        if (diLeptonMassRequirement){cuts[5]=1;}
        if (delta_yjj > 2){cuts[6]=1;}
        if (omega >= -0.4 && omega <= 1.4){cuts[7]=1;}
        if (pt_balance <= 0.15){cuts[8]=1;}

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
        tau0_ptContainer.Fill(tau_0_p4.Pt(),weight,cutsVector);
        tau1_ptContainer.Fill(tau_1_p4.Pt(),weight,cutsVector);
        n_bjetsContainer.Fill(n_bjets,weight,notFullCutsVector);
        delta_phiContainer.Fill(angle,weight,cutsVector);
        mass_jjContainer.Fill(mjj,weight,cutsVector);
        ljet0_ptContainer.Fill(ljet_0_p4.Pt(),weight,cutsVector);
        ljet1_ptContainer.Fill(ljet_1_p4.Pt(),weight,cutsVector);
        bdtContainer.Fill(VBFBDT_score,weight,notFullCutsVector);
        visibleMassContainer.Fill(inv_tautau,weight,cutsVector);
        tau1TransMassContainer.Fill(transverseMassTau1,weight,cutsVector);
        delta_R_lepjetContainer.Fill(min_dR_lep,weight,cutsVector);
        delta_R_taujetContainer.Fill(min_dR_tau,weight,cutsVector);
        Z_ptContainer.Fill(truth_z_pt,weight,notFullCutsVector);
        pt_balanceContainer.Fill(pt_balance,weight,cutsVector);
        jet_gapContainer.Fill(delta_yjj,weight,cutsVector);
        omegaContainer.Fill(omega,weight,cutsVector);
        if (omega >= 0 && omega <=1) neu0_pt_insideContainer.Fill(neus_p4.first.Pt(),weight,notFullCutsVector);
        else neu0_pt_outsideContainer.Fill(neus_p4.first.Pt(),weight,notFullCutsVector);
        if (omega >= 0 && omega <=1) neu1_pt_insideContainer.Fill(neus_p4.second.Pt(),weight,notFullCutsVector);
        else neu1_pt_outsideContainer.Fill(neus_p4.second.Pt(),weight,notFullCutsVector);

        if (omega >= 0 && omega <=1) tau0_reco_pt_inContainer.Fill(tau_0_reco_p4.Pt(),weight,notFullCutsVector);
        else tau0_reco_pt_outContainer.Fill(tau_0_reco_p4.Pt(),weight,notFullCutsVector);
        if (omega >= 0 && omega <=1) tau1_reco_pt_inContainer.Fill(tau_1_reco_p4.Pt(),weight,notFullCutsVector);
        else tau1_reco_pt_outContainer.Fill(tau_1_reco_p4.Pt(),weight,notFullCutsVector);

        if (omega >= 0 && omega <=1) mass_tt_reco_insideContainer.Fill(inv_tt_reco,weight,cutsVector);
        else mass_tt_reco_outsideContainer.Fill(inv_tt_reco,weight,cutsVector);
        met_angleContainer.Fill(MET_angle,weight,notFullCutsVector);
        z_centralityContainer.Fill(z_centrality,weight,notFullCutsVector);
        N_gapjetsContainer.Fill(N_gapjets,weight,notFullCutsVector);

        int tau0NTracks = TauNCoreTracks->at(0);
        int tau1NTracks = TauNCoreTracks->at(1);

        if (tau0NTracks==1) rnn_score_1pContainer.Fill(TauRNNJetScore->at(0),weight,notFullCutsVector);
        if (tau1NTracks==3) rnn_score_3pContainer.Fill(TauRNNJetScore->at(1),weight,notFullCutsVector);

        // Only for MC samples
        if (sampleName.substr(0,4)!="data"){
          if (tau0NTracks==1) tau_matched_1pContainer.Fill(TauRNNJetScore->at(0),weight,notFullCutsVector);
          if (tau1NTracks==3) tau_matched_3pContainer.Fill(TauRNNJetScore->at(1),weight,notFullCutsVector);
        }
      }
    }
  }
}

void CLoop::Style(double lumFactor) {
  nJets->Write();
  tau0Eta->Write();
}

void CLoop::FillTree(double weight, int z_sample, const std::string& sampleName) {
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

  //Tau loose working point
  bool passed_loose_tau_RNN{false};
  if (TauRNNJetScore->at(0) > 0.15 && TauNCoreTracks->at(0) == 1) passed_loose_tau_RNN = true;
  else if (TauRNNJetScore->at(0) > 0.25 && TauNCoreTracks->at(0) == 3) passed_loose_tau_RNN = true;

  std::size_t nTaus = TauPt->size();

  if (qtau0!=qtau1 && nTaus==2 && passed_loose_tau_RNN && n_bjets==0){
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
    

      // LEP-TAU INVARIANT MASS
      double inv_tautau=sqrt((2*tau_1_p4.Pt()*tau_0_p4.Pt())*(cosh(tau_1_p4.Eta()-tau_0_p4.Eta())-cos(tau_1_p4.Phi()-tau_0_p4.Phi())));

      // Minimum DeltaR between lepton and jets
      double min_dR_tau = min_deltaR(tau_0_p4,ljet_0_p4,ljet_1_p4);
      double min_dR_lep = min_deltaR(tau_1_p4,ljet_0_p4,ljet_1_p4);

      // Transverse mass
      double transverseMassTau1 = sqrt(2*tau_1_p4.Pt()*met_reco_p4.Pt()*(1-cos(tau_1_p4.Phi()-met_reco_p4.Phi())));

      // Pt balance
      std::vector<TLorentzVector> particles{ljet_0_p4,ljet_1_p4,tau_0_p4,tau_1_p4};
      double pt_balance = CalculatePtBalance(particles);

      // Rapidity seperation
      double delta_y = abs(ljet_0_p4.Rapidity()-ljet_1_p4.Rapidity());
    
      // Cuts vector
      std::vector<int> cuts={0,0,0,0,0,0,0,0,0};
      // CUTS
      if(tau_1_p4.Pt()>=30){cuts[0]=1;}
      if (tau_0_p4.Pt()>=35){cuts[1]=1;}
      if(ljet_0_p4.Pt()>=75){cuts[2]=1;}
      if(ljet_1_p4.Pt()>=70){cuts[3]=1;}
      if(mjj>=1000){cuts[4]=1;} 
      bool diLeptonMassRequirement =  inv_tautau >= 70;
      if (diLeptonMassRequirement){cuts[5]=1;}
      if (delta_y > 2){cuts[6]=1;}
      
      // SUM OF THE VECTOR STORING IF CUTS PASS OR NOT
      size_t sum{0};
      for(auto &j : cuts){sum=sum+j;}

      std::vector<int> cutsVector{1};
      cutsVector.insert(cutsVector.end(),cuts.begin(),cuts.end());
      bool passedAllCuts = (sum+1==cutsVector.size());
      std::vector<int> notFullCutsVector{1,static_cast<int>(passedAllCuts)};

      if (passedAllCuts){
      // FILLING TTree
      // Check if sample is VBF Ztautau
      bool isVBF = sampleName.find("VBF") != std::string::npos;
      bool isZtautau = sampleName.find("Ztautau") != std::string::npos;
      if (isVBF) 
      {
        if (isZtautau) 
        {
          m_signalTree.m_mcWeight = weight;
          m_signalTree.m_mjj = mjj;
          m_signalTree.m_deltaPhiLT = angle;
          m_signalTree.m_jetRNNScore = TauRNNJetScore->at(0);
          m_signalTree.m_transverseMassLep = transverseMassTau1;
          m_signalTree.m_massTauLep = inv_tautau;
          m_signalTree.m_tau_pT = tau_0_p4.Pt();
          m_signalTree.m_lep_pT = tau_1_p4.Pt();
          m_signalTree.m_jet0_pT = ljet_0_p4.Pt();
          m_signalTree.m_jet1_pT = ljet_1_p4.Pt();
          m_signalTree.m_met_pT = met_reco_p4.Pt();
          m_signalTree.m_event_number = eventNumber;
          m_signalTree.m_pt_balance = pt_balance;
          m_signalTree.m_jet_gap = delta_y;
          // Fill tree
          m_signalTree.FillTree();
        } 
      } else {
        m_backgroundTree.m_mcWeight = weight;
        m_backgroundTree.m_mjj = mjj;
        m_backgroundTree.m_deltaPhiLT = angle;
        m_backgroundTree.m_jetRNNScore = TauRNNJetScore->at(0);
        m_backgroundTree.m_transverseMassLep = transverseMassTau1;
        m_backgroundTree.m_massTauLep = inv_tautau;
        m_backgroundTree.m_tau_pT = tau_0_p4.Pt();
        m_backgroundTree.m_lep_pT = tau_1_p4.Pt();
        m_backgroundTree.m_jet0_pT = ljet_0_p4.Pt();
        m_backgroundTree.m_jet1_pT = ljet_1_p4.Pt();
        m_backgroundTree.m_met_pT = met_reco_p4.Pt();
        m_backgroundTree.m_event_number = eventNumber;
        m_backgroundTree.m_pt_balance = pt_balance;
        m_backgroundTree.m_jet_gap = delta_y;
        // Fill tree
        m_backgroundTree.FillTree();
      }
    }
  }
}
}