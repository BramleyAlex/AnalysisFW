#include <memory>
#include "CLoop.h"

CLoop::CLoop(TTree *tree,std::string sample_name) : fChain(0)
{
   Init(tree,sample_name);
}

CLoop::~CLoop()
{
   if (!fChain) return;
}

Int_t CLoop::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}

Long64_t CLoop::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void CLoop::Init(TTree *tree,std::string sample_name)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   TauEta = 0;
   TauPhi = 0;
   TauPt = 0;
   TauE = 0;
   TauMinMuonDR = 0;
   TauMinMuonDPhi = 0;
   TauMinMuonDEta = 0;
   TauMinMuonPhi = 0;
   TauMinMuonEta = 0;
   TauMinMuonPt = 0;
   TauMinMuonE = 0;
   TauMinMuonAuthor = 0;
   TauMinMuonQuality = 0;
   TauMinMuonType = 0;
   TauNCoreTracks = 0;
   TauCharge = 0;
   TauDecayMode = 0;
   TauBDTEleScore = 0;
   TauBDTJetScore = 0;
   TauRNNJetScore = 0;
   TauWP = 0;
   TauRecoSF = 0;
   TauTrackWidthPt1000PV = 0;
   TauTrackWidthPt500PV = 0;
   TauTrackWidthPt1000TV = 0;
   TauTrackWidthPt500TV = 0;
   Tau_trigSF_HLT_tau80L1TAU60_medium1_tracktwo = 0;
   Tau_trigSF_HLT_tau125_medium1_tracktwo = 0;
   Tau_trigSF_HLT_tau160_medium1_tracktwo = 0;
   TauMatchedTriggers = 0;
   mmc_mass_mlm = 0;
   mmc_mass_maxw = 0;
   mmc_taupt_maxw = 0;
   mmc_tauphi_maxw = 0;
   mmc_taueta_maxw = 0;
   mmc_taumass_maxw = 0;
   mmc_mass_mlnu3m = 0;
   mmc_taupt_mlnu3m = 0;
   mmc_tauphi_mlnu3m = 0;
   mmc_taueta_mlnu3m = 0;
   mmc_taumass_mlnu3m = 0;
   JetEta = 0;
   JetPhi = 0;
   JetPt = 0;
   JetE = 0;
   Jet_btag = 0;
   Jet_btSF = 0;
   Jet_JVT = 0;
   MuonEta = 0;
   MuonPhi = 0;
   MuonPt = 0;
   MuonE = 0;
   MuonCharge = 0;
   Muon_d0sig = 0;
   Muon_delta_z0 = 0;
   Muon_recoSF = 0;
   Muon_isoSF = 0;
   MuonMatchedTriggers = 0;
   Muon_ttvaSF = 0;
   EleEta = 0;
   ElePhi = 0;
   ElePt = 0;
   EleE = 0;
   EleCharge = 0;
   Ele_d0sig = 0;
   Ele_delta_z0 = 0;
   Ele_recoSF = 0;
   Ele_idSF = 0;
   Ele_isoSF = 0;
   Ele_trigSF = 0;
   EleMatchedTriggers = 0;
   PhotonEta = 0;
   PhotonPhi = 0;
   PhotonPt = 0;
   PhotonE = 0;
   TruthJetEta = 0;
   TruthJetPhi = 0;
   TruthJetPt = 0;
   TruthJetE = 0;
   TruthNeutrinoEta = 0;
   TruthNeutrinoPhi = 0;
   TruthNeutrinoPt = 0;
   TruthNeutrinoE = 0;
   TruthNeutrinoCharge = 0;
   TruthMuonEta = 0;
   TruthMuonPhi = 0;
   TruthMuonPt = 0;
   TruthMuonE = 0;
   TruthMuonCharge = 0;
   TruthEleEta = 0;
   TruthElePhi = 0;
   TruthElePt = 0;
   TruthEleE = 0;
   TruthEleCharge = 0;
   TruthTauEta = 0;
   TruthTauPhi = 0;
   TruthTauPt = 0;
   TruthTauM = 0;
   VisTruthTauEta = 0;
   VisTruthTauPhi = 0;
   VisTruthTauPt = 0;
   VisTruthTauM = 0;
   TruthTauCharge = 0;
   TruthTau_isHadronic = 0;
   TruthTau_decay_mode = 0;
   DNN_outputs = 0;
   MDN_outputs = 0;
   DS_outputs = 0;
   MDN_covariance = 0;
   PassedTriggers = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("mcWeight", &mcWeight, &b_mcWeight);
   fChain->SetBranchAddress("rwCorr", &rwCorr, &b_rwCorr);
   fChain->SetBranchAddress("prwWeight", &prwWeight, &b_prwWeight);
   fChain->SetBranchAddress("jvtSF", &jvtSF, &b_jvtSF);
   fChain->SetBranchAddress("fjvtSF", &fjvtSF, &b_fjvtSF);
   fChain->SetBranchAddress("selectionSF", &selectionSF, &b_selectionSF);
   fChain->SetBranchAddress("weight", &weight, &b_weight);
   fChain->SetBranchAddress("runNumber", &runNumber, &b_runNumber);
   fChain->SetBranchAddress("mcChannel", &mcChannel, &b_mcChannel);
   fChain->SetBranchAddress("eventNumber", &eventNumber, &b_eventNumber);
   fChain->SetBranchAddress("nVtx", &nVtx, &b_nVtx);
   fChain->SetBranchAddress("stitch_mll_pass", &stitch_mll_pass, &b_stitch_mll_pass);
   fChain->SetBranchAddress("mu", &mu, &b_mu);
   fChain->SetBranchAddress("amu", &amu, &b_amu);
   fChain->SetBranchAddress("amuScaled", &amuScaled, &b_amuScaled);
   fChain->SetBranchAddress("passTruth", &passTruth, &b_passTruth);
   fChain->SetBranchAddress("passReco", &passReco, &b_passReco);
   fChain->SetBranchAddress("passTrigger", &passTrigger, &b_passTrigger);
   fChain->SetBranchAddress("TauEta", &TauEta, &b_TauEta);
   fChain->SetBranchAddress("TauPhi", &TauPhi, &b_TauPhi);
   fChain->SetBranchAddress("TauPt", &TauPt, &b_TauPt);
   fChain->SetBranchAddress("TauE", &TauE, &b_TauE);
   fChain->SetBranchAddress("TauMinMuonDR", &TauMinMuonDR, &b_TauMinMuonDR);
   fChain->SetBranchAddress("TauMinMuonDPhi", &TauMinMuonDPhi, &b_TauMinMuonDPhi);
   fChain->SetBranchAddress("TauMinMuonDEta", &TauMinMuonDEta, &b_TauMinMuonDEta);
   fChain->SetBranchAddress("TauMinMuonPhi", &TauMinMuonPhi, &b_TauMinMuonPhi);
   fChain->SetBranchAddress("TauMinMuonEta", &TauMinMuonEta, &b_TauMinMuonEta);
   fChain->SetBranchAddress("TauMinMuonPt", &TauMinMuonPt, &b_TauMinMuonPt);
   fChain->SetBranchAddress("TauMinMuonE", &TauMinMuonE, &b_TauMinMuonE);
   fChain->SetBranchAddress("TauMinMuonAuthor", &TauMinMuonAuthor, &b_TauMinMuonAuthor);
   fChain->SetBranchAddress("TauMinMuonQuality", &TauMinMuonQuality, &b_TauMinMuonQuality);
   fChain->SetBranchAddress("TauMinMuonType", &TauMinMuonType, &b_TauMinMuonType);
   fChain->SetBranchAddress("TauNCoreTracks", &TauNCoreTracks, &b_TauNCoreTracks);
   fChain->SetBranchAddress("TauCharge", &TauCharge, &b_TauCharge);
   fChain->SetBranchAddress("TauDecayMode", &TauDecayMode, &b_TauDecayMode);
   fChain->SetBranchAddress("TauBDTEleScore", &TauBDTEleScore, &b_TauBDTEleScore);
   fChain->SetBranchAddress("TauBDTJetScore", &TauBDTJetScore, &b_TauBDTJetScore);
   fChain->SetBranchAddress("TauRNNJetScore", &TauRNNJetScore, &b_TauRNNJetScore);
   fChain->SetBranchAddress("TauWP", &TauWP, &b_TauWP);
   fChain->SetBranchAddress("TauRecoSF", &TauRecoSF, &b_TauRecoSF);
   fChain->SetBranchAddress("TauTrackWidthPt1000PV", &TauTrackWidthPt1000PV, &b_TauTrackWidthPt1000PV);
   fChain->SetBranchAddress("TauTrackWidthPt500PV", &TauTrackWidthPt500PV, &b_TauTrackWidthPt500PV);
   fChain->SetBranchAddress("TauTrackWidthPt1000TV", &TauTrackWidthPt1000TV, &b_TauTrackWidthPt1000TV);
   fChain->SetBranchAddress("TauTrackWidthPt500TV", &TauTrackWidthPt500TV, &b_TauTrackWidthPt500TV);
   fChain->SetBranchAddress("Tau_trigSF_HLT_tau80L1TAU60_medium1_tracktwo", &Tau_trigSF_HLT_tau80L1TAU60_medium1_tracktwo, &b_Tau_trigSF_HLT_tau80L1TAU60_medium1_tracktwo);
   fChain->SetBranchAddress("Tau_trigSF_HLT_tau125_medium1_tracktwo", &Tau_trigSF_HLT_tau125_medium1_tracktwo, &b_Tau_trigSF_HLT_tau125_medium1_tracktwo);
   fChain->SetBranchAddress("Tau_trigSF_HLT_tau160_medium1_tracktwo", &Tau_trigSF_HLT_tau160_medium1_tracktwo, &b_Tau_trigSF_HLT_tau160_medium1_tracktwo);
   fChain->SetBranchAddress("TauMatchedTriggers", &TauMatchedTriggers, &b_TauMatchedTriggers);
   fChain->SetBranchAddress("mll", &mll, &b_mll);
   fChain->SetBranchAddress("lbjet_pt", &lbjet_pt, &b_lbjet_pt);
   fChain->SetBranchAddress("lnbjet_pt", &lnbjet_pt, &b_lnbjet_pt);
   fChain->SetBranchAddress("dphi_Muon_MET", &dphi_Muon_MET, &b_dphi_Muon_MET);
   fChain->SetBranchAddress("mmc_mass_mlm", &mmc_mass_mlm, &b_mmc_mass_mlm);
   fChain->SetBranchAddress("mmc_mass_maxw", &mmc_mass_maxw, &b_mmc_mass_maxw);
   fChain->SetBranchAddress("mmc_taupt_maxw", &mmc_taupt_maxw, &b_mmc_taupt_maxw);
   fChain->SetBranchAddress("mmc_tauphi_maxw", &mmc_tauphi_maxw, &b_mmc_tauphi_maxw);
   fChain->SetBranchAddress("mmc_taueta_maxw", &mmc_taueta_maxw, &b_mmc_taueta_maxw);
   fChain->SetBranchAddress("mmc_taumass_maxw", &mmc_taumass_maxw, &b_mmc_taumass_maxw);
   fChain->SetBranchAddress("mmc_mass_mlnu3m", &mmc_mass_mlnu3m, &b_mmc_mass_mlnu3m);
   fChain->SetBranchAddress("mmc_taupt_mlnu3m", &mmc_taupt_mlnu3m, &b_mmc_taupt_mlnu3m);
   fChain->SetBranchAddress("mmc_tauphi_mlnu3m", &mmc_tauphi_mlnu3m, &b_mmc_tauphi_mlnu3m);
   fChain->SetBranchAddress("mmc_taueta_mlnu3m", &mmc_taueta_mlnu3m, &b_mmc_taueta_mlnu3m);
   fChain->SetBranchAddress("mmc_taumass_mlnu3m", &mmc_taumass_mlnu3m, &b_mmc_taumass_mlnu3m);
   fChain->SetBranchAddress("JetEta", &JetEta, &b_JetEta);
   fChain->SetBranchAddress("JetPhi", &JetPhi, &b_JetPhi);
   fChain->SetBranchAddress("JetPt", &JetPt, &b_JetPt);
   fChain->SetBranchAddress("JetE", &JetE, &b_JetE);
   fChain->SetBranchAddress("Jet_btag", &Jet_btag, &b_Jet_btag);
   fChain->SetBranchAddress("Jet_btSF", &Jet_btSF, &b_Jet_btSF);
   fChain->SetBranchAddress("Jet_JVT", &Jet_JVT, &b_Jet_JVT);
   fChain->SetBranchAddress("n_bjets", &n_bjets, &b_n_bjets);
   fChain->SetBranchAddress("truth_n_bjets", &truth_n_bjets, &b_truth_n_bjets);
   fChain->SetBranchAddress("MuonEta", &MuonEta, &b_MuonEta);
   fChain->SetBranchAddress("MuonPhi", &MuonPhi, &b_MuonPhi);
   fChain->SetBranchAddress("MuonPt", &MuonPt, &b_MuonPt);
   fChain->SetBranchAddress("MuonE", &MuonE, &b_MuonE);
   fChain->SetBranchAddress("MuonCharge", &MuonCharge, &b_MuonCharge);
   fChain->SetBranchAddress("Muon_d0sig", &Muon_d0sig, &b_Muon_d0sig);
   fChain->SetBranchAddress("Muon_delta_z0", &Muon_delta_z0, &b_Muon_delta_z0);
   fChain->SetBranchAddress("Muon_recoSF", &Muon_recoSF, &b_Muon_recoSF);
   fChain->SetBranchAddress("Muon_isoSF", &Muon_isoSF, &b_Muon_isoSF);
   fChain->SetBranchAddress("MuonMatchedTriggers", &MuonMatchedTriggers, &b_MuonMatchedTriggers);
   fChain->SetBranchAddress("Muon_ttvaSF", &Muon_ttvaSF, &b_Muon_ttvaSF);
   fChain->SetBranchAddress("EleEta", &EleEta, &b_EleEta);
   fChain->SetBranchAddress("ElePhi", &ElePhi, &b_ElePhi);
   fChain->SetBranchAddress("ElePt", &ElePt, &b_ElePt);
   fChain->SetBranchAddress("EleE", &EleE, &b_EleE);
   fChain->SetBranchAddress("EleCharge", &EleCharge, &b_EleCharge);
   fChain->SetBranchAddress("Ele_d0sig", &Ele_d0sig, &b_Ele_d0sig);
   fChain->SetBranchAddress("Ele_delta_z0", &Ele_delta_z0, &b_Ele_delta_z0);
   fChain->SetBranchAddress("Ele_recoSF", &Ele_recoSF, &b_Ele_recoSF);
   fChain->SetBranchAddress("Ele_idSF", &Ele_idSF, &b_Ele_idSF);
   fChain->SetBranchAddress("Ele_isoSF", &Ele_isoSF, &b_Ele_isoSF);
   fChain->SetBranchAddress("Ele_trigSF", &Ele_trigSF, &b_Ele_trigSF);
   fChain->SetBranchAddress("EleMatchedTriggers", &EleMatchedTriggers, &b_EleMatchedTriggers);
   //fChain->SetBranchAddress("PhotonEta", &PhotonEta, &b_PhotonEta);
   //fChain->SetBranchAddress("PhotonPhi", &PhotonPhi, &b_PhotonPhi);
   //fChain->SetBranchAddress("PhotonPt", &PhotonPt, &b_PhotonPt);
   //fChain->SetBranchAddress("PhotonE", &PhotonE, &b_PhotonE);
   fChain->SetBranchAddress("TruthJetEta", &TruthJetEta, &b_TruthJetEta);
   fChain->SetBranchAddress("TruthJetPhi", &TruthJetPhi, &b_TruthJetPhi);
   fChain->SetBranchAddress("TruthJetPt", &TruthJetPt, &b_TruthJetPt);
   fChain->SetBranchAddress("TruthJetE", &TruthJetE, &b_TruthJetE);
   fChain->SetBranchAddress("TruthNeutrinoEta", &TruthNeutrinoEta, &b_TruthNeutrinoEta);
   fChain->SetBranchAddress("TruthNeutrinoPhi", &TruthNeutrinoPhi, &b_TruthNeutrinoPhi);
   fChain->SetBranchAddress("TruthNeutrinoPt", &TruthNeutrinoPt, &b_TruthNeutrinoPt);
   fChain->SetBranchAddress("TruthNeutrinoE", &TruthNeutrinoE, &b_TruthNeutrinoE);
   fChain->SetBranchAddress("TruthNeutrinoCharge", &TruthNeutrinoCharge, &b_TruthNeutrinoCharge);
   fChain->SetBranchAddress("TruthMuonEta", &TruthMuonEta, &b_TruthMuonEta);
   fChain->SetBranchAddress("TruthMuonPhi", &TruthMuonPhi, &b_TruthMuonPhi);
   fChain->SetBranchAddress("TruthMuonPt", &TruthMuonPt, &b_TruthMuonPt);
   fChain->SetBranchAddress("TruthMuonE", &TruthMuonE, &b_TruthMuonE);
   fChain->SetBranchAddress("TruthMuonCharge", &TruthMuonCharge, &b_TruthMuonCharge);
   fChain->SetBranchAddress("TruthEleEta", &TruthEleEta, &b_TruthEleEta);
   fChain->SetBranchAddress("TruthElePhi", &TruthElePhi, &b_TruthElePhi);
   fChain->SetBranchAddress("TruthElePt", &TruthElePt, &b_TruthElePt);
   fChain->SetBranchAddress("TruthEleE", &TruthEleE, &b_TruthEleE);
   fChain->SetBranchAddress("TruthEleCharge", &TruthEleCharge, &b_TruthEleCharge);
   fChain->SetBranchAddress("TruthTauEta", &TruthTauEta, &b_TruthTauEta);
   fChain->SetBranchAddress("TruthTauPhi", &TruthTauPhi, &b_TruthTauPhi);
   fChain->SetBranchAddress("TruthTauPt", &TruthTauPt, &b_TruthTauPt);
   fChain->SetBranchAddress("TruthTauM", &TruthTauM, &b_TruthTauM);
   fChain->SetBranchAddress("VisTruthTauEta", &VisTruthTauEta, &b_VisTruthTauEta);
   fChain->SetBranchAddress("VisTruthTauPhi", &VisTruthTauPhi, &b_VisTruthTauPhi);
   fChain->SetBranchAddress("VisTruthTauPt", &VisTruthTauPt, &b_VisTruthTauPt);
   fChain->SetBranchAddress("VisTruthTauM", &VisTruthTauM, &b_VisTruthTauM);
   fChain->SetBranchAddress("TruthTauCharge", &TruthTauCharge, &b_TruthTauCharge);
   fChain->SetBranchAddress("TruthTau_isHadronic", &TruthTau_isHadronic, &b_TruthTau_isHadronic);
   fChain->SetBranchAddress("TruthTau_decay_mode", &TruthTau_decay_mode, &b_TruthTau_decay_mode);
   fChain->SetBranchAddress("MET_etx", &MET_etx, &b_MET_etx);
   fChain->SetBranchAddress("MET_ety", &MET_ety, &b_MET_ety);
   fChain->SetBranchAddress("MET_met", &MET_met, &b_MET_met);
   fChain->SetBranchAddress("MET_phi", &MET_phi, &b_MET_phi);
   //fChain->SetBranchAddress("MDN_tautau_mll", &MDN_tautau_mll, &b_MDN_tautau_mll);
   //fChain->SetBranchAddress("MDN_tautau_unc", &MDN_tautau_unc, &b_MDN_tautau_unc);
   //fChain->SetBranchAddress("DNN_mll", &DNN_mll, &b_DNN_mll);
   //fChain->SetBranchAddress("DSR_mll", &DSR_mll, &b_DSR_mll);
   //fChain->SetBranchAddress("DNN_score", &DNN_score, &b_DNN_score);
   //fChain->SetBranchAddress("DNN_outputs", &DNN_outputs, &b_DNN_outputs);
   //fChain->SetBranchAddress("MDN_outputs", &MDN_outputs, &b_MDN_outputs);
   //fChain->SetBranchAddress("DS_outputs", &DS_outputs, &b_DS_outputs);
   //fChain->SetBranchAddress("MDN_covariance", &MDN_covariance, &b_MDN_covariance);
   fChain->SetBranchAddress("smim_mass", &smim_mass, &b_smim_mass);
   fChain->SetBranchAddress("smim_ditaupt", &smim_ditaupt, &b_smim_ditaupt);
   fChain->SetBranchAddress("smim_tau1pt", &smim_tau1pt, &b_smim_tau1pt);
   fChain->SetBranchAddress("smim_tau2pt", &smim_tau2pt, &b_smim_tau2pt);
   fChain->SetBranchAddress("smim_tau1eta", &smim_tau1eta, &b_smim_tau1eta);
   fChain->SetBranchAddress("smim_tau2eta", &smim_tau2eta, &b_smim_tau2eta);
   fChain->SetBranchAddress("MT2", &MT2, &b_MT2);
   fChain->SetBranchAddress("ttm_mass", &ttm_mass, &b_ttm_mass);
   fChain->SetBranchAddress("collinear", &collinear, &b_collinear);
   fChain->SetBranchAddress("PassedTriggers", &PassedTriggers, &b_PassedTriggers);
   fChain->SetBranchAddress("TriggerSF", &TriggerSF, &b_TriggerSF);

   Notify();
}

Bool_t CLoop::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.
   std::cout << "Branches loaded fine" << std::endl;
   return kTRUE;
}

void CLoop::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}

Int_t CLoop::Cut(/*Long64_t entry*/)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}

void CLoop::ActivateBranches(const std::string& key){
    // Only activate relevant branches
}

void CLoop::createOutputFile(const std::string& key){
    const char*  name_root = key.c_str();
    // Create output file
    m_outputFile = std::make_unique<TFile>(name_root,"recreate");
}