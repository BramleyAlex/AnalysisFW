#pragma once

double pi=TMath::Pi();
std::vector<std::string> cutNames{"basic","t0pt","t1pt","j0pt","j1pt","mtt","yjj","ome","ptbal","rgn","mjj"};
std::vector<std::string> notFull{"basic","all"};

// Raw historgrams
TH1F* nJets = new TH1F("nJets","Number of jets",10,0,10);
TH1F* tau0Eta = new TH1F("tau0Eta","Leading Tau Eta",60,-3,3);

// Histogram containers
histogramContainer tau0_ptContainer{"tau0_pt","Leading visible Tau pT",500,0,500,cutNames,"t0pt"};
histogramContainer tau1_ptContainer{"tau1_pt","Sub-leading visible Tau pT",500,0,500,cutNames,"t1pt"};
histogramContainer n_bjetsContainer{"n_bjets","Number of b_jets",5,0,5,notFull};
histogramContainer delta_phiContainer{"delta_phi","Delta phi between taus",32,0,3.2,cutNames};
histogramContainer mass_jjContainer{"mass_jj","Invariant mass di_jet system",5000,0,5000,cutNames,"mjj"};
histogramContainer ljet0_ptContainer{"ljet0_pt","Light jet0 pT",1000,0,1000,cutNames,"j0pt"};
histogramContainer ljet1_ptContainer{"ljet1_pt","Light jet1 pT",1000,0,1000,cutNames,"j1pt"};
histogramContainer bdtContainer{"bdtScore","BDT Score",200,-1,1,notFull};
histogramContainer visibleMassContainer{"visibleMass","Visible mass tau-tau",160,0,160,cutNames};
histogramContainer tau1TransMassContainer{"tau1TransMass","Transverse mass sub-leading Tau",500,0,500,cutNames};
histogramContainer rnn_score_1pContainer{"rnn_score_1p","RNN Score 1 prong taus",100,0,1,notFull};
histogramContainer rnn_score_3pContainer{"rnn_score_3p","RNN Score 3 prong taus",100,0,1,notFull};
histogramContainer Z_ptContainer{"Z_pt","ZpT visible",1000,0,1000,notFull};
histogramContainer tau_matched_1pContainer{"tau_matched_1p","Tau truth matched 1 prong",2,0,2,notFull};
histogramContainer tau_matched_3pContainer{"tau_matched_3p","Tau truth matched 3 prong",2,0,2,notFull};
histogramContainer delta_R_lepjetContainer{"delta_R_tau1jet","Delta R tau1-jet",60,0,6,cutNames};
histogramContainer delta_R_taujetContainer{"delta_R_tau0jet","Delta R tau0-jet",60,0,6,cutNames};
histogramContainer pt_balanceContainer{"pt_balance","pT Balance",100,0,1,cutNames,"ptbal"};
histogramContainer jet_gapContainer{"jet_gap", "Jet Gap",100,0,10,cutNames,"yjj"};
histogramContainer omegaContainer{"omega", "Omega (MET angle parameter)",600,-3,3,cutNames,"ome"};

histogramContainer neu0_pt_insideContainer{"neu0_pt_in", "Neutrino 0 pT (MET inside)", 500,0,500,notFull};
histogramContainer neu1_pt_insideContainer{"neu1_pt_in", "Neutrino 1 pT (MET inside)", 500,0,500,notFull};
histogramContainer neu0_pt_outsideContainer{"neu0_pt_out", "Neutrino 0 pT (MET outside)", 500,0,500,notFull};
histogramContainer neu1_pt_outsideContainer{"neu1_pt_out", "Neutrino 1 pT (MET outside)", 500,0,500,notFull};

histogramContainer tau0_reco_pt_inContainer{"tau0_reco_pt_in","Leading reconstructed Tau pT (MET inside)",500,0,500,notFull};
histogramContainer tau1_reco_pt_inContainer{"tau1_reco_pt_in","Sub-leading reconstructed Tau pT (MET inside)",500,0,500,notFull};
histogramContainer tau0_reco_pt_outContainer{"tau0_reco_pt_out","Leading reconstructed Tau pT (MET outside)",500,0,500,notFull};
histogramContainer tau1_reco_pt_outContainer{"tau1_reco_pt_out","Sub-leading reconstructed Tau pT (MET outside)",500,0,500,notFull};

histogramContainer neutrino_pt_testContainer{"neutrino_pt_test","difference between neutrinos and MET pt",1000,-500,500,notFull};
histogramContainer neutrino_phi_testContainer{"neutrino_phi_test","difference between neutrinos and MET phi",64,-3.2,3.2,notFull};
histogramContainer neutrino_eta_testContainer{"neutrino_eta_test","difference between neutrinos and MET eta",500,-5,5,notFull};
histogramContainer neutrino_px_testContainer{"neutrino_px_test", "difference between neutrinos and MET px",1000,-500,500,notFull};
histogramContainer neutrino_py_testContainer{"neutrino_py_test", "difference between neutrinos and MET py",1000,-500,500,notFull};

histogramContainer reco_mass_tt_insideContainer{"reco_mass_tt_in","Invariant of mass reconstructed tau-tau (MET inside)",1000,0,1000,cutNames,"mtt"};
histogramContainer reco_mass_tt_outsideContainer{"reco_mass_tt_out","Invariant of mass reconstructed tau-tau (MET outside)",1000,0,1000,cutNames,"mtt"};
histogramContainer reco_mass_ttContainer{"reco_mass_tt","Invariant of mass reconstructed tau-tau",1000,0,1000,cutNames,"mtt"};

histogramContainer met_angleContainer{"met_angle","Angle of MET with closest tau",32,0,3.2,notFull};

histogramContainer z_centralityContainer{"centrality", "Z centrality",40,0,2,cutNames,"rgn"};
histogramContainer N_gapjetsContainer{"N_gapjets", "Number of gap jets",2,0,2,cutNames,"rgn"};

histogramContainer z_sampleContainer{"z_sample", "z sample",5,0,5,notFull};
histogramContainer weightContainer{"weight", "weight",200,-10,10,notFull};

histogramContainer jet_0_phiContainer{"jet_0_phi", "Leading Jet Phi",64,-3.2,3.2,cutNames};
histogramContainer jet_1_phiContainer{"jet_1_phi", "Sub-Leading Jet Phi",64,-3.2,3.2,cutNames};
histogramContainer jet_delta_phiContainer{"jet_delta_phi", "Angle between Jets",64,0,3.2,cutNames};

histogramContainer jet_0_etaContainer{"jet_0_eta", "Leading Jet Pseudo-Rapidity",60,-3,3,cutNames};
histogramContainer jet_1_etaContainer{"jet_1_eta", "Sub-Leading Jet Pseudo-Rapidity",60,-3,3,cutNames};
histogramContainer jet_delta_etaContainer{"jet_delta_eta", "Pseudo-Rapidity between jets",60,0,6,cutNames};

histogramContainer jet_0_EContainer{"jet_0_E", "Leading Jet Energy",1000,0,1000,cutNames};
histogramContainer jet_1_EContainer{"jet_1_E", "Sub-Leading Jet Energy",1000,0,1000,cutNames};