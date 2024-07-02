#pragma once

double pi=TMath::Pi();
std::vector<std::string> cutNames{"basic","t1pt","t0pt","j0pt","j1pt","mjj","mtt","yjj"};
std::vector<std::string> notFull{"basic","all"};

// Raw historgrams
TH1F* nJets = new TH1F("nJets","Number of jets",10,0,10);
TH1F* tau0Eta = new TH1F("tau0Eta","Leading Tau Eta",60,-3,3);

// Histogram containers
histogramContainer tau1_ptContainer{"tau1_pt","Sub-leading Tau pT",500,0,500,cutNames,"t1pt"};
histogramContainer tau0_ptContainer{"tau0_pt","Leading Tau pT",500,0,500,cutNames,"t0pt"};
histogramContainer n_bjetsContainer{"n_bjets","Number of b_jets",5,0,5,notFull};
histogramContainer delta_phiContainer{"delta_phi","Delta phi between taus",32,0,3.2,cutNames};
histogramContainer mass_jjContainer{"mass_jj","Invariant mass di_jet system",50,0,5000,cutNames,"mjj"};
histogramContainer mass_ttContainer{"mass_tt","Invariant of mass tau-tau system",500,0,1500,cutNames};
histogramContainer ljet0_ptContainer{"ljet0_pt","Light jet0 pT",1000,0,1000,cutNames,"j0pt"};
histogramContainer ljet1_ptContainer{"ljet1_pt","Light jet1 pT",1000,0,1000,cutNames,"j1pt"};
histogramContainer bdtContainer{"bdtScore","BDT Score",200,-1,1,notFull};
histogramContainer visibleMassContainer{"visibleMass","Visible mass tau-tau",1000,0,1000,cutNames,"mtt"};
histogramContainer tau1TransMassContainer{"tau1TransMass","Transverse mass sub-leading Tau",500,0,500,cutNames};
histogramContainer rnn_score_1pContainer{"rnn_score_1p","RNN Score 1 prong taus",100,0,1,notFull};
histogramContainer rnn_score_3pContainer{"rnn_score_3p","RNN Score 3 prong taus",100,0,1,notFull};
histogramContainer Z_ptContainer{"Z_pt","ZpT visible",1000,0,1000,notFull};
histogramContainer tau_matched_1pContainer{"tau_matched_1p","Tau truth matched 1 prong",2,0,2,notFull};
histogramContainer tau_matched_3pContainer{"tau_matched_3p","Tau truth matched 3 prong",2,0,2,notFull};
histogramContainer delta_R_lepjetContainer{"delta_R_tau1jet","Delta R tau1-jet",60,0,6,cutNames};
histogramContainer delta_R_taujetContainer{"delta_R_tau0jet","Delta R tau0-jet",60,0,6,cutNames};
histogramContainer pt_balanceContainer{"pt_balance","pT Balance",50,0,1,cutNames};
histogramContainer jet_gapContainer{"jet_gap", "Jet Gap",50,0,10,cutNames,"yjj"};
histogramContainer omegaContainer{"omega", "Omega (MET angle parameter)",30,-3,3,cutNames};