#pragma once

double pi=TMath::Pi();
// Tutorial: https://github.com/diegobaronm/AnalysisFW/wiki/Mini%E2%80%90tutorial
std::vector<std::string> cutNames{"basic","dphi","ptl","tpt","j1pt","mjj","mtl","metproj"};
std::vector<std::string> notFull{"basic","all"};

// Raw historgrams
TH1F* nJets = new TH1F("nJets","Number of jets",10,0,10);
TH1F* tauEta = new TH1F("tauEta","Tau Eta",60,-3,3);

// Tutorial: https://github.com/diegobaronm/AnalysisFW/wiki/Mini%E2%80%90tutorial
histogramContainer metProjecClosestLepContainer{"projec_closest_lep","MET projection to the closest lepton",100,-200,200,cutNames,"metproj"};

// Histogram containers
histogramContainer lep_ptContainer{"lep_pt","Lep pT",500,0,500,cutNames,"ptl"};
histogramContainer tau_ptContainer{"tau_pt","Tau pT",500,0,500,cutNames,"tpt"};
histogramContainer n_bjetsContainer{"n_bjets","Number of b_jets",5,0,5,notFull};
histogramContainer lepisoContainer{"lepiso","Lep Isolation",2,0,2,notFull};
histogramContainer delta_phiContainer{"delta_phi","Delta phi between tau and lep",32,0,3.2,cutNames,"dphi"};
histogramContainer mass_jjContainer{"mass_jj","Invariant mass di_jet system",5000,0,5000,cutNames,"mjj"};
histogramContainer ljet0_ptContainer{"ljet0_pt","Light jet0 pT",1000,0,1000,cutNames,"j1pt"};
histogramContainer ljet1_ptContainer{"ljet1_pt","Light jet1 pT",1000,0,1000,cutNames};
histogramContainer bdtContainer{"bdtScore","BDT Score",200,-1,1,notFull};
histogramContainer visibleMassContainer{"visibleMass","Visible mass tau-lep",1000,0,1000,cutNames,"mtl"};
histogramContainer lepTransMassContainer{"lepTransMass","Transverse mass lepton",500,0,500,cutNames};
histogramContainer rnn_score_1pContainer{"rnn_score_1p","RNN Score 1 prong taus",100,0,1,notFull};
histogramContainer rnn_score_3pContainer{"rnn_score_3p","RNN Score 3 prong taus",100,0,1,notFull};
histogramContainer Z_ptContainer{"Z_pt","ZpT visible",1000,0,1000,notFull};
histogramContainer tau_matched_1pContainer{"tau_matched_1p","Tau truth matched 1 prong",2,0,2,notFull};
histogramContainer tau_matched_3pContainer{"tau_matched_3p","Tau truth matched 3 prong",2,0,2,notFull};
histogramContainer delta_R_lepjetContainer{"delta_R_lepjet","Delta R lep-jet",60,0,6,cutNames};
histogramContainer delta_R_taujetContainer{"delta_R_taujet","Delta R tau-jet",60,0,6,cutNames};
