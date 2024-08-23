#include "OutputTree.h"

OutputTree::OutputTree(const char* treeName, const char* treeDescription){
    // Creating the tree
    m_Tree = std::make_unique<TTree>(treeName, treeDescription);

    // Setting tree branches
    m_Tree->Branch("mcWeight", &m_mcWeight);
    m_Tree->Branch("mjj", &m_mjj);
    m_Tree->Branch("deltaPhiLT",&m_deltaPhiLT);
    m_Tree->Branch("t0RNNScore",&m_t0RNNScore);
    m_Tree->Branch("t1RNNScore",&m_t1RNNScore);
    m_Tree->Branch("transverseMassLep",&m_transverseMassLep);
    m_Tree->Branch("massTauTau",&m_massTauTau);
    m_Tree->Branch("tau0_pT", &m_tau0_pT);
    m_Tree->Branch("tau1_pT", &m_tau1_pT);
    m_Tree->Branch("jet0_pT", &m_jet0_pT);
    m_Tree->Branch("jet1_pT", &m_jet1_pT);
    m_Tree->Branch("met_pT", &m_met_pT);
    m_Tree->Branch("eventNumber", &m_event_number);
    m_Tree->Branch("PtBalance", &m_pt_balance);
    m_Tree->Branch("JetGap", &m_jet_gap);
    m_Tree->Branch("omega", &m_omega);
    m_Tree->Branch("centrality", &m_centrality);
    m_Tree->Branch("total_jet_pT", &m_total_jet_pT);
    m_Tree->Branch("total_tau_pT", &m_total_tau_pT);
    m_Tree->Branch("total_tau_mom", &m_total_tau_mom);
    m_Tree->Branch("transMassRatio", &m_transMassRatio);
    m_Tree->Branch("transMassRatioFunc", &m_transMassRatioFunc);
}