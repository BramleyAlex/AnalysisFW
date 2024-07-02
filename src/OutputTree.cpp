#include "OutputTree.h"

OutputTree::OutputTree(const char* treeName, const char* treeDescription){
    // Creating the tree
    m_Tree = std::make_unique<TTree>(treeName, treeDescription);

    // Setting tree branches
    m_Tree->Branch("mcWeight", &m_mcWeight);
    m_Tree->Branch("mjj", &m_mjj);
    m_Tree->Branch("deltaPhiLT",&m_deltaPhiLT);
    m_Tree->Branch("jetRNNScore",&m_jetRNNScore);
    m_Tree->Branch("transverseMassLep",&m_transverseMassLep);
    m_Tree->Branch("massTauLep",&m_massTauLep);
    m_Tree->Branch("tau_p4", &m_tau_pT);
    m_Tree->Branch("lep_p4", &m_lep_pT);
    m_Tree->Branch("jet0_p4", &m_jet0_pT);
    m_Tree->Branch("jet1_p4", &m_jet1_pT);
    m_Tree->Branch("met_p4", &m_met_pT);
    m_Tree->Branch("eventNumber", &m_event_number);
    m_Tree->Branch("PtBalance", &m_pt_balance);
    m_Tree->Branch("JetGap", &m_jet_gap);
}