#pragma once
#include <TTree.h>
#include <memory>

class OutputTree {
    public:
    OutputTree() = default;

    OutputTree(const char* treeName, const char* treeDescription);

    ~OutputTree(){};

    void FillTree(){
        m_Tree->Fill();
    }

    const TTree* GetTree() const {
        return m_Tree.get();
    }

    private:
    std::unique_ptr<TTree> m_Tree = nullptr;
    public: // To be able to access this members directly and assign in FillTree.
    double m_mcWeight;
    double m_mjj;
    double m_deltaPhiLT;
    double m_jetRNNScore;
    double m_transverseMassLep;
    double m_massTauLep;
    double m_tau_pT;
    double m_lep_pT;
    double m_jet0_pT;
    double m_jet1_pT;
    double m_met_pT;
    double m_event_number;
    double m_pt_balance;
    double m_jet_gap;
};