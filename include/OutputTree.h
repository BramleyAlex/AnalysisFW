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
    // Tutorial: https://github.com/diegobaronm/AnalysisFW/wiki/Mini%E2%80%90tutorial
    double m_metProjection;
    double m_mcWeight;
    double m_mjj;
    double m_deltaRapidity;
    double m_deltaPhiLT;
    double m_deltaPhiJJ;
    double m_jetRNNScore;
    double m_ptBalance;
    double m_zCentrality;
    double m_omega;
    double m_reco_mass;
    double m_lepNuPt;
    double m_transverseMassLep;
    double m_transverseRecoMassVariable;
    double m_massTauLep;
    int m_nLightJets;
    double m_tau_pT;
    double m_lep_pT;
    double m_jet0_pT;
    double m_jet1_pT;
    double m_met_pT;
    double m_event_number;
};