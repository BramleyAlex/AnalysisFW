#pragma once
#include <TMVA/Reader.h>
#include <memory>
#include <string>
// Class to manage the BDT
class VBFBDT {
  public:
    VBFBDT() = default;
    
    VBFBDT(const std::string& weightsFilePath) {
        m_weightsFilePath = weightsFilePath.c_str();
        m_reader = std::make_unique<TMVA::Reader>("Silent");
        m_reader->AddVariable("mjj",&m_bdt_mjj);
        m_reader->AddVariable("JetGap",&m_bdt_drap);
        m_reader->AddVariable("PtBalance",&m_bdt_ptbal);
        m_reader->AddVariable("centrality",&m_bdt_zcen);
        m_reader->AddSpectator("eventNumber", &m_bdt_eventNumber); // For deterministic split
        m_reader->BookMVA("VBF_BDT", weightsFilePath.c_str());
    }

    ~VBFBDT() {}

    VBFBDT& operator=(const VBFBDT& other) {
        if (this != &other) {
            m_reader = std::make_unique<TMVA::Reader>("Silent");
            m_bdt_mjj = other.m_bdt_mjj;
            m_bdt_drap = other.m_bdt_drap;
            m_bdt_ptbal = other.m_bdt_ptbal;
            m_bdt_zcen = other.m_bdt_zcen;
            m_bdt_eventNumber = other.m_bdt_eventNumber;
            m_weightsFilePath = other.m_weightsFilePath;

            m_reader->AddVariable("mjj",&m_bdt_mjj);
            m_reader->AddVariable("JetGap",&m_bdt_drap);
            m_reader->AddVariable("PtBalance",&m_bdt_ptbal);
            m_reader->AddVariable("centrality",&m_bdt_zcen);
            m_reader->AddSpectator("eventNumber", &m_bdt_eventNumber); // For deterministic split
            m_reader->BookMVA("VBF_BDT",  other.m_weightsFilePath);
        }
        return *this;
    }

    void update(float mjj, float drap, float ptbal, float zcen, long eventNumber) {
        m_bdt_mjj = mjj;
        m_bdt_drap = drap;
        m_bdt_ptbal = ptbal;
        m_bdt_zcen = zcen;
        m_bdt_eventNumber = eventNumber % 10;
    }

    double evaluate() {
        double bdtScore = m_reader->EvaluateMVA("VBF_BDT");
        reset();
        return bdtScore;
    }

    void reset() {
        m_bdt_mjj = 0;
        m_bdt_drap = 0;
        m_bdt_ptbal = 0;
        m_bdt_zcen = 0;
        m_bdt_eventNumber = 0;
    }

  private:
    std::unique_ptr<TMVA::Reader> m_reader;
    float m_bdt_mjj;
    float m_bdt_drap;
    float m_bdt_ptbal;
    float m_bdt_zcen;
    int m_bdt_eventNumber;
    const char* m_weightsFilePath;
};