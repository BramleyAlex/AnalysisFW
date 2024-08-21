
/**
 * @brief Wrapper classes are defined here to load them in python
 * Wrapper cannot return references or custom classes.
 */
#include <memory>
#include <TTree.h>
#include <string>
#include "CLoop.h"
#include "CLoopConfig.h"

class VBFWrapper{
    public:
        explicit VBFWrapper(long long unsigned int TTreePointer, std::string sampleName)
        {m_cloop = std::make_shared<LowMass>(reinterpret_cast<TTree*>(TTreePointer), sampleName);}
        ~VBFWrapper() = default;

        void Loop(float lumFactor, int z_sample, std::string key, CLoopConfig config){
            m_cloop->Loop(lumFactor, z_sample, key, config);
        }
    protected:
        std::shared_ptr<CLoop> m_cloop;
};

class MJWrapper{
    public:
        explicit MJWrapper(long long unsigned int TTreePointer, std::string sampleName)
        {m_cloop = std::make_shared<MultiJetBG>(reinterpret_cast<TTree*>(TTreePointer), sampleName);}
        ~MJWrapper() = default;

    void Loop(float lumFactor, int z_sample, std::string key, CLoopConfig config){
        m_cloop->Loop(lumFactor, z_sample, key, config);
    }
    protected:
        std::shared_ptr<CLoop> m_cloop;
};

class HMWrapper{
    public:
        explicit HMWrapper(long long unsigned int TTreePointer, std::string sampleName)
        {m_cloop = std::make_shared<HighMass>(reinterpret_cast<TTree*>(TTreePointer), sampleName);}
        ~HMWrapper() = default;

    void Loop(float lumFactor, int z_sample, std::string key, CLoopConfig config){
        m_cloop->Loop(lumFactor, z_sample, key, config);
    }
    protected:
        std::shared_ptr<CLoop> m_cloop;
};

class BDTWrapper{
    public:
        explicit BDTWrapper(long long unsigned int TTreePointer, std::string sampleName)
        {m_cloop = std::make_shared<BDT_Sample>(reinterpret_cast<TTree*>(TTreePointer), sampleName);}
        ~BDTWrapper() = default;

    void Loop(float lumFactor, int z_sample, std::string key, CLoopConfig config){
        m_cloop->Loop(lumFactor, z_sample, key, config);
    }
    protected:
        std::shared_ptr<CLoop> m_cloop;
};

class BDT_CutWrapper{
    public:
        explicit BDT_CutWrapper(long long unsigned int TTreePointer, std::string sampleName)
        {m_cloop = std::make_shared<BDTCuts>(reinterpret_cast<TTree*>(TTreePointer), sampleName);}
        ~BDT_CutWrapper() = default;

    void Loop(float lumFactor, int z_sample, std::string key, CLoopConfig config){
        m_cloop->Loop(lumFactor, z_sample, key, config);
    }
    protected:
        std::shared_ptr<CLoop> m_cloop;
};
