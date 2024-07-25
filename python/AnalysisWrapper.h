
/**
 * @brief Wrapper classes are defined here to load them in python
 * Wrapper cannot return references or custom classes.
 */
#include <memory>
#include <TTree.h>
#include <string>
#include "CLoop.h"
#include "CLoopConfig.h"


class CLoopWrapper {
    public:
    explicit CLoopWrapper(long long unsigned int TTreePointer, std::string sampleName) :
        m_cloop(std::make_shared<CLoop>(reinterpret_cast<TTree*>(TTreePointer), sampleName)) {};

    ~CLoopWrapper() = default;

    void Loop(float lumFactor, int z_sample, std::string key, CLoopConfig config){
        m_cloop->Loop(lumFactor, z_sample, key, config);
    }

    private:
    std::shared_ptr<CLoop> m_cloop;
};
