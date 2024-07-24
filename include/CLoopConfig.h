/*
@brief Config structure for CLoop. Stores options like if the user wants to save histograms or events.
*/

#pragma once

#include <string>

enum MassRegion
{
    LOW_M,
    MEDIUM_M,
    HIGH_M,
};

struct CLoopConfig
{
    CLoopConfig() = default;
    CLoopConfig(bool saveHistograms, bool saveEvents, bool reweightMjj, std::string bdtWeightsPath, std::string region, int massRegion) : 
        m_saveHistograms{saveHistograms}, m_saveEvents{saveEvents}, m_reweightMjj{reweightMjj}, m_bdtWeightsPath{bdtWeightsPath}, m_region{region},
        m_massRegion{massRegion} 
        {}
    ~CLoopConfig() = default;
    
    bool m_saveHistograms{true};
    bool m_saveEvents{false};
    bool m_reweightMjj{true};
    std::string m_bdtWeightsPath{""};
    std::string m_region{"all"};
    int m_massRegion{LOW_M};
};