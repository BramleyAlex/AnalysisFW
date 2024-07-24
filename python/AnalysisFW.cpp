#include <iostream>
#include <string>
#include "CLoop.h"
#include "TFile.h"
#include "TTree.h"

#include <boost/python.hpp>
#include "AnalysisWrapper.h"
#include "CLoopConfig.h"
using namespace boost::python;

BOOST_PYTHON_MODULE(AnalysisFW)
{
    class_<CLoopWrapper>("CLoop", init<long long unsigned int, std::string>())
        .def("Loop", &CLoopWrapper::Loop)
        .def("GetNEntries", &CLoopWrapper::GetNEntries)
    ;

    class_<CLoopConfig>("CLoopConfig", init<bool, bool, bool, std::string, std::string, int>())
        .def_readwrite("m_saveHistograms", &CLoopConfig::m_saveHistograms)
        .def_readwrite("m_saveEvents", &CLoopConfig::m_saveEvents)
        .def_readwrite("m_reweightMjj", &CLoopConfig::m_reweightMjj)
        .def_readwrite("m_bdtWeightsPath", &CLoopConfig::m_bdtWeightsPath)
        .def_readwrite("m_region", &CLoopConfig::m_region)
        .def_readwrite("m_massRegion", &CLoopConfig::m_massRegion)
        .enable_pickling()
    ;
}