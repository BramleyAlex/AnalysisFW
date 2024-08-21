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
    class_<CLoopConfig>("CLoopConfig", init<bool, bool, bool, std::string, std::string, int>())
        .def_readwrite("m_saveHistograms", &CLoopConfig::m_saveHistograms)
        .def_readwrite("m_saveEvents", &CLoopConfig::m_saveEvents)
        .def_readwrite("m_reweightMjj", &CLoopConfig::m_reweightMjj)
        .def_readwrite("m_bdtWeightsPath", &CLoopConfig::m_bdtWeightsPath)
        .def_readwrite("m_region", &CLoopConfig::m_region)
        .def_readwrite("m_massRegion", &CLoopConfig::m_massRegion)
        .enable_pickling()
    ;

    class_<VBFWrapper>("LM", init<long long unsigned int, std::string>())
        .def("Loop", &VBFWrapper::Loop)
    ;

    class_<MJWrapper>("MJ", init<long long unsigned int, std::string>())
        .def("Loop", &MJWrapper::Loop)
    ;

    class_<HMWrapper>("HM", init<long long unsigned int, std::string>())
        .def("Loop", &HMWrapper::Loop)
    ;

    class_<BDTWrapper>("BDT", init<long long unsigned int, std::string>())
        .def("Loop", &BDTWrapper::Loop)
    ;

    class_<BDT_CutWrapper>("BDT_Cut", init<long long unsigned int, std::string>())
        .def("Loop", &BDT_CutWrapper::Loop)
    ;
}