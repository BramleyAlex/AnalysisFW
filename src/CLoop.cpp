#include <cmath>
#include <TMacro.h>
#include "ReweightingTools.h"
#include "CLoop.h"
#include "CLoopConfig.h"
#include <TLorentzVector.h>
#include <memory>

std::vector<std::string> split(const std::string& s, char delimiter);
TLorentzVector& toGeV(TLorentzVector &v);

void CLoop::Loop(float lumFactor, int z_sample, std::string key, const CLoopConfig& config)
{
    clock_t startTime = clock(); // get start time

    if (fChain == 0) return;

    Long64_t nentries = fChain->GetEntries();

    // if in fast mode only loop over 1% of the entries
    Long64_t nLoop = nentries;

    std::cout<<"Analysing "<<nLoop<<" Events!"<<std::endl;

    Long64_t nbytes = 0, nb = 0;

    ActivateBranches(key);
    // Create output file
    key = key+".root";
    createOutputFile(key);

    // Create BDT
    m_vbfBDT = VBFBDT(config.m_bdtWeightsPath);

    // Create TTree
    bool saveHistograms = config.m_saveHistograms;
    bool saveEvents = config.m_saveEvents;   
    // loop over number of entries
    int unordered_jet{0};
    int unordered_tau{0};
    for (Long64_t jentry=0; jentry<nLoop;jentry++) {
        Long64_t ientry = LoadTree(jentry);
        if (ientry < 0) break;
        nb = fChain->GetEntry(jentry,0);    nbytes += nb;

        // First, check that we have at least two jets and two taus
        if(TauPt->size() < 2 || JetPt->size() < 2) continue;
        // Jet vectors
        TLorentzVector ljet_0_p4;
        TLorentzVector ljet_1_p4;
        ljet_1_p4.SetPtEtaPhiE(JetPt->at(1),JetEta->at(1),JetPhi->at(1),JetE->at(1));
        ljet_0_p4.SetPtEtaPhiE(JetPt->at(0),JetEta->at(0),JetPhi->at(0),JetE->at(0));
        ljet_0_p4 = toGeV(ljet_0_p4);
        ljet_1_p4 = toGeV(ljet_1_p4);
        // Tau vectors
        TLorentzVector tau_0_p4;
        TLorentzVector tau_1_p4;
        tau_0_p4.SetPtEtaPhiE(TauPt->at(0),TauEta->at(0),TauPhi->at(0),TauE->at(0));
        tau_1_p4.SetPtEtaPhiE(TauPt->at(1),TauEta->at(1),TauPhi->at(1),TauE->at(1));
        tau_0_p4 = toGeV(tau_0_p4);
        tau_1_p4 = toGeV(tau_1_p4);

        if (ljet_0_p4.Pt() < ljet_1_p4.Pt()) unordered_jet += 1;
        if (tau_0_p4.Pt() < tau_1_p4.Pt()) unordered_jet += 1;

        // Variable defining regions
        // DELTA RAPIDITY 2-JETS
        double delta_y = abs(ljet_0_p4.Rapidity()-ljet_1_p4.Rapidity());
        // Z BOSON CENTRALITY
        double lepton_xi=(tau_0_p4+tau_1_p4).Rapidity();
        double dijet_xi=ljet_0_p4.Rapidity()+ljet_1_p4.Rapidity();
        double z_centrality=abs(lepton_xi-0.5*dijet_xi)/delta_y;

        Region region = Region::DefaultNoRW;
        if (z_centrality<0.5){region = Region::SR;}
        else if (z_centrality<=1.0){region = Region::CR;}

        double mjj = sqrt(2*(ljet_0_p4.Dot(ljet_1_p4)));
        double mjj_w = 1.0;

        // mjj reweighting
        bool reweight_mjj = config.m_reweightMjj;
        if (reweight_mjj){
            MC mcSample = static_cast<MC>(z_sample);
            if(mcSample == MC::PowHegPythia){
                mjj_w = 1.0;
            } else if (mcSample == MC::SHERPA){
                mjj_w = mjj_rw(mjj,parametersSHERPA[region]); 
            } else if (mcSample == MC::MadGraph){ 
                mjj_w = mjj_rw(mjj,parametersMadGraph[region]);
            }
        }
        double eventWeight = 1;
        // check if event is from real data
        if (!(key.substr(0,4)=="data")) {
            // take product of all scale factors
            eventWeight = weight*lumFactor*mjj_w;
    
        }

        // fill histograms
        //std::cout << eventWeight<< "   ";
        if (saveHistograms) Fill(eventWeight, z_sample, key);
        if (saveEvents) FillTree(eventWeight, z_sample, key);
        // end filling

    }
    // end style and writing
    if (unordered_jet > 0) std::cout<<"Warning: Sample contains unordered jets in "<<unordered_jet<<" events"<<std::endl;
    if (unordered_tau > 0) std::cout<<"Warning: Sample contains unordered jets in "<<unordered_tau<<" events"<<std::endl;


    if (saveHistograms) Style(lumFactor);   
    if (saveEvents) {
        m_outputFile->WriteObject(m_signalTree.GetTree(),"SIGNAL");
        m_outputFile->WriteObject(m_backgroundTree.GetTree(),"BACKGROUND");
    }
    // Add the code to the file if it is the first sample of the kind
    std::vector<std::string> tokens = split(key,'_');
    std::string lastToken = tokens[tokens.size()-1];
    bool isFirstSample = lastToken=="0NOMINAL.root";
    if (isFirstSample) {
        TMacro sourceCode("Analysis.cpp");
        sourceCode.Write();
    }

    clock_t endTime = clock(); // get end time
    // calculate time taken and print it
    double time_spent = (endTime - startTime) / CLOCKS_PER_SEC;
    std::cout << "Time processing == " <<time_spent << std::endl;
}

int CLoop::GetNEntries()
{
  Long64_t nentries = fChain->GetEntries();
  int n = static_cast<int>(nentries);
  return n;
}