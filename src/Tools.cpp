// Header file with common tools used in analysis.
#include <vector>
#include <TLorentzVector.h>
#include <sstream>
#include <string>
#include <memory>
#include <cmath>

// Function to split a string by a delimiter and return a vector of strings.
// Like Python's split function.
std::vector<std::string> split(const std::string& s, char delimiter) {
    std::vector<std::string> tokens;
    std::string token;
    std::istringstream tokenStream(s);
    while (std::getline(tokenStream, token, delimiter)) {
        tokens.push_back(token);
    }
    return tokens;
}

// Run numbers
const int run2015Begin = 276262;
const int run2015End   = 284484;

const int run2016Begin = 297730;
const int run2016End   = 311481;

const int run2017Begin = 323427;
const int run2017End   = 341649;

const int run2018Begin = 341649;
const int run2018End   = 364292;


// Function to calculate delta phi between two angles
// @param phi_1: angle 1
// @param phi_2: angle 2
double del_phi(double phi_1, double phi_2){
    double pi=TMath::Pi();
    double phi_1_norm, phi_2_norm;
    if (phi_1<0.0){
        phi_1_norm=phi_1+2*pi;
    }else {
        phi_1_norm=phi_1;
    }

    if (phi_2<0.0){
        phi_2_norm=phi_2+2*pi;
    }else {
        phi_2_norm=phi_2;
    }
    double delta=std::abs(phi_1_norm-phi_2_norm);
    if (delta>pi){
        delta=2*pi-delta;
        delta=std::abs(delta);
    }

    return delta;
}

// Function to calculate the minimum delta R between a test particle and a container of particles
// @param test_particle: particle to test
// @param bool_vector_container: container of booleans to select particles
// @param jet_container: container of particles to test against
double min_deltaR(const TLorentzVector& test_particle, const TLorentzVector& jet1, const TLorentzVector& jet2){

    double delta_R1=test_particle.DeltaR(jet1);
    double delta_R2=test_particle.DeltaR(jet2);

    double min_dR=std::min(delta_R1,delta_R2);
    return min_dR;
}

TLorentzVector& toGeV(TLorentzVector &v) {
    v.SetPtEtaPhiE(v.Pt()/1000., v.Eta(), v.Phi(), v.E()/1000.);
    return v;
}

// Function to calculate the pT balance of the jets and taus. NOTE - CURRENTLY MISSING NEUTRINO ENERGY
//@param particles: vector of particle four momenta
double CalculatePtBalance(const std::vector<TLorentzVector> particles)
{
    double vector_sum_x{0};
    double vector_sum_y{0};
    double scalar_sum{0};

    std::vector<TLorentzVector>::const_iterator particle_p4;
    for (particle_p4 = particles.begin(); particle_p4 < particles.end(); particle_p4++)
    {
        vector_sum_x += particle_p4->Px();
        vector_sum_y += particle_p4->Py();
        scalar_sum += particle_p4->Pt();
    }

    double vector_sum_mag = std::sqrt(std::pow(vector_sum_x,2) + std::pow(vector_sum_y,2));
    return vector_sum_mag / scalar_sum;
}

// Function to calculate omega parameter
double CalculateOmega(const TLorentzVector& tau_0_p4, const TLorentzVector& tau_1_p4, const TLorentzVector& met_p4)
{
    double angle_tau_tau = del_phi(tau_0_p4.Phi(), tau_1_p4.Phi());
    double angle_tau0_met = del_phi(tau_0_p4.Phi(), met_p4.Phi());
    double angle_tau1_met = del_phi(tau_1_p4.Phi(), met_p4.Phi());
    if (angle_tau0_met < angle_tau1_met)
    {
        double omega = angle_tau0_met / angle_tau_tau;
        if (angle_tau1_met < angle_tau_tau) return omega;
        else return omega * -1;
    }
    else
    {
        double omega = angle_tau1_met / angle_tau_tau;
        if (angle_tau0_met < angle_tau_tau) return 1 - omega;
        else return omega + 1;
    }
}