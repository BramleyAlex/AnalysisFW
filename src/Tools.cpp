// Header file with common tools used in analysis.
#include <vector>
#include <TLorentzVector.h>
#include <sstream>
#include <string>
#include <memory>
#include <cmath>
#include <utility>

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

// Function to calculate the pT balance of the jets and taus.
//@param particles: vector of particle four momenta
double CalculatePtBalance(const std::vector<TLorentzVector>& particles)
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
    if (angle_tau0_met <= angle_tau1_met)
    {
      double omega = angle_tau0_met / angle_tau_tau;
      if (angle_tau0_met + angle_tau1_met <= (angle_tau_tau+0.00000001)) return omega;
      else return omega * -1;
    }
    else
    {
      double omega = angle_tau1_met / angle_tau_tau;
      if (angle_tau0_met + angle_tau1_met <= (angle_tau_tau+0.00000001)) return 1 - omega;
      else return omega + 1;
    }
}

std::pair<TLorentzVector,TLorentzVector> CalculateNeutrinoVector(const TLorentzVector& met_p4, const TLorentzVector& tau_0_p4, const TLorentzVector& tau_1_p4)
{
  double omega = CalculateOmega(tau_0_p4, tau_1_p4, met_p4);
  TLorentzVector neu_0_p4;
  TLorentzVector neu_1_p4;

  TVector3 tau_0_dir = tau_0_p4.Vect();
  tau_0_dir.SetZ(0);
  tau_0_dir = tau_0_dir.Unit();
  TVector3 tau_1_dir = tau_1_p4.Vect();
  tau_1_dir.SetZ(0);
  tau_1_dir = tau_1_dir.Unit();

  if (omega >= 0 && omega <= 1)
  {
    double numerator_0 = met_p4.Vect().Dot((tau_0_dir - (tau_0_dir.Dot(tau_1_dir)) * tau_1_dir));
    double denominator = 1 - std::pow(tau_0_dir.Dot(tau_1_dir),2);
    double neu_0_pt = numerator_0/denominator;
    neu_0_p4.SetPtEtaPhiM(neu_0_pt,tau_0_p4.Eta(),tau_0_p4.Phi(),0);

    double numerator_1 = met_p4.Vect().Dot((tau_1_dir - (tau_0_dir.Dot(tau_1_dir)) * tau_0_dir));
    double neu_1_pt = numerator_1/denominator;
    neu_1_p4.SetPtEtaPhiM(neu_1_pt,tau_1_p4.Eta(),tau_1_p4.Phi(),0);
  }
  else if (omega < 0)
  {
    double neu_0_pt = met_p4.Vect() * tau_0_dir;
    if (neu_0_pt >=0) neu_0_p4.SetPtEtaPhiM(neu_0_pt,tau_0_p4.Eta(),tau_0_p4.Phi(),0);
    else neu_0_p4.SetPtEtaPhiM(-1*neu_0_pt,tau_0_p4.Eta(),tau_0_p4.Phi()+TMath::Pi(),0);
    neu_1_p4.SetXYZM(0, 0, 0, 0);
  }
  else
  {
    double neu_1_pt = met_p4.Vect() * tau_1_dir;
    neu_0_p4.SetXYZM(0, 0, 0, 0);
    if (neu_1_pt >=0) neu_1_p4.SetPtEtaPhiM(neu_1_pt,tau_1_p4.Eta(),tau_1_p4.Phi(),0);
    else neu_1_p4.SetPtEtaPhiM(-1*neu_1_pt,tau_1_p4.Eta(),tau_1_p4.Phi()+TMath::Pi(),0);
  }
  return std::pair<TLorentzVector,TLorentzVector> (neu_0_p4,neu_1_p4);
}

void PrintTLorentzVectorPtEtaPhi(const TLorentzVector& input)
{
  std::cout<<"Pt: "<<input.Pt()<<" Eta: "<<input.Eta()<<" Phi: "<<input.Phi()<< " E: "<<input.E()<<std::endl;
}

void PrintTLorentzVectorXYZT(const TLorentzVector& input)
{
  std::cout<<"X: "<<input.X()<<" Y: "<<input.Y()<<" Z: "<<input.Z()<< " T: "<<input.E()<<std::endl;
}

bool CalculateNGapJets(const double &ljet_0_rapidity, const double &ljet_1_rapidity, const std::vector<float> *JetEta)
{
  if (JetEta->size()>=3)
  {
    double lower_bound;
    double upper_bound;
    if (ljet_0_rapidity < ljet_1_rapidity)
    {
      lower_bound = ljet_0_rapidity;
      upper_bound = ljet_1_rapidity;
    }
    else
    {
      lower_bound = ljet_1_rapidity;
      upper_bound = ljet_0_rapidity;
    }
    double ljet_2_rapidity = JetEta->at(2);
    if (ljet_2_rapidity >= lower_bound && ljet_2_rapidity <= upper_bound) 
    {
      return true;
    }
  }
  return false;
}

bool Region(std::string region_id, const bool& N_gap_jets, const double& Z_centrality)
{
  if (region_id == "SR") {if (N_gap_jets == 0 && Z_centrality < 0.5) return true;}
  else if (region_id == "CRa") {if (N_gap_jets == 1 && Z_centrality < 0.5) return true;}
  else if (region_id == "CRb") {if (N_gap_jets == 1 && Z_centrality >= 0.5 && Z_centrality <=1) return true;}
  else if (region_id == "CRc") {if (N_gap_jets == 0 && Z_centrality >= 0.5 && Z_centrality <=1) return true;}
  else if (region_id == "CR") {if (Z_centrality >= 0.5 || N_gap_jets == 1) return true;}
  else if (region_id == "all") {return true;}
  return false;
}
