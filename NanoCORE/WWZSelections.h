#ifndef WWZSELECTIONS_H
#define WWZSELECTIONS_H
#include "Nano.h"
#include "Base.h"


std::vector<pair<int,int>> getCandPairs_el(std::vector<int> el_input);
std::vector<pair<int,int>> getCandPairs_mu(std::vector<int> mu_input);
std::vector<pair<float,int>> get_el_resonance(std::vector<pair<int,int>> cand_el_pairs);
std::vector<pair<float,int>> get_mu_resonance(std::vector<pair<int,int>> cand_mu_pairs);
std::vector<pair<int,int>> get_W_cands(bool is_ee, std::vector<pair<int,int>> Z_cand_vector, std::vector<int> el_W_cands, std::vector<int> mu_W_cands);
