//#include "MuonSelections.h"
//#include "IsolationTools.h"
#include "Config.h"

using namespace tas;

std::vector<pair<int,int>> getCandPairs_el(std::vector<int> el_input) {
    //if (Muon_pt().at(idx) < 5.) { return false; }
    std::vector<pair<int,int>> el_pair;
    for ( int i = 0; i < el_input.size(); i++ ){
	  for ( int k = i+1; k < el_input.size(); k++ ){
	        if ( Electron_pdgId().at(el_input[i]) + Electron_pdgId().at(el_input[k]) != 0) continue;
		else{
		     float m_ee = (Electron_p4().at(el_input[i])+Electron_p4().at(el_input[k])).M();
		     
		     if ( m_ee < 12. ) continue;
		     el_pair.emplace_back(i,k);
		}
	  }	  
    }
    return el_pair;
}

std::vector<pair<int,int>> getCandPairs_mu(std::vector<int> mu_input) {
    std::vector<pair<int,int>> mu_pair;
    for ( int i = 0; i < mu_input.size(); i++ ){
          for ( int k = i+1; k < mu_input.size(); k++ ){
                if ( Muon_pdgId().at(mu_input[i]) + Muon_pdgId().at(mu_input[k]) != 0) continue;
                else{
                     float m_mm = (Muon_p4().at(mu_input[i])+Muon_p4().at(mu_input[k])).M();                         
		     if ( m_mm < 12. ) continue;
                     mu_pair.emplace_back(i,k);
                }
          }
    }
    return mu_pair;
}

std::vector<pair<float,int>> get_el_resonance(std::vector<pair<int,int>> cand_el_pairs){
    std::vector<pair<float,int>> el_mass_index;
    float m_Z = 91.2;
    float el_Z_cand_mass = 999999999999999.;
    for ( int j = 0; j < cand_el_pairs.size(); j++ ){
	  float m_ee = (Electron_p4().at(cand_el_pairs[j].first)+Electron_p4().at(cand_el_pairs[j].second)).M();
	  if ( std::abs(m_ee - m_Z) < std::abs(el_Z_cand_mass - m_Z ){
	       el_mass_index.clear();
               el_Z_cand_mass = m_ee;
	       el_mass_index.emplace_back(m_ee,j);
	  }
    }
    return el_mass_index;
}

std::vector<pair<float,int>> get_mu_resonance(std::vector<pair<int,int>> cand_mu_pairs){
    std::vector<pair<float,int>> mu_mass_index;
    float m_Z = 91.2;
    float mu_Z_cand_mass = 999999999999999.;
    for ( int j = 0; j < cand_mu_pairs.size(); j++ ){
          float m_mm = (Muon_p4().at(cand_mu_pairs[j].first)+Muon_p4().at(cand_mu_pairs[j].second)).M();
          if ( std::abs(m_mm - m_Z) < std::abs(mu_Z_cand_mass - m_Z ){
               mu_mass_index.clear();
               mu_Z_cand_mass = m_mm;
               mu_mass_index.emplace_back(m_mm,j);
          }
    }
    return mu_mass_index;
}

std::vector<pair<int,int>> get_W_cands(bool is_ee, std::vector<pair<int,int>> Z_cand_vector, std::vector<int> el_W_cands, std::vector<int> mu_W_cands){
   std::vector<pair<int,int>> W_cand;
   if ( is_ee ){
      if ( el_W_cands.size() > 0 && mu_W_cands.size() == 0 ){
           for ( int i = 0; i < Z_cand_vector.size(); i++ ){
		for ( int k = 0; k < el_W_cands.size(); k++ ){
		    if ( (el_W_cands[k] == Z_cand_vector[i].first) || (el_W_cands[k] == Z_cand_vector[i].second )) continue;
                    W_cand.emplace_back(el_W_cands[k],Electron_pdgId().at(el_W_cands[k])));	 	 
		    
		}
	   }
      }
      if ( el_W_cands.size() > 0 && mu_W_cands.size() > 0 ){
           for ( int i = 0; i < Z_cand_vector.size(); i++ ){
		for ( int k = 0; k < el_W_cands.size(); k++ ){
                    if ( (el_W_cands[k] == Z_cand_vector[i].first) || (el_W_cands[k] == Z_cand_vector[i].second )) continue;
                    W_cand.emplace_back(el_W_cands[k],Electron_pdgId().at(el_W_cands[k])));        
                }
           }
           for ( int k = 0; k < mu_W_cands.size(); k++ ){
                    W_cand.emplace_back(mu_W_cands[k],Muon_pdgId().at(mu_W_cands[k]))); 
                }
      }    
   }
   else{ 
     if ( mu_W_cands.size() > 0 && el_W_cands.size == 0 ){
	  for ( int p = 0; p < Z_cand_vector.size(); p++ ){
	        for ( int q = 0; q < mu_W_cands.size(); q++ ){
                    if ( (mu_W_cands[q] == Z_cand_vector[p].first) || (mu_W_cands[q] == Z_cand_vector[p].second) ) continue;
		    W_cand.emplace_back(mu_W_cands[q],Muon_pdgId().at(mu_W_cands[q])));
		}
	  }
     }
     if ( mu_W_cands.size() > 0 && el_W_cands.size() > 0 ){
	  for ( int p = 0; p < Z_cand_vector.size(); p++ ){
                for ( int q = 0; q < mu_W_cands.size(); q++ ){
                    if ( (mu_W_cands[q] == Z_cand_vector[p].first) || (mu_W_cands[q] == Z_cand_vector[p].second) ) continue;
                    W_cand.emplace_back(mu_W_cands[q],Muon_pdgId().at(mu_W_cands[q])));
                }
          }
          for ( int r = 0; r < el_W_cands.size(); r++ ){
		W_cand.emplace_back(el_W_cands[r],Electron_pdgId().at(el_W_cands[r]));
	  }
     }

   }

   return W_cand;
  
}

